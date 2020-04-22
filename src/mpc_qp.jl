using Test
using JuMP
using MathProgBase
using SparseArrays
using Printf

function mpc_check_KKT(m::Model, F::Array{Float64})

    # --------------------------------------------------------------------
    # Compute the KKT error using the normal map equation:
    #
    #   1. Set y = z - F(z) --> y should be in the normal cone to B at z,
    #                           where B is a box constraint.
    #   2. if z = lb, error += max(0, y)   --> y should be nonpositive.
    #      if z = ub, error += max(0, -y)  --> y should be nonnegative.
    #      if lb < z < ub, error += abs(y) --> y should be zero.
    # --------------------------------------------------------------------

    prim_err = 0
    dual_err= 0

    lb = m.colLower
    ub = m.colUpper
    x = internalmodel(m).inner.x

    for i=1:length(lb)
	if x[i] == lb[i]
	    dual_err += max(0, -F[i])
	elseif x[i] == ub[i]
	    dual_err += max(0, F[i])
	else
	    dual_err += abs(F[i])
	end
    end

    if MathProgBase.numconstr(m) > 0
	n = MathProgBase.numvar(m)
	g_lb, g_ub = JuMP.constraintbounds(m)
	μ = internalmodel(m).inner.mult_g

	for i=1:length(g_lb)
	    if g_lb[i] == g_ub[i] # g should be satisfied as equality.
		prim_err += abs((F[n+i] - g_lb[i]))
	    elseif g_lb[i] != -Inf && μ[i] <= -1e-6 # g should be at its lower bound.
		prim_err += abs(F[n+i] - g_lb[i])
	    elseif g_ub[i] != Inf && μ[i] >= 1e-6   # g should be at its upper bound.
		prim_err += abs(F[n+i] - g_ub[i])
	    else # g should be in between lower/upper bounds.
		if g_lb[i] != -Inf
		    prim_err += max(0, -(F[n+i] - g_lb[i]))
		end

		if g_ub[i] != Inf
		    prim_err += max(0, F[n+i] - g_ub[i])
		end
	    end
	end
    end

    return prim_err, dual_err
end

function mpc_qp(m::Model; optfile="ipopt.op2", check_KKT=false)
    nvar = MathProgBase.numvar(m)
    nconstr = MathProgBase.numconstr(m)

    inner = internalmodel(m).inner

    println("Instantiating NLPEvaluator...")
    flush(stdout)

    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:Grad, :Jac, :Hess])

    println("Evaluating gradient, Jacobian, and Hessian...")
    flush(stdout)

    # --------------------------------------------------------------------
    # Evaluate Hessian: σ∇^2f(x) + ∑ μi*∇^2gi(x) with σ = 1.0.
    # We make it upper triangular and merge duplicates.
    # --------------------------------------------------------------------

    I, J = MathProgBase.hesslag_structure(d)
    Hv_tmp = zeros(length(I))
    MathProgBase.eval_hesslag(d, Hv_tmp, inner.x, 1.0, inner.mult_g)

    # Make it upper-triangular.
    for i=1:length(I)
	if I[i] > J[i]
	    I[i], J[i] = J[i], I[i]
	end
    end

    # Merge duplicates.
    Ih, Jh, Vh = findnz(sparse(I, J, [Int[i] for i=1:length(I)],
			       nvar, nvar, vcat))
    Hv = zeros(length(Ih))
    for i=1:length(Ih)
	Hv[i] = sum(Hv_tmp[Vh[i]])
    end

    # --------------------------------------------------------------------
    # Evaluate Jacobian: J = [∇gi(x)] for i=1..m.
    # We do not need to merge duplicates. They'll be added in the computation.
    # --------------------------------------------------------------------

    if nconstr > 0
	Ij, Jj = MathProgBase.jac_structure(d)
	Jv = zeros(length(Ij))
	MathProgBase.eval_jac_g(d, Jv, inner.x)
    end

    # --------------------------------------------------------------------
    # Evaluate gradient of f: ∇f(x).
    # --------------------------------------------------------------------

    gradf = zeros(length(inner.x))
    MathProgBase.eval_grad_f(d, gradf, inner.x)

    # --------------------------------------------------------------------
    # Evaluate constraints: g(x).
    # --------------------------------------------------------------------

    if nconstr > 0
	g = zeros(nconstr)
	MathProgBase.eval_g(d, g, inner.x)
    end

    if check_KKT
	# --------------------------------------------------------------------
	# Evaluate Jacobian-vector product: [μi*∇gi(x)] for i=1..m.
	# --------------------------------------------------------------------

	if nconstr > 0
	    jw = zeros(nvar)

	    for i=1:length(Ij)
		jw[Jj[i]] += Jv[i]*inner.mult_g[Ij[i]]
	    end
	end

	if nconstr > 0
	    c = [ gradf .+ jw; g ]
	else
	    c = gradf
	end

	# Check if the current solution satisfies the KKT conditions.
	prim_err, dual_err = mpc_check_KKT(m, c)
	@test prim_err <= 1e-6
	@test dual_err <= 1e-6
    end

    # ---------------------------------------------------------------------
    # Construct a QP:
    #
    #    min 0.5*(x-x0)*H*(x-x0) + ∇f(x0)(x-x0) + f(x0)
    #    s.t. g(x0) + ∇g(x0)(x-x0) ∈ K
    #
    # We eliminate the constants from the objective:
    #
    # => min 0.5*xHx - x0Hx + ∇f(x0)x
    #    s.t. g(x0) + ∇g(x0)(x-x0) ∈ K
    # ---------------------------------------------------------------------

    println("Building a QP model...")
    flush(stdout)

    m_qp = Model()

    lb = m.colLower; ub = m.colUpper
    g_lb, g_ub = JuMP.constraintbounds(m)

    x0 = inner.x
    @variable(m_qp, lb[i] <= x[i=1:nvar] <= ub[i], start=x0[i])

    Jj_vec = [Int[] for i=1:nconstr]
    Jv_vec = [Float64[] for i=1:nconstr]
    for i=1:length(Ij)
	push!(Jj_vec[Ij[i]], Jj[i])
	push!(Jv_vec[Ij[i]], Jv[i])
    end

    @constraint(m_qp, linconstr[i=1:nconstr],
		g_lb[i] <= g[i] + sum(Jv_vec[i][j]*(x[Jj_vec[i][j]] - x0[Jj_vec[i][j]]) for j=1:length(Jj_vec[i])) <= g_ub[i])

    D_Hv = zeros(nvar)
    off_Ih = Int[]
    off_Jh = Int[]
    off_Hv = Float64[]
    for i=1:length(Ih)
	if Ih[i] == Jh[i]
	    D_Hv[Ih[i]] = Hv[i]
	else
	    push!(off_Ih, Ih[i])
	    push!(off_Jh, Jh[i])
	    push!(off_Hv, Hv[i])
	end
    end

    @NLobjective(m_qp, Min,
		 0.5*(sum(D_Hv[i]*x[i]*x[i] for i=1:nvar)
		      + 2*sum(off_Hv[i]*x[off_Ih[i]]*x[off_Jh[i]] for i=1:length(off_Ih))
		      )
		 - (sum(D_Hv[i]*x[i]*x0[i] for i=1:nvar)
		    + sum(off_Hv[i]*x[off_Ih[i]]*x0[off_Jh[i]] for i=1:length(off_Ih))
		    + sum(off_Hv[i]*x[off_Jh[i]]*x0[off_Ih[i]] for i=1:length(off_Ih))
		    )
		 + sum(gradf[i]*x[i] for i=1:nvar)
		 )

    setsolver(m_qp, IpoptSolver(option_file_name=optfile))
    JuMP.build(m_qp)
    inner_m_qp = internalmodel(m_qp).inner
    copyto!(inner_m_qp.x, 1, inner.x, 1, nvar)
    copyto!(inner_m_qp.mult_x_L, 1, inner.mult_x_L, 1, nvar)
    copyto!(inner_m_qp.mult_x_U, 1, inner.mult_x_U, 1, nvar)

    if nconstr > 0
	copyto!(inner_m_qp.mult_g, 1, inner.mult_g, 1, nconstr)
    end

    MathProgBase.setwarmstart!(internalmodel(m_qp), inner_m_qp.x)
    MathProgBase.optimize!(internalmodel(m_qp))

    status = MathProgBase.status(internalmodel(m_qp))
    if status != :Infeasible && status != :Unbounded
    else
	error("QP failed: ", status)
    end

    return d, m_qp
end
