using JuMP
using Ipopt

function get_mpcmodel(circuit, demand; has_ramping=true,
                      phase1=false, piecewise=false,
                      prev_val=nothing)

    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)

    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    if piecewise
        @variable(m, Cg[t=1:T,g=1:num_gens])
        @NLobjective(m, Min, sum(Cg[t,g] for t=1:T,g=1:num_gens))
        @constraint(m, plcurve[t=1:T,g=1:num_gens,p=1:gen[g].n-1],
		    Cg[t,g] - (((gen[g].coeff[2*p+2] - gen[g].coeff[2*p])/(gen[g].coeff[2*p+1] - gen[g].coeff[2*p-1]))*(baseMVA*Pg[t,g] - gen[g].coeff[2*p-1]) + gen[g].coeff[2*p]) >= 0
                   )
    else
        @NLobjective(m, Min, sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2
                                + gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g])
                                + gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens))
    end

    # Ramping up/down constraints
    if has_ramping
        if phase1 == false
	    @constraint(m, ramping[t=1:T-1,g=1:num_gens],
		        -gen[g].ramp_agc <= Pg[t+1,g] - Pg[t,g] <= gen[g].ramp_agc)
        else
	    @constraint(m, ramping[t=1:T,g=1:num_gens],
		        -gen[g].ramp_agc <= Pg[t,g] - prev_val[g] <= gen[g].ramp_agc)
        end
    end

    # Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
                  == 0)

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    # Line limits
    rateA = getfield.(line, :rateA)
    limind = findall((rateA .!= 0) .& (rateA .< 1.0e10))
    num_linelimits = length(limind)

    Yff_abs2 = zeros(num_linelimits)
    Yft_abs2 = zeros(num_linelimits)
    Yre = zeros(num_linelimits)
    Yim = zeros(num_linelimits)
    flowmax = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (from bus)
        l = limind[i]
        flowmax[i] = (line[l].rateA / baseMVA)^2
        Yff_abs2[i] = yline[l].YffR^2 + yline[l].YffI^2
        Yft_abs2[i] = yline[l].YftR^2 + yline[l].YftI^2
        Yre[i] = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
        Yim[i] = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR
    end

    @NLconstraint(m, flowmaxfrom[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].from]]^2 *
                  (Yff_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Yft_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    - Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <= 0)

    Ytf_abs2 = zeros(num_linelimits)
    Ytt_abs2 = zeros(num_linelimits)

    for i in 1:num_linelimits
        # Apparent power limits (to bus)
        l = limind[i]
        Ytf_abs2[i] = yline[l].YtfR^2 + yline[l].YtfI^2
        Ytt_abs2[i] = yline[l].YttR^2 + yline[l].YttI^2
        Yre[i] = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
        Yim[i] = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR
    end

    @NLconstraint(m, flowmaxto[t=1:T,i=1:num_linelimits],
                  Vm[t,busdict[line[limind[i]].to]]^2 *
                  (Ytf_abs2[i]*Vm[t,busdict[line[limind[i]].from]]^2
                   + Ytt_abs2[i]*Vm[t,busdict[line[limind[i]].to]]^2
                   + 2*Vm[t,busdict[line[limind[i]].from]]*Vm[t,busdict[line[limind[i]].to]]*
                   (Yre[i]*cos(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]])
                    -Yim[i]*sin(Va[t,busdict[line[limind[i]].from]] - Va[t,busdict[line[limind[i]].to]]))
                   ) - flowmax[i] <=0)
#=
    for l in 1:num_lines
        num_linelimits = 0

        if line[l].rateA != 0 && line[l].rateA < 1.0e10
            num_linelimits += 1
            flowmax = (line[l].rateA / baseMVA)^2

            # Apparent power limits (from bus)
            Yff_abs2 = yline[l].YffR^2 + yline[l].YffI^2
            Yft_abs2 = yline[l].YftR^2 + yline[l].YftI^2
            Yre = yline[l].YffR*yline[l].YftR + yline[l].YffI*yline[l].YftI
            Yim = -yline[l].YffR*yline[l].YftI + yline[l].YffI*yline[l].YftR

            @NLconstraint(m, [t=1:T],
                          Vm[t,busdict[line[l].from]]^2 *
                          (Yff_abs2*Vm[t,busdict[line[l].from]]^2
                           + Yft_abs2*Vm[t,busdict[line[l].to]]^2
                           + 2*Vm[t,busdict[line[l].from]]*Vm[t,busdict[line[l].to]]*
                           (Yre*cos(Va[t,busdict[line[l].from]] - Va[t,busdict[line[l].to]])
                            - Yim*sin(Va[t,busdict[line[l].from]] - Va[t,busdict[line[l].to]]))
                           ) - flowmax <= 0)

            # Apparent power limits (to bus)
            Ytf_abs2 = yline[l].YtfR^2 + yline[l].YtfI^2
            Ytt_abs2 = yline[l].YttR^2 + yline[l].YttI^2
            Yre = yline[l].YtfR*yline[l].YttR + yline[l].YtfI*yline[l].YttI
            Yim = -yline[l].YtfR*yline[l].YttI + yline[l].YtfI*yline[l].YttR

            @NLconstraint(m, [t=1:T],
                          Vm[t,busdict[line[l].to]]^2 *
                          (Ytf_abs2*Vm[t,busdict[line[l].from]]^2
                           + Ytt_abs2*Vm[t,busdict[line[l].to]]^2
                           + 2*Vm[t,busdict[line[l].from]]*Vm[t,busdict[line[l].to]]*
                           (Yre*cos(Va[t,busdict[line[l].from]] - Va[t,busdict[line[l].to]])
                            -Yim*sin(Va[t,busdict[line[l].from]] - Va[t,busdict[line[l].to]]))
                           ) - flowmax <=0)
        end
    end
=#

    return m
end

function get_mpcpfmodel(circuit, demand)
    m = Model()

    # Shortcuts
    baseMVA = circuit.baseMVA
    busref = circuit.busref
    bus = circuit.bus
    line = circuit.line
    gen = circuit.gen
    yline = circuit.yline
    ybus = circuit.ybus
    busdict = circuit.busdict
    frombus = circuit.frombus
    tobus = circuit.tobus
    bus2gen = circuit.bus2gen

    Pd = demand.pd
    Qd = demand.qd
    T = size(Pd,2)

    num_buses = length(bus)
    num_gens = length(gen)
    num_lines = length(line)

    @variable(m, gen[g].Pmin <= Pg[t=1:T,g=1:num_gens] <= gen[g].Pmax)
    @variable(m, gen[g].Qmin <= Qg[t=1:T,g=1:num_gens] <= gen[g].Qmax)
    @variable(m, bus[b].Vmin <= Vm[t=1:T,b=1:num_buses] <= bus[b].Vmax)
    @variable(m, Va[t=1:T,b=1:num_buses])

    for t in 1:T
        setlowerbound(Va[t,busref], bus[busref].Va)
        setupperbound(Va[t,busref], bus[busref].Va)
    end

    @objective(m, Min, 0)

    # Power flow constraints: real part
    @NLconstraint(m, pfreal[t=1:T,b=1:num_buses],
                  (sum(yline[l].YffR for l in frombus[b])
                   + sum(yline[l].YttR for l in tobus[b])
                   + ybus[b].YshR)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (yline[l].YftR*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftI*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (yline[l].YtfR*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfI*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Pg[t,g] for g in bus2gen[b]) - Pd[b,t]) / baseMVA
                  == 0)

    # Power flow constraints: imaginary part
    @NLconstraint(m, pfimag[t=1:T,b=1:num_buses],
                  (sum(-yline[l].YffI for l in frombus[b])
                   + sum(-yline[l].YttI for l in tobus[b])
                   - ybus[b].YshI)*Vm[t,b]^2
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].to]]*
                        (-yline[l].YftI*cos(Va[t,b]-Va[t,busdict[line[l].to]])
                         + yline[l].YftR*sin(Va[t,b]-Va[t,busdict[line[l].to]]))
                        for l in frombus[b])
                  + sum(Vm[t,b]*Vm[t,busdict[line[l].from]]*
                        (-yline[l].YtfI*cos(Va[t,b]-Va[t,busdict[line[l].from]])
                         + yline[l].YtfR*sin(Va[t,b]-Va[t,busdict[line[l].from]]))
                        for l in tobus[b])
                  - (sum(baseMVA*Qg[t,g] for g in bus2gen[b]) - Qd[b,t]) / baseMVA
                  == 0)

    return m
end

function get_mpcobjectivevalue(m, circuit, T)
    baseMVA = circuit.baseMVA
    gen = circuit.gen

    num_gens = length(gen)


    if gen[1].gentype == 1
	Cg = getvalue(m[:Cg])
	objval = sum(Cg[t,g] for t=1:T,g=1:num_gens)
    else
	Pg = getvalue(m[:Pg])
	objval = sum(gen[g].coeff[gen[g].n-2]*(baseMVA*Pg[t,g])^2 +
		     gen[g].coeff[gen[g].n-1]*(baseMVA*Pg[t,g]) +
		     gen[g].coeff[gen[g].n] for t=1:T,g=1:num_gens)
    end

    return objval
end

function get_mpcconstrviolation(m, d)
    viol = 0

    # Check constraints.
    if MathProgBase.numconstr(m) > 0
	g = zeros(MathProgBase.numconstr(m))
	MathProgBase.eval_g(d, g, m.colVal)

	g_lb, g_ub = JuMP.constraintbounds(m)
	for i=1:length(g_lb)
	    if g_lb[i] != -Inf
		err = max(0, -(g[i] - g_lb[i]))
		if viol < err
		    viol = err
		end
	    end

	    if g_ub[i] != Inf
		err = max(0, g[i] - g_ub[i])
		if viol < err
		    viol = err
		end
	    end
	end
    end

    # Check bound constraints.
    for i=1:MathProgBase.numvar(m)
	if m.colLower[i] != -Inf
	    err = max(0, -(m.colVal[i] - m.colLower[i]))
	    if viol < err
		viol = err
	    end
	end

	if m.colUpper[i] != Inf
	    err = max(0, m.colVal[i] - m.colUpper[i])
	    if viol < err
		viol = err
	    end
	end
    end

    return viol
end
