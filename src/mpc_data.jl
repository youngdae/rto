using Printf

mutable struct Bus
    bus_i::Int
    bustype::Int
    Pd::Float64
    Qd::Float64
    Gs::Float64
    Bs::Float64
    area::Int
    Vm::Float64
    Va::Float64
    baseKV::Float64
    zone::Int
    Vmax::Float64
    Vmin::Float64
end

mutable struct Line
    from::Int
    to::Int
    r::Float64
    x::Float64
    b::Float64
    rateA::Float64
    rateB::Float64
    rateC::Float64
    ratio::Float64 #TAP
    angle::Float64 #SHIFT
    status::Int
    angmin::Float64
    angmax::Float64
end

mutable struct Gen
    # .gen fields
    bus::Int
    Pg::Float64
    Qg::Float64
    Qmax::Float64
    Qmin::Float64
    Vg::Float64
    mBase::Float64
    status::Int
    Pmax::Float64
    Pmin::Float64
    Pc1::Float64
    Pc2::Float64
    Qc1min::Float64
    Qc1max::Float64
    Qc2min::Float64
    Qc2max::Float64
    ramp_agc::Float64
    # .gencost fields
    gentype::Int
    startup::Float64
    shutdown::Float64
    n::Int
    coeff::Array{Float64}
end

struct Yline
    from::Int
    to::Int
    YffR::Float64
    YffI::Float64
    YttR::Float64
    YttI::Float64
    YtfR::Float64
    YtfI::Float64
    YftR::Float64
    YftI::Float64
end

struct Ybus
    bus::Int
    YshR::Float64
    YshI::Float64
end

struct Circuit
    baseMVA::Float64
    busref::Int
    bus::Array{Bus}
    line::Array{Line}
    gen::Array{Gen}
    yline::Array{Yline}
    ybus::Array{Ybus}
    busdict::Dict{Int,Int}
    frombus::Array
    tobus::Array
    bus2gen::Array
end

struct Load
    pd
    qd
end

function get_busmap(bus)
    busdict = Dict{Int,Int}()

    for i in 1:length(bus)
        @assert !haskey(busdict,bus[i].bus_i)
        busdict[bus[i].bus_i] = i
    end

    return busdict
end

function get_linetobusmap(bus, line, busdict)
    num_buses = length(bus)
    from = [Int[] for i in 1:num_buses]
    to = [Int[] for i in 1:num_buses]

    for i in 1:length(line)
        idx = busdict[line[i].from]
        @assert 1 <= idx <= num_buses
        push!(from[idx], i)

        idx = busdict[line[i].to]
        @assert 1 <= idx <= num_buses
        push!(to[idx], i)
    end

    return from, to
end

function get_bustogenmap(bus, gen, busdict)
    bus2gen = [Int[] for i in 1:length(bus)]

    for i in 1:length(gen)
        idx = busdict[gen[i].bus]
        push!(bus2gen[idx], i)
    end

    return bus2gen
end

# -------------------------------------------------------------------------
# Compute admittances.
# -------------------------------------------------------------------------
function getY(case, line, bus, baseMVA)
    dim1 = size(line,1)
    Ys = complex(zeros(dim1, 1))
    tap = complex(ones(dim1, 1))
    Ytt = complex(zeros(dim1, 1))
    Yff = complex(zeros(dim1, 1))
    Yft = complex(zeros(dim1, 1))
    Ytf = complex(zeros(dim1, 1))

    # ---------------------------------------------------------------------
    # bus f: tap bus, bus t: impedance bus or Z bus
    #
    # Ys: the admittance between bus f and bus t.
    #     It is the reciprocal of a series impedance between bus f and t.
    #
    # When there is a off-nominal transformer, the admittance matrix is
    # defined as follows:
    #
    #   / Ift \ = / Yff  Yft \ / Vf \
    #   \ Itf / = \ Ytf  Ytt / \ Vt /
    #
    # where
    #
    #    Yff = ((Ys + j*bft/2) / |a|^2)
    #    Yft = (-Ys / conj(a))
    #    Ytf = (-Ys / a)
    #    Ytt = (Ys + j*bft/2)
    #
    # When we compute If or It (total current injection at bus f or t),
    # we need to add the bus shunt, YshR and YshI.
    # ---------------------------------------------------------------------

    Ys = 1 ./ (line[:,3] .+ line[:,4] .* im)  # r
    ind = findall(line[:,9] .!= 0)
    tap[ind] = line[ind,9]                    # ratio
    tap .*= exp.(line[:,10] .* pi/180 .* im)  # angle
    Ytt = Ys .+ line[:,5] ./ 2 .* im          # b
    Yff = Ytt ./ (tap.*conj.(tap))
    Yft = -Ys ./ conj.(tap)
    Ytf = -Ys ./ tap

    yline = zeros(dim1, 10)
    yline[:,1:2] = line[:,1:2]
    yline[:,3] = real.(Yff); yline[:,4] = imag.(Yff)
    yline[:,5] = real.(Ytt); yline[:,6] = imag.(Ytt)
    yline[:,7] = real.(Ytf); yline[:,8] = imag.(Ytf)
    yline[:,9] = real.(Yft); yline[:,10] = imag.(Yft)

    ybus = zeros(size(bus,1), 3)
    YshR = bus[:,5] ./ baseMVA      # Gs
    YshI = bus[:,6] ./ baseMVA      # Bs
    ybus = [ bus[:,1] YshR YshI ]

    @assert 0==length(findall(isnan.(yline[:,3:10])))
    @assert 0==length(findall(isinf.(yline[:,3:10])))
    @assert 0==length(findall(isnan.(ybus[:,2:3])))
    @assert 0==length(findall(isinf.(ybus[:,2:3])))

    nylines = size(yline,1)
    nybuses = size(ybus,1)
    Ylines = Array{Yline}(undef, nylines)
    Ybuses = Array{Ybus}(undef, nybuses)

    for i in 1:nylines
        Ylines[i] = Yline(yline[i,1:end]...)
    end

    for i in 1:nybuses
        Ybuses[i] = Ybus(ybus[i,1:end]...)
    end

    return Ylines, Ybuses
end

# -------------------------------------------------------------------------
# Get circuit and power demand data for OPF computation.
# -------------------------------------------------------------------------
function getcircuit(case, baseMVA, ramp_scaling; neg_line=-1, neg_gen=-1)
    bus_mat = readdlm(case*".bus")
    bus_mat[:,9] *= pi/180  # multiply Va by pi/180.

    branch_mat = readdlm(case*".branch")
    active_line_ind = findall(branch_mat[:,11] .> 0)

    if neg_line != -1
        l = active_line_ind[neg_line]
        @printf("The line (%5d,%5d) has been cut.\n",
                 branch_mat[l,1], branch_mat[l,2])
        deleteat!(active_line_ind, neg_line)
    end

    line_mat = branch_mat[active_line_ind,:]

    gen_mat = readdlm(case*".gen")
    gencost_mat = readdlm(case*".gencost")

    # Allocate space for coefficients of a quadratic objective function.
    gens_on = findall(x-> x != 0, gen_mat[:,8])
    if neg_gen != -1
        @printf("Generator at %d connected to bus %d has been turned off.\n",
                 gens_on[neg_gen], gen_mat[gens_on[neg_gen],1])
        deleteat!(gens_on, neg_gen)
    end
    num_on = length(gens_on)

    # Compute admittances.
    yline, ybus = getY(case, line_mat, bus_mat, baseMVA)

    num_buses = size(bus_mat,1)
    num_lines = size(line_mat,1)
    num_gens = size(gen_mat,1)

    bus = Array{Bus}(undef, num_buses)
    line = Array{Line}(undef, num_lines)
    gen = Array{Gen}(undef, num_on)

    busref = -1
    for i in 1:num_buses
        bus[i] = Bus(bus_mat[i,1:end]...)
        if bus[i].bustype == 3
            if busref > 0
                error("More than one reference bus present")
            else
                busref = i
            end
        end
    end

    for i in 1:num_lines
        line[i] = Line(line_mat[i,1:end]...)
    end

    j = 0
    for i in gens_on
        j += 1

        gen[j] = Gen(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Array{Int}(undef,0))
        gen[j].bus = gen_mat[i,1]
        gen[j].Pg = gen_mat[i,2] / baseMVA
        gen[j].Qg = gen_mat[i,3] / baseMVA
        gen[j].Qmax = gen_mat[i,4] / baseMVA
        gen[j].Qmin = gen_mat[i,5] / baseMVA
        gen[j].Vg = gen_mat[i,6]
        gen[j].mBase = gen_mat[i,7]
        gen[j].status = gen_mat[i,8]
        @assert gen[j].status == 1
        gen[j].Pmax = gen_mat[i,9] / baseMVA
        gen[j].Pmin = gen_mat[i,10] / baseMVA
        gen[j].Pc1 = gen_mat[i,11]
        gen[j].Pc2 = gen_mat[i,12]
        gen[j].Qc1min = gen_mat[i,13]
        gen[j].Qc1max = gen_mat[i,14]
        gen[j].Qc2min = gen_mat[i,15]
        gen[j].Qc2max = gen_mat[i,16]
        if gencost_mat[i,1] == 1
            gen[j].ramp_agc = gen_mat[i,17] / baseMVA
        else
            gen[j].ramp_agc = gen[j].Pmax * ramp_scaling
        end
        gen[j].gentype = gencost_mat[i,1]
        gen[j].startup = gencost_mat[i,2]
        gen[j].shutdown = gencost_mat[i,3]
        gen[j].n = gencost_mat[i,4]
        gen[j].coeff = gencost_mat[i,5:end]
    end

    # Create dictionaries due to the lack of set data structure in JuMP.
    busdict = get_busmap(bus)
    frombus, tobus = get_linetobusmap(bus, line, busdict)
    bus2gen = get_bustogenmap(bus, gen, busdict)

    circuit = Circuit(baseMVA, busref, bus, line, gen, yline, ybus,
                      busdict, frombus, tobus, bus2gen)

    return circuit
end

function getload(scen, load_scale)
    pd_mat = readdlm(scen*".Pd")
    qd_mat = readdlm(scen*".Qd")

    load = Load(pd_mat.*load_scale, qd_mat.*load_scale)

    return load
end
