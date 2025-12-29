parentchemical(m::AbstractChemical) = m
parentchemical(sil::IsotopiclabeledChemical) = sil.chemical
parentchemical(fg::Substituent) = parentchemical(fg.chemical)
parentchemical(sm::SubstitutedChemical) = sm.chemical
parentchemical(::FunctionalGroup{M}) where M = M()
parentchemical(x::XLinkedFunctionalGroup) = parentchemical(x.functionalgroup)
parentchemical(::M) where {M <: UnknownGroup} = UnknownChemical{M}()

leavinggroup(::FunctionalGroup{M, S}) where {M, S} = S()
function dehydrogengroup(m; position = nothing) 
    if isnothing(position) || ismissing(position)
        position = dehydrogenposition(m)
    end
    if isnothing(position) || ismissing(position)
        position = Linkageposition(0x00)
    elseif position isa Number 
        position = Linkageposition(position)
    end
    Substituent(Dehydrogen, m, position)
end
function dehydroxygroup(m; position = nothing) 
    if isnothing(position) || ismissing(position)
        position = dehydroxyposition(m)
    end
    if isnothing(position) || ismissing(position)
        position = Linkageposition(0x00)
    elseif position isa Number 
        position = Linkageposition(position)
    end
    Substituent(Dehydroxy, m, position)
end
conjugation(m::AbstractChemical) = try
        dehydrogengroup(m) 
    catch
        dehydroxygroup(m)
    end

includeSIL(T) = Union{<: T, <: IsotopiclabeledChemical{<: T}}

leavinggroupelements(::Dehydrogen) = ["H" => -1]
leavinggroupelements(::Odehydrogen) = ["H" => -1]
leavinggroupelements(::Ndehydrogen) = ["H" => -1]
leavinggroupelements(::Didehydrogen) = ["H" => -2]
leavinggroupelements(::Dehydroxy) = ["O" => -1, "H" => -1]
leavinggroupelements(::Deamine) = ["N" => -1, "H" => -2]
leavinggroupelements(::Demethine) = ["C" => -1, "H" => -1]

function losswaterelements(m::AbstractChemical, p, q)
    e = chemicalelements(m)
    if p != 1 
        i = findfirst(x -> first(x) == "H", e)
        e[i] = "H" => (last(e[i]) - 1)
    end
    if q != 0
        i = findlast(x -> first(x) == "O", e)
        e[i] = "O" => (last(e[i]) - 1)
        i = findlast(x -> first(x) == "H", e)
        e[i] = "H" => (last(e[i]) - 1)
    end
    filter!(x -> last(x) != 0, e)
end

dehydrogenposition(a::AbstractChemical) = 0x00
dehydrogenposition(a::DehydratedChemical) = dehydrogenposition(first(getchaincomponent(a)))
dehydroxyposition(a::AbstractChemical) = 0x00
dehydroxyposition(a::DehydratedChemical) = dehydroxyposition(last(getchaincomponent(a)))

dehydrogenfunctionalgroup(a::AbstractChemical; position = nothing) = Hydroxy()
dehydrogenfunctionalgroup(a::DehydratedChemical; position = nothing) = dehydrogenfunctionalgroup(first(getchaincomponent(a)))
dehydrogenfunctionalgroup(a::Substituent) = dehydrogenfunctionalgroup(a.chemical; position = a.position)
dehydrogenfunctionalgroup(a::FunctionalGroup{M, <: Union{Dehydrogen, Odehydrogen, Ndehydrogen}}) where M = dehydrogenfunctionalgroup(parentchemical(a); position = dehydrogenposition(parentchemical(a)))
dehydroxyfunctionalgroup(a::AbstractChemical; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(a::DehydratedChemical; position = nothing) = dehydroxyfunctionalgroup(first(getchaincomponent(a)))
dehydroxyfunctionalgroup(a::Substituent) = dehydroxyfunctionalgroup(a.chemical; position = a.position)
dehydroxyfunctionalgroup(a::FunctionalGroup{M, Dehydroxy}) where M = dehydroxyfunctionalgroup(parentchemical(a); position = dehydroxyposition(parentchemical(a)))

isnulllinkage(::Type{<: DehydratedChemical}, l) = l == lk(0x00)
isnulllinkage(::Type{<: DehydratedChemical}, l::Pair) = all(==(lk(0x00)), l)
isnulllinkage(::Type{<: ChainedChemical}, l) = first(l) == lk(0x00)
isnulllinkage(::Type{<: ChainedChemical}, l::Pair) = all(x -> ==(first(x), lk(0x00)), l)

makelinkage(::Type{<: ChainedChemical}, a, b) = (lk(dehydroxyposition(a)), Dehydroxy()) => (lk(dehydrogenposition(b)), Dehydrogen())
makelinkage(::Type{<: DehydratedChemical}, a, b) = lk(dehydroxyposition(a)) => lk(dehydrogenposition(b))

transformlinkage(::Type{<: AbstractChemical}, m::AbstractChemical) = getchainlinkage(m)
function transformlinkage(::Type{<: DehydratedChemical}, m::ChainedChemical)
    Tuple(first(first(ls)) => first(last(ls)) for ls in getchainlinkage(m))
end
function transformlinkage(::Type{<: ChainedChemical}, m::DehydratedChemical)
    Tuple((first(ls), Dehydroxy()) => (last(ls), Dehydrogen()) for ls in getchainlinkage(m))
end

nlinkage(::AbstractFunctionalGroup) = 1
nlinkage(::FunctionalGroup) = 1
nlinkage(::FunctionalGroup{T, Didehydrogen}) where T = 2
nlinkage(::XLinkedFunctionalGroup) = 1
nbridge(::AbstractFunctionalGroup) = 0
ntotallinkage(x) = nlinkage(x) + nbridge(x)

composition(m::AbstractChemical) = Dict{AbstractChemical, UInt8}(m => 0x01)
function composition(m::ChainedChemical) 
    d = Dict{AbstractChemical, UInt8}()
    for c in getchaincomponent(m)
        # iterate composition
        get!(d, c, 0x00)
        d[c] += 0x01 
    end
    d 
end
function composition(m::DehydratedChemical) 
    d = Dict{AbstractChemical, UInt8}()
    for c in getchaincomponent(m)
        # iterate composition
        get!(d, c, 0x00)
        d[c] += 0x01 
    end
    d 
end
getchaincomponent(m::AbstractChemical) = (m, )
getchaincomponent(m::ChainedChemical) = m.chain
getchaincomponent(m::DehydratedChemical) = m.chain
getchainlinkage(m::AbstractChemical) = ()
getchainlinkage(m::ChainedChemical) = m.linkage
getchainlinkage(m::DehydratedChemical) = m.linkage
getchainconfig(m::AbstractChemical) = nothing 
getchainconfig(m::ChainedChemical) = m.config
getchainconfig(m::DehydratedChemical) = m.config
ischainedchemical(m::Type{<: AbstractChemical}) = false
ischainedchemical(m::Type{<: ChainedChemical}) = true
ischainedchemical(m::Type{<: DehydratedChemical}) = true
ischainedchemical(m::AbstractChemical) = false
ischainedchemical(m::ChainedChemical) = true
ischainedchemical(m::DehydratedChemical) = true
requirelinkage(m::Type{<: AbstractChemical}) = false
requirelinkage(m::Type{<: ChainedChemical}) = true
requirelinkage(m::Type{<: DehydratedChemical}) = true
requireconfig(m::Type{<: AbstractChemical}) = false
requireconfig(m::Type{<: ChainedChemical}) = true
requireconfig(m::Type{<: DehydratedChemical}) = true

isdissociated(fg) = false
deisomerize(m::AbstractChemical) = m
deisomerize(m::DehydratedChemical) = DehydratedChemical(map(deisomerize, m), nothing, nothing)
deisomerize(m::ChainedChemical) = ChainedChemical(map(deisomerize, m), nothing, nothing)
decompose(m::AbstractChemical) = Dict(m => 0x01)
function decompose(m::DehydratedChemical)
    d = Dict{AbstractChemical, UInt8}(Hydroxy() => 0x01)
    p = map(decompose, getchaincomponent(m))
    prev = nothing
    for ((i, a), c) in zip(enumerate(getchaincomponent(m)), p)
        for (k, v) in c 
            get!(d, k, 0x00)
            d[k] += v 
        end
        if isnothing(m.linkage) || isnothing(prev)
            d[Hydroxy()] -= 0x01
        elseif prev == CarboxylicAcidGroup()
            curr = dehydrogenfunctionalgroup(a; position = first(m.linkage[i]))
            d[curr] -= 0x01
            d[CarboxylicAcidGroup()] -= 0x01
            if curr == Hydroxy()
                d[Ester()] += 0x01
            elseif curr == Amino() || curr == MethylAmino() 
                d[Amide()] += 0x01
            elseif curr == Sulfo()
                d[Thioester()] += 0x01
            else
                d[Oxo()] += 0x01
            end
        else
            d[Hydroxy()] -= 0x01
        end
        prev = dehydroxyfunctionalgroup(a; position = first(m.linkage[i]))
    end
    d
end

function decompose(fg::Substituent)
    d = decompose(df.chemical)
    if S == Dehydroxy 
        d[Hydroxy()] -= 0x01 
    end
    d
end
function decompose(fg::XLinkedFunctionalGroup)
    d = decompose(fg.functionalgroup)
    get!(d, fg.linkage, 0x00)
    d[fg.linkage] += 0x01
    d
end
function decompose(fg::XLinkedFunctionalGroup{CarboxylicLinkage}) 
    d = decompose(fg.functionalgroup)
    curr = dehydrogenfunctionalgroup(fg.functionalgroup)
    if curr == Hydroxy()
        d[Ester()] += 0x01
    elseif curr == Amino() || curr == MethylAmino() 
        d[Amide()] += 0x01
    elseif curr == Sulfo()
        d[Thioester()] += 0x01
    else
        d[Oxo()] += 0x01
    end
    d
end
