module MassSpecBioChemicals

using MassSpecChemicals
using MassSpecChemicals: AbstractChemical
using IterTools
import Base: show
import MassSpecChemicals: chemicalname, chemicalformula, chemicalabbr, repr_smiles
export FunctionalGroup, UnknownGroup, UnknownChemical, LeavingGroup, Dehydrogen, Odehydrogen, Ndehydrogen, Didehydrogen, Odidehydrogen, Ndidehydrogen, Deamine, Dehydroxy, Demethyl,
        Substituent, SubstitutedChemical, XLinkedFunctionalGroup, ChainedChemical, DehydratedChemical, IsotopiclabeledChemical,
        AbstractLinkageposition, Linkageposition, 
        originalmolecule, leavinggroup, conjugation, includeSIL, getchainlinkage, getchaincomponent, ischainedchemical

abstract type AbstractFunctionalGroup end
abstract type FunctionalGroup{M, S} <: AbstractFunctionalGroup end
abstract type UnknownGroup{M, S} <: AbstractFunctionalGroup end
struct UnknownChemical{F <: UnknownGroup} <: AbstractChemical end

abstract type LeavingGroup end
struct Dehydrogen <: LeavingGroup end # Dehydrogen
struct Odehydrogen <: LeavingGroup end # Dehydroxy
struct Ndehydrogen <: LeavingGroup end # Dehydroxy
struct Didehydrogen <: LeavingGroup end # Didehydrogen/Dehydrogen,Dehydrogen
# struct Odidehydrogen <: LeavingGroup end
# struct Ndidehydrogen <: LeavingGroup end
struct Dehydroxy <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen
struct Deamine <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen
struct Demethyl <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen
abstract type AbstractLinkageposition end
struct Linkageposition <: AbstractLinkageposition
    position::Union{Nothing, UInt8}
end
lk(x) = Linkageposition(UInt8(x))
lk(x::Nothing) = Linkageposition(x)

struct Substituent{M <: AbstractChemical, S <: LeavingGroup} <: FunctionalGroup{M, S}
    molecule::M
    position::AbstractLinkageposition
end
struct SubstitutedChemical{M <: AbstractChemical, S, T} <: AbstractChemical
    molecule::M
    substituted::S
    substituent::T
end 
#= 
S
LeavingGroup
Vector{Pair{LeavingGroup, UInt8}}: LeavingGroup => Number
Vector{Pair{AbstractLinkageposition, LeavingGroup}}: Position => LeavingGroup
T
Vector{Pair{FunctionalGroup, UInt8}}: FunctionalGroup => Number
Vector{Pair{AbstractLinkageposition, FunctionalGroup}}: Position => FunctionalGroup
=#
struct XLinkedFunctionalGroup{M, L} <: AbstractFunctionalGroup
    xlinkage::M
    functionalgroup::L
end
# isdissociated check adjacent C=O

struct ChainedChemical{M, L} <: AbstractChemical
    chain::M
    linkage::L
end
#=
L
Vector{Pair{LeavingGroup, LeavingGroup}}: LeavingGroup => LeavingGroup
Vector{Vector{Pair{LeavingGroup, LeavingGroup}}, Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}}: LeavingGroup => LeavingGroup Position([a/b]p for hex) => Position
=#
struct DehydratedChemical{M, L} <: AbstractChemical
    chain::M
    linkage::L
end
#=
L
Nothing
Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}: Position([a/b]p for hex) => Position
=#
struct IsotopiclabeledChemical{M, I} <: AbstractChemical
    molecule::M
    isotopiclabel::I
end
Substituent(::Type{S}, m::M, p = Linkageposition(0x00)) where {M, S} = Substituent{M, S}(m, p)
originalmolecule(fg::Substituent) = fg.molecule
originalmolecule(sm::SubstitutedChemical) = sm.molecule
originalmolecule(::FunctionalGroup{M}) where M = M()
originalmolecule(x::XLinkedFunctionalGroup) = originalmolecule(x.functionalgroup)
originalmolecule(::M) where {M <: UnknownGroup} = UnknownChemical{M}()
leavinggroup(::FunctionalGroup{M, S}) where {M, S} = S()
conjugation(m::AbstractChemical) = Substituent(Ndehydrogen, m)
includeSIL(T) = Union{<: T, <: IsotopiclabeledChemical{<: T}}

function makemolecule(m; sil)
    isnothing(sil) ? m : IsotopiclabeledChemical(m, s)
end

function makemolecule(::Type{T}, m::T; linkage = nothing, sil = nothing) where T
    isnothing(sil) ? m : IsotopiclabeledChemical(m, s)
end

function makemolecule(::Type{T}, ms...; linkage = nothing, sil = nothing) where T
    molecules = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makemolecule(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    length(last(molecules)) == 1 || throw(ArgumentError("Last chemical should be a single chemical"))
    ids = zeros(Int, length(molecules) - 1)
    i = lastindex(molecules) - 1
    j = lastindex(molecules)
    while i >= firstindex(molecules)
        ids[i] = j
        if length(molecules[i]) == 1
            j = i
            i -= 1
        else
            i -= 1
        end
    end
    lso = collect(Any, map(x -> collect(transformlinkage(T, x)), ms))
    ls = Any[makelinkage(T, last(getchaincomponent(a)), molecules[b]) for (a, b) in zip(molecules[begin:end - 1], ids)]
    for i in eachindex(ls)
        if !isnothing(linkage) && !isnothing(linkage[i])
            ls[i] = linkage[i]
        else 
            l = length(lso[i]) > 0 ? last(lso[i]) : nothing
            if isnothing(l) || isnulllinkage(T, l)
                continue
            elseif isnulllinkage(T, first(l))
                ls[i] = first(ls[i]) => last(l)
            elseif isnulllinkage(T, last(l))
                ls[i] = first(l) => last(ls[i])
            else
                ls[i] = l
            end
        end
    end
    T(molecules, ls)
end

function concatmolecule(m::T; linkage = nothing, sil) where T
    concatmolecule(T, m; linkage, sil)
end

function concatmolecule(::Type{T}, m::T; linkage = nothing, sil = nothing) where T
    isnothing(sil) && isnothing(linkage) && return m
    ischainedchemical(m) || return IsotopiclabeledChemical(m, s)
    ms = getchaincomponent(m)
    ls = getchainlinkage(m)
    if length(ls) < length(ms) - 1
        throw(ArgumentError("`linkage` for $(m) must have $(length(ms)) or $(length(ms) - 1) elements"))
    end
    molecules = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makemolecule(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    T(molecules, isnothing(linkage) ? ls : [first(ls, length(molecules) - 1)..., first(linkage)])
end

function concatmolecule(::Type{T}, ms...; linkage = nothing, sil = nothing) where T
    # T is chain ?
    mss = collect(Any, map(collect ∘ getchaincomponent, ms))
    lss = collect(Any, map(x -> collect(transformlinkage(T, x)), ms))
    lastlk = false
    if !isnothing(linkage)
        if length(linkage) == length(mss) - 1
            lastlk = false
        elseif length(linkage) == length(mss) && isnothing(last(linkage))
            lastlk = false
        elseif length(linkage) == length(mss)
            lastlk = true
        else
            throw(ArgumentError("`linkage` must have $(length(mss)) or $(length(mss) - 1) elements"))
        end
    end
    # linkage ?
    for i in eachindex(mss)
        if length(lss[i]) < length(mss[i]) - 1
            throw(ArgumentError("`linkage` for $(ms[i]) must have $(length(mss[i])) or $(length(mss[i]) - 1) elements"))
        elseif i == lastindex(mss)
            if lastlk
                lss[i] = [first(lss[i], length(mss[i]) - 1)..., last(linkage)]
            end
        elseif !isnothing(linkage) && !isnothing(linkage[i])
            lss[i] = [first(lss[i], length(mss[i]) - 1)..., linkage[i]]
        elseif length(lss[i]) >= length(mss[i])
            ls = lss[i][length(mss[i])]
            ni = findfirst(x -> length(x) == 1, mss[i + 1])
            nm = mss[i + 1][ni]
            if isnulllinkage(T, ls)
                lss[i] = [first(lss[i], length(mss[i]) - 1)..., makelinkage(T, last(mss[i]), nm)]
            elseif isnulllinkage(T, first(ls))
                lss[i] = [first(lss[i], length(mss[i]) - 1)..., first(makelinkage(T, last(mss[i]), nm)) => last(ls)]
            elseif isnulllinkage(T, last(ls))
                lss[i] = [first(lss[i], length(mss[i]) - 1)..., first(ls) => last(makelinkage(T, last(mss[i]), nm))]
            else
                lss[i] = first(lss[i], length(mss[i]))
            end
        else
            ni = findfirst(x -> length(x) == 1, mss[i + 1])
            nm = mss[i + 1][ni]
            lss[i] = [lss[i]..., makelinkage(T, last(mss[i]), first(mss[i + 1]))]
        end
    end
    ms = vcat(mss...)
    ls = vcat(lss...)
    molecules = if isnothing(sil)
        (ms..., )
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makemolecule(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    T(molecules, ls)
end

dehydrogenposition(a::AbstractChemical) = 0x00
dehydrogenposition(a::DehydratedChemical) = dehydrogenposition(first(getchaincomponent(a)))
dehydroxyposition(a::AbstractChemical) = 0x00
dehydroxyposition(a::DehydratedChemical) = dehydroxyposition(last(getchaincomponent(a)))
isnulllinkage(::Type{DehydratedChemical}, l) = l == lk(0x00)
isnulllinkage(::Type{DehydratedChemical}, l::Pair) = all(==(lk(0x00)), l)
isnulllinkage(::Type{ChainedChemical}, l) = first(l) == lk(0x00)
isnulllinkage(::Type{ChainedChemical}, l::Pair) = all(x -> ==(first(x), lk(0x00)), l)
makelinkage(::Type{ChainedChemical}, a, b) = (lk(dehydroxyposition(a)), Dehydroxy()) => (lk(dehydrogenposition(b)), Dehydrogen())
makelinkage(::Type{DehydratedChemical}, a, b) = lk(dehydroxyposition(a)) => lk(dehydrogenposition(b))
transformlinkage(::Type{<: AbstractChemical}, m::AbstractChemical) = getchainlinkage(m)
function transformlinkage(::Type{DehydratedChemical}, m::ChainedChemical)
    Tuple(first(first(ls)) => first(last(ls)) for ls in getchainlinkage(m))
end
function transformlinkage(::Type{ChainedChemical}, m::DehydratedChemical)
    Tuple((first(ls), Dehydroxy()) => (last(ls), Dehydrogen()) for ls in getchainlinkage(m))
end

getchaincomponent(m::AbstractChemical) = (m, )
getchaincomponent(m::ChainedChemical) = m.chain
getchaincomponent(m::DehydratedChemical) = m.chain
getchainlinkage(m::AbstractChemical) = ()
getchainlinkage(m::ChainedChemical) = m.linkage
getchainlinkage(m::DehydratedChemical) = m.linkage
ischainedchemical(m::AbstractChemical) = false
ischainedchemical(m::ChainedChemical) = true
ischainedchemical(m::DehydratedChemical) = true


include("io.jl")
include("BasicCompounds.jl")
include("Metabolites.jl")
include(joinpath("protein", "Proteins.jl"))
include(joinpath("glycan", "Glycans.jl"))
include(joinpath("nucleotide", "Nucleotides.jl"))
include(joinpath("lipid", "Lipids.jl"))

end