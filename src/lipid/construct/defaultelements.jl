carbonchaindefaultelements() = Pair{String, Int}[]

function carbonchaindefaultelements(chain::CarbonChain{T}) where T
    ec = chainelements(chain)
    es, lh = chemicalelements_nlinkage(chain.substituent)
    isnothing(es) && return ec
    ec[2] = "H" => (last(ec[2]) - lh)
    for e in es
        i = findall(x -> last(x) < 0, e)
        if !isempty(i)
            for j in i
                k = findfirst(x -> first(x) == first(e[j]), ec)
                ec[k] = first(ec[k]) => (last(ec[k]) + last(e[j]))
                filter!(x -> last(x) > 0, ec)
            end
            deleteat!(e, i)
        end
    end
    T <: Acyl ? vcat(ec[begin:end - 2], es..., ec[end - 1:end]) : vcat(ec, es...)
end

function chainelements(chain::CarbonChain{T}) where T
    ec = T <: Acyl ?
        ["C" => ncarbon(chain) - 1, "H" => (ncarbon(chain) * 2 + 2 - ndoublebond(chain) * 2 - 1), "C" => 1] :
        ["C" => ncarbon(chain), "H" => (ncarbon(chain) * 2 + 2 - ndoublebond(chain) * 2 - 1)]
    S = T <: Tuple ? T : Tuple{T}
    for C in Iterators.reverse(S.parameters)
        for (e, n) in chainelements(C)
            if n < 0
                i = findfirst(x -> first(x) == e, ec)
                ec[i] = e => (last(ec[i]) + n)
                filter!(x -> last(x) > 0, ec)
            else
                push!(ec, e => n)
            end
        end
    end
    ec
end

chainelements(::Type{<:Acyl}) = ["H" => -2, "O" => 1]
chainelements(::Type{<:Alkyl}) = Pair{String, Int}[]
chainelements(::Type{<:Alkenyl}) = ["H" => -2]
chainelements(::Type{<:SPB}) = ["N" => 1, "H" => 1]
chainelements(::Type{<:SulfoSPB}) = ["H" => -1, "S" => 1, "O" => 3, "H" => 1, "N" => 1, "H" => 1]

function chemicalelements_nlinkage(sub)
    if ispropertyempty(sub)
        (nothing, 0)
    elseif sub isa Number
        (["O" => Int(sub)], 0)
    elseif ispropertynumber(sub)
        sub = filter(x -> first(x) != Hydrogen(), sub)
        (map(sub) do x
            repeat(chemicalelements(first(x)), last(x))
        end, sum(x -> ntotallinkage(first(x)) * last(x), sub))
    else
        sub = filter(x -> last(x) != Hydrogen(), sub)
        sort!(sub; by = first)
        (chemicalelements.(last.(sub)), sum(ntotallinkage ∘ last, sub))
    end
end

function CarbonChain{T}(carbon::UInt8, doublebond, substituent, chirality, isotopiclabel = nothing;
        name = nothing, elements = carbonchaindefaultelements()) where {T}
    D, S, C, I = typeof(doublebond), typeof(substituent), typeof(chirality), typeof(isotopiclabel)
    chain = CarbonChain{T, D, S, C, I}(carbon, doublebond, substituent, chirality, isotopiclabel, "", Pair{String, Int}[])
    chainname = isnothing(name) ? repr_singlechain(chain) : String(name)
    CarbonChain{T, D, S, C, I}(carbon, doublebond, substituent, chirality, isotopiclabel,
        chainname, collect(Pair{String, Int}, elements))
end

_lipidchain_elements(chain::CarbonChain) = carbonchaindefaultelements(chain)
_lipidchain_elements(chain) = collect(Pair{String, Int}, chemicalelements(chain))
_lipidchain_elements(chains::Tuple) = reduce(vcat, _lipidchain_elements.(chains); init = Pair{String, Int}[])

function _adjust_backbone_hydrogens(elements::Vector{Pair{String, Int}}, count; fromlast = true)
    result = copy(elements)
    for _ in 1:count
        i = fromlast ? findlast(x -> first(x) == "H", result) : findfirst(x -> first(x) == "H", result)
        isnothing(i) && break
        result[i] = "H" => (last(result[i]) - 1)
        filter!(x -> last(x) > 0, result)
    end
    result
end
_adjust_backbone_hydrogens(backbone, count; kwargs...) =
    _adjust_backbone_hydrogens(chemicalelements(backbone), count; kwargs...)

function lipiddefaultelements(::Type{<:LipidChemical}, chemical)
    collect(Pair{String, Int}, chemicalelements(chemical))
end
function lipiddefaultelements(::Type{<:MonoFattyAcyl}, backbone, chain)
    b = collect(Pair{String, Int}, chemicalelements(backbone))
    if length(b) == 1
        elements = carbonchaindefaultelements(chain)
        i = findfirst(x -> first(x) == "H", elements)
        !isnothing(i) && (elements[i] = "H" => (last(elements[i]) + 1))
        return elements
    end
    vcat(_lipidchain_elements(chain), _adjust_backbone_hydrogens(b, 1; fromlast = false))
end
function lipiddefaultelements(::Type{<:NacylAmine}, backbone, chain)
    vcat(_lipidchain_elements(chain), _adjust_backbone_hydrogens(backbone, 1; fromlast = false))
end
function lipiddefaultelements(::Type{<:NacylAlkylAmine}, backbone, chain)
    b = reverse(chemicalelements(backbone))
    vcat(_lipidchain_elements(chain), _adjust_backbone_hydrogens(b, 1; fromlast = false))
end
function lipiddefaultelements(::Type{<:FattyAcylEster}, backbone, chain, position)
    b = reverse(chemicalelements(backbone))
    vcat(_lipidchain_elements(chain), _adjust_backbone_hydrogens(b, 1; fromlast = false))
end

lipiddefaultelements(::Type{<:Monoradylglycerol}, backbone, chain, args...) =
    vcat(_adjust_backbone_hydrogens(chemicalelements(backbone), 1), _lipidchain_elements(chain))
lipiddefaultelements(::Type{<:Diradylglycerol}, backbone, chain, args...) =
    vcat(_adjust_backbone_hydrogens(chemicalelements(backbone), 2), _lipidchain_elements(chain))
function lipiddefaultelements(::Type{<:Triradylglycerol}, backbone, chain, args...)
    b = _adjust_backbone_hydrogens(chemicalelements(backbone), 1; fromlast = false)
    vcat(_adjust_backbone_hydrogens(b, 2), _lipidchain_elements(chain))
end
lipiddefaultelements(::Type{<:Union{Omodifiedmonoradylglycerol, Monoradylglycerophosphate,
        Lysophosphatidylserine, LysophosphatidylNmodifiedserine, Lysophosphatidylglycerol,
        Lysophosphatidylglycerolphosphate}}, backbone, chain, args...) =
    vcat(_adjust_backbone_hydrogens(chemicalelements(backbone), 1), _lipidchain_elements(chain))
lipiddefaultelements(::Type{<:Union{Omodifieddiradylglycerol, Diradylglycerophosphate,
        Phosphatidylserine, PhosphatidylNmodifiedserine, Phosphatidylglycerol,
        Phosphatidylglycerolphosphate}}, backbone, chain, args...) =
    vcat(_adjust_backbone_hydrogens(chemicalelements(backbone), 2), _lipidchain_elements(chain))

function lipiddefaultelements(::Type{T}, backbone, chain, args...) where
        {T <: Union{Bisradylglycerophosphate, Bisradylglycerophosphoglycerol}}
    b = _adjust_backbone_hydrogens(chemicalelements(backbone), 1)
    b = _adjust_backbone_hydrogens(b, 1; fromlast = false)
    ncarbonchain(chain) > 2 && (b = _adjust_backbone_hydrogens(b, 1))
    ncarbonchain(chain) > 3 && (b = _adjust_backbone_hydrogens(b, 1; fromlast = false))
    vcat(b, _lipidchain_elements(chain))
end
lipiddefaultelements(::Type{<:GlycerophosphoNacylethanolamine}, backbone, chain, args...) =
    vcat(_lipidchain_elements(chain),
        _adjust_backbone_hydrogens(chemicalelements(backbone), 1; fromlast = false))

function _sphingolipid_chain_elements(chain)
    chains = tuplize(chain)
    elements = _lipidchain_elements(chains)
    firstchain = first(chains)
    chain_species = firstchain isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} ||
        (firstchain isa CarbonChain && !isnothing(firstchain.substituent) &&
         any(x -> first(x) isa OxygenAtom, firstchain.substituent))
    elements, firstchain, chain_species
end
function lipiddefaultelements(::Type{<:SphingoBone}, headgroup, chain, position, chirality)
    elements, firstchain, chain_species = _sphingolipid_chain_elements(chain)
    if firstchain isa CarbonChain{<:AbstractSPB}
        i = findfirst(x -> first(x) == "N", elements)
        elements[i + 1] = "H" => (last(elements[i + 1]) + 1)
    end
    if isnothing(headgroup) || !chain_species
        return elements
    end
    elements[2] = "H" => (last(elements[2]) - 1)
    vcat(elements, chemicalelements(Substituent(Dehydroxy, headgroup)))
end
function lipiddefaultelements(::Type{<:MixSphingoBone}, headgroup, chain, position, chirality)
    elements, firstchain, chain_species = _sphingolipid_chain_elements(chain)
    if chain isa CarbonChain{<:AbstractSPB}
        i = findfirst(x -> first(x) == "N", elements)
        elements[i + 1] = "H" => (last(elements[i + 1]) + 1)
    end
    if isnothing(headgroup) || !chain_species
        return elements
    end
    elements[2] = "H" => (last(elements[2]) - length(headgroup))
    vcat(elements, [chemicalelements(Substituent(Dehydroxy, x)) for x in headgroup]...)
end

lipiddefaultelements(::Type{<:SterolBone}, chain) = _lipidchain_elements(chain)
lipiddefaultelements(::Type{<:SubstitutedSterol}, backbone) = chemicalelements(backbone)
lipiddefaultelements(::Type{<:Retinylester}, chain) = _lipidchain_elements(chain)
lipiddefaultelements(::Type{<:CoenzymeQ}, chain) = _lipidchain_elements(chain)

lipiddefaultelements(::Type{T}, backbone, chain, args...) where {T <: Lipid} =
    vcat(_lipidchain_elements(chain), collect(Pair{String, Int}, chemicalelements(backbone)))
lipiddefaultelements(::Type{T}, chain) where {T <: Lipid} = _lipidchain_elements(chain)
