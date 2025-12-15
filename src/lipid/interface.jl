# CarbonChain
# MassSpecChemicals
function getchemicalattr(chain::CarbonChain{T}, ::Val{:elements}; kwargs...) where T
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
    if T <: Acyl 
        vcat(ec[begin:end - 2], es..., ec[end - 1:end])
    else
        vcat(ec, es...)
    end
end

function chainelements(chain::CarbonChain{T}) where T 
    ec = T <: Acyl ? ["C" => ncarbon(chain) - 1, "H" => (ncarbon(chain) * 2 + 2 - ndoublebond(chain) * 2 - 1), "C" => 1] : ["C" => ncarbon(chain), "H" => (ncarbon(chain) * 2 + 2 - ndoublebond(chain) * 2 - 1)]
    S = T <: Tuple ? T : Tuple{T}
    for C in Iterators.reverse(S.parameters)
        for (e, n) in chainelements(C)
            if n < 0
                i = findfirst(x -> ==(first(x), e), ec)
                ec[i] = e => (last(ec[i]) + n)
                filter!(x -> last(x) > 0, ec)
            else
                push!(ec, e => n)
            end
        end
    end
    ec
end

chainelements(::Type{Acyl}) = ["H" => - 2, "O" => 1]
chainelements(::Type{Alkyl}) = []
chainelements(::Type{Alkenyl}) = ["H" => - 2]
chainelements(::Type{SPB}) = ["N" => 1, "H" => 1]
chainelements(::Type{SulfoSPB}) = ["H" => -1, "S" => 1, "O" => 3, "H" => 1, "N" => 1, "H" => 1]

function chemicalelements_nlinkage(sub)
    if ispropertyempty(sub) 
        (nothing, 0)
    elseif isa(sub, Number) 
        return (["O" => Int(sub)], 0)
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

function getchemicalattr(lipid::MonoFattyAcyl, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    if length(b) == 1
        ec = chemicalelements(lipid.chain)
        i = findfirst(x -> ==(first(x), "H"), ec)
        ec[i] = "H" => (last(ec[i]) + 1)
        return ec
    end
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(chemicalelements(lipid.chain), b)
end
function getchemicalattr(lipid::NacylAmine, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(chemicalelements(lipid.chain), b)
end
function getchemicalattr(lipid::NacylAlkylAmine, ::Val{:elements}; kwargs...)
    b = reverse(chemicalelements(lipid.backbone))
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(chemicalelements(lipid.chain), b)
end
function getchemicalattr(lipid::FattyAcylEster, ::Val{:elements}; kwargs...)
    b = reverse(chemicalelements(lipid.backbone))
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(chemicalelements(lipid.chain), b)
end
function getchemicalattr(lipid::Monoradylglycerol, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(b, chemicalelements(lipid.chain))
end
function getchemicalattr(lipid::Diradylglycerol, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(b, chemicalelements.(lipid.chain)...)
end
function getchemicalattr(lipid::Triradylglycerol, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(b, chemicalelements.(lipid.chain)...)
end
function getchemicalattr(lipid::Union{Omodifiedmonoradylglycerol, Monoradylglycerophosphate, Lysophosphatidylserine, LysophosphatidylNmodifiedserine, Lysophosphatidylglycerol, Lysophosphatidylglycerolphosphate}, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(b, chemicalelements.(lipid.chain)...)
end
function getchemicalattr(lipid::Union{Omodifieddiradylglycerol, Diradylglycerophosphate, Phosphatidylserine, PhosphatidylNmodifiedserine, Phosphatidylglycerol, Phosphatidylglycerolphosphate}, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(b, chemicalelements.(lipid.chain)...)
end
function getchemicalattr(lipid::Union{Bisradylglycerophosphate, Bisradylglycerophosphoglycerol}, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findlast(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    if ncarbonchain(lipid) > 2 
        i = findlast(x -> ==(first(x), "H"), b)
        b[i] = "H" => (last(b[i]) - 1)
        filter!(x -> last(x) > 0, b)
    end
    if ncarbonchain(lipid) > 3 
        j = findfirst(x -> ==(first(x), "H"), b)
        i = findnext(x -> ==(first(x), "H"), b, j + 1)
        b[i] = "H" => (last(b[i]) - 1)
        filter!(x -> last(x) > 0, b)
    end
    vcat(b, chemicalelements.(lipid.chain)...)
end
function getchemicalattr(lipid::GlycerophosphoNacylethanolamine, ::Val{:elements}; kwargs...)
    b = chemicalelements(lipid.backbone)
    i = findfirst(x -> ==(first(x), "H"), b)
    b[i] = "H" => (last(b[i]) - 1)
    filter!(x -> last(x) > 0, b)
    vcat(chemicalelements(lipid.chain), b)
end
function getchemicalattr(lipid::SphingoidBaseBone, ::Val{:elements}; kwargs...)
    ec = chemicalelements(lipid.chain)
    i = findfirst(x -> ==(first(x), "N"), ec)
    ec[i + 1] = "H" => (last(ec[i + 1]) + 1)
    fc = lipid.chain
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if isnothing(lipid.headgroup) || !chain_species
        ec
    else
        ec[2] = "H" => (last(ec[2]) - 1)
        vcat(ec, chemicalelements(Substituent(Dehydroxy, lipid.headgroup)))
    end
end

function getchemicalattr(lipid::CeramideBone, ::Val{:elements}; kwargs...)
    ec = vcat(chemicalelements.(lipid.chain)...)
    chain = tuplize(lipid.chain)
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if isnothing(lipid.headgroup) || !chain_species
        ec
    else
        ec[2] = "H" => (last(ec[2]) - 1)
        vcat(ec, chemicalelements(Substituent(Dehydroxy, lipid.headgroup)))
    end
end

function getchemicalattr(lipid::MixSphingoBone, ::Val{:elements}; kwargs...)
    ec = vcat(chemicalelements.(lipid.chain)...)
    if lipid.chain isa CarbonChain{<: AbstractSPB}
        i = findfirst(x -> ==(first(x), "N"), ec)
        ec[i + 1] = "H" => (last(ec[i + 1]) + 1)
    end
    chain = tuplize(lipid.chain)
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if isnothing(lipid.headgroup) || !chain_species
        ec
    else
        ec[2] = "H" => (last(ec[2]) - length(lipid.headgroup))
        vcat(ec, [chemicalelements(Substituent(Dehydroxy, x)) for x in lipid.headgroup]...)
    end
end

getchemicalattr(lipid::Lipid, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(lipid); unique, kwargs...)

# MassSpecBioChemicals
parentchemical(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Acyl}) = FattyAcid(chain)
parentchemical(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00, AChirality())
leavinggroup(::CarbonChain{Acyl}) = Dehydroxy()
leavinggroup(::CarbonChain{SPB}) = Dehydrogen()
parentchemical(chain::CarbonChain{<: AbstractSTRing}) = SterolBone(chain)
leavinggroup(::CarbonChain{<: AbstractSTRing}) = Dehydrogen()

# FattyAcyl
# MassSpecBioChemicals
dehydroxyposition(::FattyAcid) = nothing
dehydrogenposition(::FattyAcid) = missing
dehydroxyposition(::FattyAlcohol) = nothing
dehydrogenposition(::FattyAlcohol) = nothing
dehydroxyposition(::FattyAmine) = nothing
dehydrogenposition(::FattyAmine) = nothing

# Lipid
# MassSpecBioChemicals.Lipids
chainposition(::Type{<: Hydrocarbon}) = ["hydrocarbon"]
chainposition(::Type{<: FattyAcid}) = ["acid"]
chainposition(::Type{<: FattyAlcohol}) = ["alcohol"]
chainposition(::Type{<: FattyAldehyde}) = ["acyl"]
chainposition(::Type{<: FattyAmide}) = ["nacyl"]
chainposition(::Type{<: FattyAmine}) = ["amine"]
chainposition(::Type{<: FattyAcylCarnitine}) = ["oacyl"]
chainposition(::Type{<: FattyAcylCoA}) = ["sacyl"]
chainposition(::Type{<: NacylAmine}) = ["amine", "nacyl"]
chainposition(::Type{<: FattyAcylEster}) = ["alcohol", "oacyl"]
chainposition(::Type{<: WaxEster}) = ["alcohol", "oacyl"]
chainposition(::Type{<: Glycerophospholipid}) = ["sn-1", "sn-2"]
chainposition(::Type{<: Sphingolipid}) = ["lcb", "nacyl"]
chainposition(::Type{<: Sphingolipid{H, <: CarbonChain{SPB}}}) where H = ["lcb"]
chainposition(::Type{<: Glycerolipid}) = ["sn-1", "sn-2", "sn-3"]
chainposition(::Type{<: Omodifiedradylglycerol}) = ["sn-1", "sn-2"]
chainposition(::Type{<: Bisradylglycerophosphoglycerol}) = ["sn-1", "sn-2", "sn-1'", "sn-2'"]
chainposition(::Type{<: Radyldiglycerol}) = ["sn-1", "sn-2", "sn-1'", "sn-2'"]
chainposition(::Type{<: Bisradylglycerophosphate}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
chainposition(::Type{<: Radyltriglycerol}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
chainposition(lipid::T) where {T <: Lipid} = chainposition(T)
nchainposition(::Type{T}) where {T <: Lipid} = length(chainposition(T))
nchainposition(lipid::T) where {T <: Lipid} = length(chainposition(T))
ncarbonchain(lipid::Lipid) = ncarbonchain(lipid.chain)

hascycloheadgroup(::Lipid) = false
hascycloheadgroup(::Sphingolipid) = false
hascycloheadgroup(::CeramidePhosphate) = true
hascycloheadgroup(::SphingoidBasePhosphate) = true

getlipidchain(lipid::Lipid) = tuplize(lipid.chain)
getlipidbody(lipid::Lipid) = tuplize(lipid.backbone)
getlipidbody(lipid::Sphingolipid) = tuplize(lipid.headgroup)

# CarbonChain
# MassSpecBioChemicals.Lipids
ncarbonchain(chains::Tuple) = sum(ncarbonchain, chains)
ncarbonchain(::CarbonChain{<: CarbonChainType}) = 1
ncarbonchain(::CarbonChain{<: T}) where {T <: Tuple} = length(T.parameters)
ncarbon(chains::Tuple) = sum(ncarbon, chains)
ncarbon(chain::CarbonChain) = chain.carbon
ndoublebond(chain::CarbonChain{S, UInt8}) where S = Int(chain.doublebond)
ndoublebond(chain::CarbonChain{S}) where S = length(chain.doublebond)

# interface MassSpecBioChemicals.Lipids for new lipid subclass
# getindex, length, size, eltype, getproperty, propertynames, setproperty!, keys, values, pairs
# interface for database
function Base.getindex(db::LipidDataBase, i)
    LipidDataBaseObject(db, i)
end

function Base.getindex(dbo::LipidDataBaseObject, i)
    dbo.database.records[i][dbo.id]
end

function Base.getindex(db::LipidDataBase, i, j)
    db.database.records[j][i]
end

Base.propertynames(db::LipidDataBase, private = false) = private ? (:records, (Symbol(k) for k in keys(db.records))...) : Tuple(Symbol(k) for k in keys(db.records))
Base.propertynames(dbo::LipidDataBaseObject, private = false) = private ? (:database, :id, (Symbol(k) for k in keys(dbo.database.records))...) : Tuple(Symbol(k) for k in keys(dbo.database.records))


function Base.getproperty(db::LipidDataBase, i)
    if i == :records
        getfield(db, i)
    else
        getfield(db, :records)[String(i)]
    end
end

function Base.getproperty(dbo::LipidDataBaseObject, i)
    if i == :database 
        getfield(dbo, i)
    elseif i == :id 
        getfield(dbo, i)
    else
        getfield(dbo, :database).records[String(i)][db.id]
    end
end

# Lazy iteration
Base.iterate(db::LipidDataBase, state = 1) = state > length(db) ? nothing : (LipidDataBaseObject(db, i), state + 1)

Base.keys(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.values(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.pairs(dbo::LipidDataBaseObject) = (k => getproperty(dbo, k) for k in propertynames(dbo))
function Base.iterate(dbo::LipidDataBaseObject, state = 1)
    p = propertynames(dbo)
    state > length(p) ? nothing : (getproperty(dbo, p[state]), state + 1)
end