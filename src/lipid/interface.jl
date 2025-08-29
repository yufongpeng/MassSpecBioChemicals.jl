# CarbonChain
# MassSpecChemicals

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