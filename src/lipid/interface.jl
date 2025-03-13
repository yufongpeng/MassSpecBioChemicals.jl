originalmolecule(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Acyl}) = FattyAcid(chain)
originalmolecule(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00)
leavinggroup(::CarbonChain{Acyl}) = Dehydroxy()
leavinggroup(::CarbonChain{SPB}) = Dehydrogen()
originalmolecule(chain::CarbonChain{<: AbstractSTRing}) = SterolBone(chain)
leavinggroup(::CarbonChain{<: AbstractSTRing}) = Dehydrogen()

# interface MassSpecBioChemicals.Lipids for new lipid subclass
# getindex, length, size, eltype, getproperty, propertynames, setproperty!, keys, values, pairs
# interface for database
function Base.getindex(db::LipidDatabase, i)
    LipidDataBaseObject(db, i)
end

function Base.getindex(dbo::LipidDataBaseObject, i)
    dbo.database.records[i][dbo.id]
end

function Base.getindex(db::LipidDatabase, i, j)
    db.database.records[j][i]
end

Base.propertynames(db::LipidDatabase, private = false) = private ? (:records, (Symbol(k) for k in keys(db.records))...) : Tuple(Symbol(k) for k in keys(db.records))
Base.propertynames(dbo::LipidDataBaseObject, private = false) = private ? (:database, :id, (Symbol(k) for k in keys(dbo.database.records))...) : Tuple(Symbol(k) for k in keys(dbo.database.records))


function Base.getproperty(db::LipidDatabase, i)
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