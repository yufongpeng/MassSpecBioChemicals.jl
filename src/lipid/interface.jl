# CarbonChain
# MassSpecChemicals

# MassSpecBioChemicals
parentchemical(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Acyl}) = FattyAcid(chain)
parentchemical(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00)
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
# MassSpecChemicals
parse_chemical(::Type{<: Lipid}, s) = parse_lipid(s)
getchemicalattr(lipid::T, ::Val{:abbreviation}; kwargs...) where {T <: Lipid} = getchemicalattr(lipid, Val(:name); kwargs...)
function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Glycerolipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(class_abbr(lipid), " ", join(rp, "/"))
end
getchemicalattr(lipid::GlycerophosphoNacylethanolamine, ::Val{:name}; kwargs...) = string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Glycerophospholipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(class_abbr(lipid), " ", join(rp, "/"))
end

# function getchemicalattr(lipid::T) where {T <: FattyAcyl{B, <: Tuple}} where B
#     string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
# end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: FattyAcyl}
    string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Union{NacylAlkylAmine, FattyAcylEster}}
    if ncarbon(lipid.backbone.chain) == 0
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    else
        string(replace(chemicalname(lipid.backbone), class_abbr(lipid.backbone) => class_abbr(lipid)), "/", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    end
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: FattyAcylEstolid}
    if ncarbon(lipid.backbone.chain) == 0
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    elseif iszero(lipid.position)
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", replace(chemicalname(lipid.backbone), string(class_abbr(lipid.backbone), " ") => ""))
    else
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", Int(lipid.position), "O(", chemicalname(lipid.backbone), ")")
    end
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Sphingolipid{H, <: Tuple}} where H
    isnothing(lipid.headgroup) && return string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
    lv = annotationlevel(first(getlipidchain(lipid)); partial = true, additional = true, pass = true)
    if specieslevel in lv || molecularspecieslevel in lv
        return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    end
    position = decode_position(lipid)
    pos = any(==(0), position) ? "" : replace(string(position), " " => "", ",)" => ")")
    fc = deepcopy(first(getlipidchain(lipid)))
    del = Int[]
    position = vcat(collect.(position)...)
    if any(>=(structurepositionpartiallevel), lv)
        for i in position
            id = findfirst(==(i => Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            while id in del
                id = findnext(==(i => Hydroxy()), fc.substituent)
                isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
            push!(del, id)
        end
    elseif any(>=(structuredefinedpartiallevel), lv)
        for i in position
            id = findfirst(x -> ==(first(x) == Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            n = last(fc.substituent[id])
            if n > 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
            elseif n == 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
                push!(del, id)
            else
                throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
        end
    end
    deleteat!(fc.substituent, unique!(del))
    c = collect(getlipidchain(lipid))
    c[begin] = fc
    string(class_abbr(lipid), pos, " ", join(repr_singlechain.(c), "/"))
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Sphingolipid{H, <: CarbonChain}} where H
    isnothing(lipid.headgroup) && return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    lv = annotationlevel(first(getlipidchain(lipid)); partial = true, additional = true, pass = true)
    if specieslevel in lv || molecularspecieslevel in lv
        return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    end
    position = decode_position(lipid)
    pos = any(iszero, position) ? "" : replace(string(position), " " => "", ",)" => ")")
    fc = deepcopy(first(getlipidchain(lipid)))
    del = Int[]
    position = vcat(collect.(position)...)
    if any(>=(structurepositionpartiallevel), lv)
        for i in position
            id = findfirst(==(i => Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            while id in del
                id = findnext(==(i => Hydroxy()), fc.substituent)
                isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
            push!(del, id)
        end
    elseif any(>=(structuredefinedpartiallevel), lv)
        for i in position
            id = findfirst(x -> first(x) == Hydroxy(), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            n = last(fc.substituent[id])
            if n > 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
            elseif n == 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
                push!(del, id)
            else
                throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
        end
    end
    deleteat!(fc.substituent, unique!(del))
    string(class_abbr(lipid), pos, " ", repr_singlechain(fc))
end


function getchemicalattr(lipid::Monoradylglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    p = only(decode_sn(lipid))
    c = chemicalsmiles_carbonchain(lipid.chain)
    smi = chemicalsmiles(lipid.backbone)
    r = collect(eachmatch(r"\(O\)", smi))
    next = r[lastindex(r) - p + 1].match.offset
    string(smi[begin:next + 2], c, smi[next + 3:end])
end
function getchemicalattr(lipid::Diradylglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    position = decode_sn(lipid)
    p = sortperm(position; rev = true)
    position = position[p]
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain[p]]
    smi = chemicalsmiles(lipid.backbone)
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[[length(r) - t + 1 for t in position]], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

function getchemicalattr(lipid::Ceramide, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain]
    id = findfirst("(N)", first(cs))
    string(first(cs)[begin:id[2]], cs[begin + 1], first(cs)[last(id):end])
end

getchemicalattr(lipid::SphingoidBase, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles_carbonchain(lipid.chain) : ""
getchemicalattr(lipid::CeramideBone, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_headgroup(lipid)) : ""
getchemicalattr(lipid::SphingoidBaseBone, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_headgroup(lipid)) : ""

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
chainposition(::Type{<: Bisradylglycerophosphate}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
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