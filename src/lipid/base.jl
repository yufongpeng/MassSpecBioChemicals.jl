function Base.getindex(db::LipidDataBase, i)
    LipidDataBaseObject(db, i)
end

function Base.getindex(dbo::LipidDataBaseObject, i)
    dbo.database.records[i][dbo.id]
end

function Base.getindex(db::LipidDataBase, i, j)
    db.database.records[j][i]
end

Base.propertynames(db::LipidDataBase, private = false) =
    private ? (:records, (Symbol(k) for k in keys(db.records))...) : Tuple(Symbol(k) for k in keys(db.records))
Base.propertynames(dbo::LipidDataBaseObject, private = false) =
    private ? (:database, (Symbol(k) for k in keys(dbo.database.records))...) : Tuple(Symbol(k) for k in keys(dbo.database.records))

function Base.getproperty(db::LipidDataBase, i)
    if i == :records
        getfield(db, i)
    else
        getfield(db, :records)[String(i)]
    end
end

function Base.getproperty(dbo::LipidDataBaseObject, i)
    if i == :database || i == :id
        getfield(dbo, i)
    else
        getfield(dbo, :database).records[String(i)][dbo.id]
    end
end

Base.iterate(db::LipidDataBase, state = 1) =
    state > length(db) ? nothing : (LipidDataBaseObject(db, state), state + 1)
Base.length(db::LipidDataBase) = length(first(values(getfield(db, :records))))

Base.keys(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.values(dbo::LipidDataBaseObject) = (k for k in propertynames(dbo))
Base.pairs(dbo::LipidDataBaseObject) = (k => getproperty(dbo, k) for k in propertynames(dbo))
function Base.iterate(dbo::LipidDataBaseObject, state = 1)
    p = propertynames(dbo)
    state > length(p) ? nothing : (getproperty(dbo, p[state]), state + 1)
end

const UNKNOWN_CHAIN_POSITION = typemax(UInt8)

_chain_type(::CarbonChain{T}) where T = T
_chain_type(::CarbonChain{T}) where {T <: Tuple} = first(T.parameters)

function _chain_type_key(T)
    T <: AbstractSTRing ? (0, string(nameof(T))) :
    T === SulfoSPB ? (1, "") :
    T === SPB ? (2, "") :
    T === Alkenyl{ZConfiguration} ? (3, "") :
    T === Alkenyl{EConfiguration} ? (4, "") :
    T <: Alkenyl ? (5, "") :
    T <: Alkyl ? (6, "") :
    T === Acyl ? (7, "") : (8, string(nameof(T)))
end

_db_config_rank(::ZConfiguration) = 0
_db_config_rank(::EConfiguration) = 1
_db_config_rank(_) = 2
_rs_config_rank(::RChirality) = 0
_rs_config_rank(::SChirality) = 1
_rs_config_rank(_) = 2

# (minimal position, configuration rank); unknown position is the largest
function _min_db(db::Vector{<: Pair})
    ps = filter(x -> first(x) != 0x00, db)
    isempty(ps) && return UNKNOWN_CHAIN_POSITION, 2
    p = argmin(first, ps)
    first(p), _db_config_rank(last(p))
end
_min_db(_) = UNKNOWN_CHAIN_POSITION, 2

_is_chain_sub(x) = !(last(x) isa Hydrogen)
# (type: 0 known position, 1 unknown position; number; minimal position)
_sub_info(::Nothing) = 0, 0, UNKNOWN_CHAIN_POSITION
_sub_info(sub::Integer) = sub == 0 ? 0 : 1, Int(sub), UNKNOWN_CHAIN_POSITION
function _sub_info(sub::Vector{<: Pair{<: Integer}})
    ps = filter(_is_chain_sub, sub)
    known = filter(x -> first(x) != 0x00, ps)
    isempty(ps) && return 0, 0, UNKNOWN_CHAIN_POSITION
    length(known) == length(ps) ? 0 : 1, length(ps), isempty(known) ? UNKNOWN_CHAIN_POSITION : UInt8(minimum(first, known))
end
function _sub_info(sub::Vector{<: Pair})
    n = sum(x -> Int(last(x)), sub; init = 0)
    n == 0 ? 0 : 1, n, UNKNOWN_CHAIN_POSITION
end

function _sub_config_rank(chirality, p)
    p == UNKNOWN_CHAIN_POSITION && return 2
    i = findfirst(x -> first(x) == p, chirality)
    isnothing(i) ? 2 : _rs_config_rank(last(chirality[i]))
end

function _chain_sort_key(chain::CarbonChain)
    tr, tn = _chain_type_key(_chain_type(chain))
    dbp, dbc = _min_db(chain.doublebond)
    st, sn, sp = _sub_info(chain.substituent)
    sc = chain.chirality isa Vector ? _sub_config_rank(chain.chirality, sp) : 2
    (ncarbonchain(chain), tr, tn, Int(chain.carbon), ndoublebond(chain), Int(dbp), dbc, st, sn, Int(sp), sc)
end

"""
    isless(chain1::CarbonChain, chain2::CarbonChain)

Compare chains by number of chain types (fewer first), type of the first chain, chain length, number of double bonds, minimal double bond position,
its configuration, substitution type, substitution number, minimal substitution position and its configuration.
"""
Base.isless(a::CarbonChain, b::CarbonChain) = isless(_chain_sort_key(a), _chain_sort_key(b))

function ==(x::T, y::T) where {T <: AbstractCarbonChain}
    for f in fieldnames(T)
        a = getfield(x, f)
        b = getfield(y, f)
        eq = f === :elements ? chemicalelements_equal(a, b) : a == b
        eq || return false
    end
    true
end

function hash(x::T, h::UInt) where {T <: AbstractCarbonChain}
    h = hash(T, h)
    for f in fieldnames(T)
        v = getfield(x, f)
        h = f === :elements ? hash(sum(hash, sort_unique_elements(v); init = zero(h)), h) : hash(v, h)
    end
    h
end

copy(x::T) where {T <: AbstractCarbonChain} = fieldcount(T) == 0 ? x : T((safecopy(getfield(x, f)) for f in fieldnames(T))...)

function Base.show(io::IO, level::AnnotationLevel)
    print(io, @match level begin
        &specieslevel                   => "specieslevel"
        &molecularspecieslevel          => "molecularspecieslevel"
        &phosphatepositionlevel         => "phosphatepositionlevel"
        &snpositionlevel                => "snpositionlevel"
        &dbpositionpartiallevel         => "dbpositionpartiallevel"
        &dbpositionlevel                => "dbpositionlevel"
        &dbconfigpartiallevel           => "dbconfigpartiallevel"
        &dbconfiglevel                  => "dbconfiglevel"
        &structuredefinedpartiallevel   => "structuredefinedpartiallevel"
        &structuredefinedlevel          => "structuredefinedlevel"
        &structurepositionpartiallevel  => "structurepositionpartiallevel"
        &structurepositionlevel         => "structurepositionlevel"
        &structureconfigpartiallevel    => "structureconfigpartiallevel"
        &structureconfiglevel           => "structureconfiglevel"
        &fullstructurelevel             => "fullstructurelevel"
        &completestructurelevel         => "completestructurelevel"
    end)
end

function Base.show(io::IO, level::PassLevel)
    print(io, @match level begin
        &passphosphatepositionlevel     => "passphosphatepositionlevel"
        &passsnpositionlevel            => "passsnpositionlevel"
    end)
end