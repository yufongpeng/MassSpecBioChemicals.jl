# Chain position handling shared by `parse_carbonchain` and lipid constructors
# A position specification is
# - `:start` or `:noorder` => 0x00
# - `:inorder` => index of the chain
# - a string such as "(sn-1)" or "sn-1" => index in `chainposition(T)`
# - an integer => position value

is_explicit_sn(position) = position isa AbstractString || (position isa Integer && position != 0)

function uint8ize_sn(::Type{T}, position::Integer, index, label) where {T <: Lipid}
    0 <= position <= nchainposition(T) || throw(ArgumentError("Invalid chain position, \"$label\""))
    UInt8(position)
end
function uint8ize_sn(::Type{T}, position::Symbol, index, label) where {T <: Lipid}
    position === :start || position === :noorder ? 0x00 :
        position === :inorder ? UInt8(index) : throw(ArgumentError("Invalid chain position, \"$label\""))
end
function uint8ize_sn(::Type{T}, position::AbstractString, index, label) where {T <: Lipid}
    m = match(r"^\(?\s*(sn-\d+'*)\s*\)?$", position)
    isnothing(m) && throw(ArgumentError("Invalid chain position, \"$label\""))
    i = findfirst(==(first(m.captures)), chainposition(T))
    isnothing(i) && throw(ArgumentError("Invalid chain position, \"$label\""))
    0 <= i <= nchainposition(T) || throw(ArgumentError("Invalid chain position, \"$label\""))
    UInt8(i)
end
uint8ize_sn(::Type{T}, positions, labels) where {T <: Lipid} =
    UInt8[uint8ize_sn(T, p, i, labels[i]) for (i, p) in enumerate(positions)]

function reorder_lipid_chains!(order, arrays...)
    for a in arrays
        a[:] = a[order]
    end
end

# Chains with explicit positions are moved to the front, and the others keep their relative order
function move_explicit_lipid_chains!(positions, arrays...)
    explicit = findall(is_explicit_sn, positions)
    isempty(explicit) && return explicit
    order = vcat(explicit, setdiff!(collect(eachindex(positions)), explicit))
    reorder_lipid_chains!(order, positions, arrays...)
    explicit
end

function resolve_chain_order!(positions)
    length(positions) > 1 && any(==(:inorder), @view positions[begin + 1:end]) || return
    allequal(positions[begin + 1:end]) ||
        throw(ArgumentError("Chains should be all in order or not in order"))
    positions[begin] = :inorder
end

function sort_explicit_lipid_chains!(values, explicit, chains, arrays...)
    nexplicit = length(explicit)
    comparable(i) = chains[i] isa CarbonChain
    function lt(i, j)
        (i <= nexplicit) != (j <= nexplicit) && return i <= nexplicit
        values[i] != values[j] && return values[i] < values[j]
        values[i] == 0x00 && comparable(i) && comparable(j) && return isless(chains[i], chains[j])
        false
    end
    reorder_lipid_chains!(sortperm(eachindex(values); lt, alg = MergeSort), values, chains, arrays...)
end

# `positions` is either the merged code or position values of each chain
function encode_sn(::Type{T}, positions) where {T <: Lipid}
    npositions = nchainposition(T)
    code = npositions > 3 ? UInt16(0) : UInt8(0)
    if positions isa Integer
        0 <= positions <= typemax(code) || throw(ArgumentError("Invalid chain position code, $positions"))
        return typeof(code)(positions)
    end
    any(>(npositions), positions) && throw(ArgumentError("Chain position exceeds the number of positions for $T"))
    allunique(filter(!=(0x00), positions)) || throw(ArgumentError("Overlapped chain position, $positions"))
    base = convert(typeof(code), npositions + 1)
    for position in positions
        code = code * base + position
    end
    code
end

"""
    arrange_lipid_chains(T, chain, sn)

Sort chains by sn-position(s); return the arranged chain(s) and encoded position code.

`sn` can be the merged code, a symbol, a string like "(sn-1)", an integer, or a collection of them for each chain.
"""
function arrange_lipid_chains(::Type{T}, chain, sn) where {T <: Lipid}
    chains = collect(Any, tuplize(chain))
    if sn isa Integer
        # the merged code is relative to the given chains, so empty chains are kept
        return chain, encode_sn(T, sn)
    end
    n = length(chains)
    specs = collect(Any, tuplize(sn))
    length(specs) == n || throw(ArgumentError("Expected $n chain position(s), got $(length(specs))"))
    # sn isa AbstractString && n > 1 && throw(ArgumentError("Expected one position for each chain, got \"$sn\""))
    explicit = move_explicit_lipid_chains!(specs, chains)
    resolve_chain_order!(specs)
    values = uint8ize_sn(T, specs, chains)
    sort_explicit_lipid_chains!(values, explicit, chains)
    # drop empty chains and their positions as `parse_carbonchain`
    keep = findall(c -> ncarbon(c) > 0, chains)
    chains, values = chains[keep], values[keep]
    code = encode_sn(T, values)
    (length(chains) == 1 ? first(chains) : Tuple(chains)), code
end

function encode_position(position::Integer)
    position >= 0 || throw(ArgumentError("Invalid headgroup position: $position"))
    UInt8(position)
end
# encode_headgroup_position(position::Symbol) =
#     position in (:start, :noorder) ? 0x00 : throw(ArgumentError("Invalid headgroup position: $position"))

mix_lipid_chirality(chirality) = RSSystem[chirality...]

function encode_position_chirality(::Type{<:SphingoBone}, bone, pos, chirality)
    (isnothing(pos) || isnothing(bone) || bone isa GlyComp) && return Pair{UInt8, RSSystem}[0x00 => RSChirality()], true
    if pos isa Integer 
        pos = (pos, )
    end
    if chirality isa RSSystem || chirality isa AbstractString
        chirality = (chirality, )
    end
    if !isnothing(chirality)
        length(pos) == length(chirality) || 
            length(chirality) == 1 || 
            throw(ArgumentError("Mismatched number of positions and chiralities"))
        encode_position_chirality(pos, chirality), false
    else
        encode_position_chirality(pos), false
    end
end

function encode_position_chirality(::Type{<:MixSphingoBone}, bone, pos, chirality)
    isnothing(bone) && return UInt8[], true
    isnothing(pos) && return [[0x00 => RSChirality()] for i in eachindex(bone)], true
    if !isnothing(chirality)
        length(pos) == length(chirality) || throw(ArgumentError("Mismatched number of positions and chiralities"))
        [encode_position_chirality(p, c) for (p, c) in zip(pos, chirality)], false
    else
        [encode_position_chirality(p) for p in pos], false
    end
end

encode_position_chirality(p::Integer, c::Union{RSSystem, AbstractString}) = [encode_position(p) => parse_chirality(c)]
encode_position_chirality(ps, c::Union{RSSystem, AbstractString}) = [encode_position(p) => parse_chirality(c) for p in ps]
# encode_position_chirality(p::Integer, cs) = [encode_position(p) => parse_chirality(c) for c in cs]
encode_position_chirality(ps, cs) = [encode_position(p) => parse_chirality(c) for (p, c) in zip(ps, cs)]

encode_position_chirality(p::Union{Integer, AbstractString}) = [parse_headposition(p)]
encode_position_chirality(ps) = [parse_headposition(p) for p in ps]