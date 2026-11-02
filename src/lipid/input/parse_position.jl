
function parse_headposition(::Type{<: SphingoBone}, bone, pos)
    if isnothing(pos) || isnothing(bone) || bone isa GlyComp
        Pair{UInt8, RSSystem}[0x00 => RSChirality()], true
    else
        Pair{UInt8, RSSystem}[parse_headposition(x) for x in eachmatch(r"(\d+)([RS]*)", pos)], false
    end
end

function parse_headposition(::Type{<: MixSphingoBone}, bone, pos)
    if isnothing(bone)
        UInt8[], true
    elseif isnothing(pos)
        [[0x00 => RSChirality()] for i in eachindex(bone)], true
    else
        [[parse_headposition(x) for x in eachmatch(r"(\d+)([RS]*)", p)] for p in split_class_position(pos)], false
    end
end

parse_headposition(::Type{<:Lipid}, bone, pos) = pos

"""
    split_class_position(s)
Split string into position(s) for headgroup(s)
"""
function split_class_position(pos)
    s = nextind(pos, firstindex(pos))
    e = firstindex(pos)
    dp = -1
    v = Any[]
    for i in eachindex(pos)
        @inbounds x = pos[i]
        if x == '('
            dp += 1
        elseif dp == 0 && x != ',' && x != ')'
            e = i
            continue
        elseif dp == 0
            push!(v, pos[s:e])
            s = nextind(pos, i)
            continue
        elseif x == ')'
            dp -= 1
        end
        e = i
    end
    v
end

function parse_headposition(x)
    p, c = x 
    parse(UInt8, p) => parse_chirality(c)
end
parse_headposition(x::Integer) = encode_position(x) => RSChirality()

parse_chirality(x) = isempty(x) ? RSChirality() : x == "R" ? RChirality() : x == "S" ? SChirality() : throw(ArgumentError("Invalid chirality."))
parse_chirality(x::RSChirality) = x