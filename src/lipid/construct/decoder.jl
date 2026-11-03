"""
    decode_sn(::Type{<:Lipid}, nchain, snpositioncode)

Decode sn-position code into a vector of sn-positions.
"""
function decode_sn(::Type{T}, count, sn) where {T <: Lipid}
    positions = zeros(Int, count)
    base = nchainposition(T) + 1
    for i in eachindex(positions)
        place = base ^ (count - i)
        positions[i], sn = divrem(sn, place)
    end
    positions
end

"""
    decode_headgroup_position(lipid)

Decode headgroup-position code into headgroup position representation.
"""
function decode_headgroup_position(::Type{T}, position, chirality) where {T <: SphingoBone}
    p = Int(position)
    if length(tuplize(chirality)) == 2
        p1, p2 = divrem(p, 32)
        if iszero(p1) && iszero(p2)
            ""
        elseif iszero(p1) 
            rc = repr_chirality(last(chirality); bracket = false)
            string("(", p2, rc, ")") 
        else
            rc = repr_chirality.(chirality; bracket = false)
            string("(", p1, rc[1], ",", p2, rc[2], ")")
        end
    elseif iszero(p)
        ""
    else
        rc = repr_chirality(chirality; bracket = false)
        string("(", p, rc, ")") 
    end
end

decode_headgroup_position(::Type{T}, position, chirality) where {T <: SphingoBone{Nothing}} = ""

function decode_headgroup_position(::Type{T}, position, chirality) where {T <: MixSphingoBone}
    ps = Int.(position)
    s = Vector{String}(undef, length(ps))
    for (i, (p, c)) in enumerate(zip(ps, chirality))
        if length(tuplize(c)) == 2
            p1, p2 = divrem(p, 32)
            s[i] = if iszero(p1) && iszero(p2)
                ""
            elseif iszero(p1) 
                rc = repr_chirality(c; bracket = false)
                string(p2, rc) 
            else
                rc = repr_chirality.(c; bracket = false)
                string(p1, rc[1], ",", p2, rc[2])
            end
        elseif iszero(p)
            s[i] = ""
        else
            rc = repr_chirality(c; bracket = false)
            s[i] = string(p, rc) 
        end
    end
    if all(isempty, s)
        "" 
    else
        string("(", join(s, ","), ")")
    end
end