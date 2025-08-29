
"""
    split_chainmodification(mod)

Split chain modification into vector by ";"
"""
split_chainmodification(::Nothing) = String[]
function split_chainmodification(mod)
    s = Int[]
    e = Int[]
    dp = 0
    db = 0
    for i in eachindex(mod)
        @inbounds x = mod[i]
        if dp == 0 && db == 0 && x == ';'
            push!(s, i)
            push!(e, i)
            continue
        end
        @match x begin
            '(' => (dp += 1)
            '[' => (db += 1)
            ')' => (dp -= 1)
            ']' => (db -= 1)
            _   => () 
        end
        e[end] = i
    end
    @inbounds [mod[i:j] for (i, j) in zip(s, e)]
end

"""
    split_chainmodification_c(mod)

Split chain modification into vector by ","
"""
split_chainmodification_c(::Nothing) = String[]
function split_chainmodification_c(mod)
    s = Int[]
    e = Int[]
    dp = 0
    db = 0
    for i in eachindex(mod)
        @inbounds x = mod[i]
        if dp == 0 && db == 0 && (x == ',' || x == ';')
            push!(s, i)
            push!(e, i)
            continue
        end
        @match x begin
            '(' => (dp += 1)
            '[' => (db += 1)
            ')' => (dp -= 1)
            ']' => (db -= 1)
            _   => () 
        end
        e[end] = i
    end
    @inbounds [mod[i:j] for (i, j) in zip(s, e)]
end

"""
    parse_chainmodification(T, mod, ox = nothing)

Parse chain modification for full structure level
"""
function parse_chainmodification(::Type{<: Radyl}, mod, ox = nothing)
    isempty(mod) && isnothing(ox) && return (nothing, Pair{UInt8, RSSystem}[])
    i = findfirst(x -> !isnothing(match(r"^;O\d*$", x)), mod)
    if isnothing(ox) && isnothing(i)
        split_fg_chirality!(vcat((parse_definedmodification(x) for x in mod)...))
    elseif isnothing(ox)
        (sort_chainmodification!([parse_speciesmodification(x) for x in mod]), Pair{UInt8, RSSystem}[])
    elseif isnothing(i)
        push!(mod, ox)
        (sort_chainmodification!([parse_speciesmodification(x) for x in mod]), Pair{UInt8, RSSystem}[])
    else
        throw(ArgumentError("Duplicated oxygen atom in chain modification"))
    end
end

function parse_chainmodification(::Type{<: AbstractSPB}, mod, ox)
    # except O-
    isempty(mod) && isnothing(ox) && return (nothing, Pair{UInt8, RSSystem}[])
    i = findfirst(x -> !isnothing(match(r"^;O\d*$", x)), mod)
    if isnothing(ox) && isnothing(i)
        split_fg_chirality!(vcat((parse_definedmodification(x) for x in mod)...))
    elseif isnothing(ox)
        (sort_chainmodification!([parse_speciesmodification(x; onlinked = false) for x in mod]), Pair{UInt8, RSSystem}[])
    elseif isnothing(i)
        push!(mod, ox)
        (sort_chainmodification!([parse_speciesmodification(x; onlinked = false) for x in mod]), Pair{UInt8, RSSystem}[])
    else
        throw(ArgumentError("Duplicated oxygen atom in chain modification"))
    end
    # mod = first(mod).match
    # m = match(r"^;O(\d*)$", mod)
    # isnothing(m) || return (isnothing(first(m)) ? 0x01 : parse(UInt8, first(m)))
    # m = match(r"^;((\d+)OH,?)+$", mod)
    # if isnothing(m)
    #     m = match(r"^;\(+OH\)+(\d*)$", mod)
    #     isnothing(m) && throw(ArgumentError("Invalid SPB modification."))
    #     n, = m
    #     [Hydroxy() => isnothing(n) ? 0x01 : parse(UInt8, n)]
    # else
    #     m = eachmatch(r"(\d+)OH", mod)
    #     [parse(UInt8, first(x.captures)) => Hydroxy() for x in m]
    # end
end

function parse_chainmodification(::Type{T}, mod) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_chainmodification` not implemented for `$T`"))
end

"""
    parse_speciesmodification(mod; onlinked = true)

Parse chain modification of species level
`onlinked`: allow O-linked or N-linked modification, e.g. "OMe", "O(FA 16:0)".
"""
function parse_speciesmodification(mod; onlinked = true)
    FG_SP = onlinked ? FG_SPECIES : FG_SPECIES_CLINKED
    for (k, v) in FG_SP
        m = match(Regex(string("^;\\(?", k, "\\)?(\\d*)\$")), mod)
        isnothing(m) && continue
        return v => (isempty(first(m)) ? 0x01 : parse(UInt8, first(m)))
    end
    # multiple layer?
    m = startswith(mod, ";(") ? match(r"^;\((C?[ON]?)\((.*)\)\)(\d*)$", mod) : match(r"^;(C?[ON]?)\((.*)\)\)$", mod) 
    if isnothing(m)
        m = startswith(mod, ";(") ? match(r"^;\((C?[ON]?)(.*)\)(\d*)$", mod) : match(r"^;(C?[ON]?)(.*)\)$", mod) 
    end
    isnothing(m) && throw(ArgumentError("Invalid chain modification, \"$mod\$"))
    !onlinked && !isempty(l) && throw(ArgumentError("O-linked or N-linked modification not allowed, \"$mod\$"))
    l, m, n = m
    m = startswith(m, r"\d+:\d+") ? string("FOH ", m) : m
    c = parse_lipid(m)
    c = isempty(l) ? Substituent(Dehydroxy, c, lk(dehydroxyposition(c))) : 
        l == "O" ? XLinkedFunctionalGroup(OLinkage(), Substituent(Dehydroxy, c, lk(dehydroxyposition(c)))) : 
            XLinkedFunctionalGroup(NLinkage(), Substituent(Dehydroxy, c, lk(dehydroxyposition(c))))
    c => (isempty(n) ? 0x01 : parse(UInt8, n))
end

parse_rschirality(c::AbstractString) = c == "R" ? RChirality() : c == "S" ? SChirality() : throw(ArgumentError("Unknown chirality `$c`"))
function parse_structureposition(v, x) 
    p1, p2, c = x.captures
    p1 = parse(UInt8, p1)
    p2 = (isnothing(p2) || isempty(p2)) ? p1 : parse(UInt8, p2)
    c = isnothing(c) ? RSChirality() : parse_rschirality(c)
    ([p1 => v], [p2 => c])
end
function parse_structureposition(v::Union{Epoxy, Peroxy}, x) 
    p1, p2, c = x.captures
    p1 = parse(UInt8, p1)
    p2, p3 = (isnothing(p2) || isempty(p2)) ? (p1, p1 + 0x01) : (p1, parse(UInt8, p2))
    c = isnothing(c) ? RSChirality() : parse_rschirality(c)
    ([p1 => v, p3 => Hydrogen()], [p2 => c, p3 => RSChirality()])
end
function parse_structureposition(v::Oxo, x) 
    p1, p2, c = x.captures
    p1 = parse(UInt8, p1)
    p2 = (isnothing(p2) || isempty(p2)) ? p1 : parse(UInt8, p2)
    isnothing(c) || @warn "Oxo group has no chiral center."
    ([p1 => v], [p2 => AChirality()])
end
# nCOX[n-1R/S] = C_(n-1)-COX
"""
    parse_definedmodification(mod)

Parse chain modification of defined or full structure level, eg ";2OH,3OH"
"""
function parse_definedmodification(mod; onlinked = true)
    FG_SP = onlinked ? FG : FG_CLINKED 
    ms = split_chainmodification_c(mod)
    if startswith(first(ms), r";\d")
        # position 
        # config
        for (k, v) in FG_SP
            m = match.(Regex(string("^[;,](\\d+)", k, "(?:\\[(\\d*)([RS]*)\\])*\$")), ms)
            all(isnothing, m) && continue
            any(isnothing, m) && throw(ArgumentError("Some chain modification does not have position, \"$mod\""))
            return [parse_structureposition(v, x) for x in m]
        end
        mm = match(r"^[;,]\d+(C?[ON]?)\((.*)\)(?:\[\d*[RS]*\])*$", first(ms))
        v = ""
        if isnothing(mm)
            mm = match(r"^[;,]\d+(C?[ON]?)(.*[^\]])(?:\[\d*[RS]*\])*$", first(ms))
        else
            v = string(first(mm.captures), "(", last(mm.captures), ")")
        end
        isnothing(mm) && throw(ArgumentError("Invalid chain modification, \"$mod\""))
        v = isempty(v) ? string(mm...) : v
        v = replace(v, "(" => "\\(", ")" => "\\)", "[" => "\\[", "]" => "\\]")
        c = parse_tailsubstituent(mm...; onlinked)
        m = match.(Regex(string("^[;,](\\d+)", v, "(?:\\[(\\d*)([RS]*)\\])*\$")), ms)
        any(isnothing, m) && throw(ArgumentError("Some chain modification does not have position, \"$mod\""))
        return [parse_structureposition(c, x) for x in m]
    elseif length(ms) == 1
        for (k, v) in FG_SP
            m = match(Regex(string("^;\\(?", k, "\\)?(\\d*)\$")), mod)
            isnothing(m) ? continue : return [v => (isempty(first(m)) ? 0x01 : parse(UInt8, first(m)))]
        end
        m = startswith(first(ms), ";(") ? match(r"^;\((C?[ON]?)\((.*)\)\)(\d*)$", first(ms)) : match(r"^;(C?[ON]?)\((.*)\)()$", first(ms)) 
        isnothing(m) && throw(ArgumentError("Invalid chain modification, \"$mod\""))
        l, fg, n = m
        return [parse_tailsubstituent(l, fg) => (isempty(n) ? 0x01 : parse(UInt8, n))]
    else
        throw(ArgumentError("Invalid chain modification, \"$mod\""))
    end
    # [cyc], multiple layer?
end

sort_chainmodification!(::Nothing) = nothing

function sort_chainmodification!(mod::Vector{<: Tuple})
    sort!(mod; by = x -> (sub_abbr(last(first(x))), first(first(x))))
end  

# delete duplicate
function sort_chainmodification!(mod::AbstractVector{<: Pair{<: Number}})
    sort!(mod; by = x -> (sub_abbr(last(x)), first(x)))# p => fg
    del = Int[]
    p = 0x00
    for (i, x) in enumerate(mod)
        if last(x) == Hydrogen() && p == first(x)
            push!(del, i)
        elseif last(x) == Hydrogen()
            p = first(x)
        end
    end
    deleteat!(mod, del)
end

function sort_chainmodification!(mod)        
    # fg => n
    m = Dict{eltype(mod).parameters...}()
    for (k, v) in mod 
        get!(m, k, 0x00)
        m[k] += v 
    end
    sort!(filter!(x -> last(x) > 0, collect(pairs(m))); by = x -> (sub_abbr(first(x)), last(x)))
    # del = Int[]
    # for (i, v) in enumerate(mod)
    #     n = m[first(v)]
    #     if n == 0 
    #         push!(del, i)
    #         continue
    #     elseif n > 0 && last(v) != n
    #         mod[i] = first(v) => n 
    #     end 
    #     m[first(v)] = 0x00
    # end
    # sort!(deleteat!(mod, del); by = x -> (sub_abbr(first(x)), last(x)))
end

split_fg_chirality!(x::Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}}) = sort_chainmodification!(convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, x)), Pair{UInt8, RSSystem}[]
#split_fg_chirality(x::Vector{<: Tuple}) = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, first.(x)), convert(Vector{Pair{UInt8, RSSystem}}, last.(x))
function split_fg_chirality!(x::Vector{<: Tuple}) 
    cd = Dict{UInt8, RSSystem}()
    for cs in last.(x)
        for (p, c) in cs 
            get!(cd, p, RSChirality())
            if cd[p] == RSChirality()
                cd[p] = c 
            elseif cd[p] != c && c != RSChirality()
                throw(ArgumentError("Conflicts of chirality at position `$(Int(p))`"))
            end
        end
    end
    sort_chainmodification!(convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, vcat(first.(x)...))), collect(pairs(cd))
end
