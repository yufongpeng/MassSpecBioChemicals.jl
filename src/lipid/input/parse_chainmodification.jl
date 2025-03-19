
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
    isempty(mod) && isnothing(ox) && return nothing
    i = findfirst(x -> !isnothing(match(r"^;O\d*$", x)), mod)
    if isnothing(ox) && isnothing(i)
        pc = vcat((parse_definedmodification(x) for x in mod)...)
        if (eltype(pc) <: Pair{UInt8, <: AbstractFunctionalGroup}) || (eltype(pc) <: Pair{<: AbstractFunctionalGroup, UInt8})
            sort_chainmodification!(pc)
        else
            throw(ArgumentError("Chain modification should be in either defined structure level or full structure level"))
        end
    elseif isnothing(ox)
        sort_chainmodification!([parse_speciesmodification(x) for x in mod])
    elseif isnothing(i)
        push!(mod, ox)
        sort_chainmodification!([parse_speciesmodification(x) for x in mod])
    else
        throw(ArgumentError("Duplicated oxygen atom in chain modification"))
    end
end

function parse_chainmodification(::Type{<: SPB}, mod, ox)
    # except O-
    isempty(mod) && isnothing(ox) && return nothing
    i = findfirst(x -> !isnothing(match(r"^;O\d*$", x)), mod)
    if isnothing(ox) && isnothing(i)
        sort_chainmodification!(vcat((parse_definedmodification(x) for x in mod)...))
    elseif isnothing(ox)
        sort_chainmodification!([parse_speciesmodification(x; onlinked = false) for x in mod])
    elseif isnothing(i)
        push!(mod, ox)
        sort_chainmodification!([parse_speciesmodification(x; onlinked = false) for x in mod])
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

"""
    parse_definedmodification(mod)

Parse chain modification of defined or full structure level, eg ";2OH,3OH"
"""
function parse_definedmodification(mod; onlinked = true)
    FG_SP = onlinked ? FG : FG_CLINKED 
    ms = split_chainmodification_c(mod)
    if startswith(first(ms), r";\d")
        for (k, v) in FG_SP
            m = match.(Regex(string("^[;,](\\d+)", k, "\$")), ms)
            all(isnothing, m) && continue
            any(isnothing, m) && throw(ArgumentError("Some chain modification does not have position, \"$mod\""))
            return [parse(UInt8, first(x.captures)) => v for x in m]
        end
        mm = match(r"^[;,]\d+(C?[ON]?)\((.*)\)$", first(ms))
        if isnothing(mm)
            mm = match(r"^[;,]\d+(C?[ON]?)(.*)$", first(ms))
        end
        isnothing(mm) && throw(ArgumentError("Invalid chain modification, \"$mod\""))
        c = parse_tailsubstituent(mm...; onlinked)
        return [parse(UInt8, first(x.captures)) => c for x in match.(r"^[;,](\d+)C?[ON]?", ms)]
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

function sort_chainmodification!(mod)  
    if first(first(mod)) isa Number
        sort!(mod; by = x -> (sub_abbr(last(x)), first(x)))
    else
        sort!(mod; by = x -> (sub_abbr(first(x)), last(x)))
    end
end

