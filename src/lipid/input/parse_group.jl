"""
"""
function parse_glycomp2(s::AbstractString)
    parse_glycomp(replace(s, r"^SHex" => "Hex3S"))
end

# function iter_glycomp!(result::Vector{Pair{Monosaccharide, UInt8}}, s::AbstractString)
#     @inbounds for i in Iterators.reverse(eachindex(s))
#         try
#             mono = parse_monoglycomp(s[begin:i])
#             if i == lastindex(s)
#                 push!(result, mono => 0x01)
#                 return result
#             else
#                 s = string(s[nextind(s, i):end])
#                 n = match(r"^\d+", s)
#                 if isnothing(n)
#                     push!(result, mono => 0x01)
#                     return iter_glycomp!(result, s)
#                 else
#                     i = n.match.offset + n.match.ncodeunits
#                     n = parse(UInt8, n.match)
#                     push!(result, mono => n)
#                     return i == lastindex(s) ? result : iter_glycomp!(result, s[nextind(s, i):end])
#                 end
#             end
#         catch
#             continue
#         end
#     end
#     throw(ArgumentError("Invalid monosaccharide(s), \"$s\""))
# end

# function parse_monoglycomp(s::AbstractString)
#     try
#         parse_monosaccharide(s)
#     catch
#         s == "SHex" && return Hexose([Sulfo() => 0x01])
#         throw(ArgumentError("Invalid monosaccharide, \"$s\""))
#     end
# end

function parse_gsl_isomer(cls, post)
    parse_series_glycan(string(cls, isnothing(post) ? "" : post))
end

"""
    parse_headgroup(s)

Parse string into chemical as headgroup
"""
function parse_headgroup(s::AbstractString)
    s = string(s)
    m = match(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])+$", s)
    if !isnothing(m)
        ms = [replace(m.match, r"^\[" => "", r"\]" => "") for m in eachmatch(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])", m.match)]
        pushfirst!(ms, s[begin:m.match.offset])
        return ntuple(i -> parse_headgroup(ms[i]), length(ms))
    end
    p = firstindex(s)
    cs = AbstractChemical[]
    ls = Pair{AbstractLinkageposition, AbstractLinkageposition}[]
    linktocheck = nothing
    trace = false
    @inbounds for m in eachmatch(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])", s)
        n = m.match.offset
        if n > p
            rs = collect(eachmatch(r"((?:[DL]+-)*[^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:n]))
            ms = get_headgroup_compound.(rs)
            dmp = [dehydroxyposition(a) => dehydrogenposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
            push!(dmp, dehydroxyposition(last(ms)) => 0x00)
            trace = false
            for (x, c, d) in zip(rs, ms, dmp)
                trace = push_new_headgroup!(cs, ls, x, c, d)
            end
            linktocheck = trace ? lastindex(ls) : nothing
        end
        push!(cs, parse_headgroup(first(match(r"^\[?(.*?)\]?$", m.match))))
        push!(ls, last(last(cs).linkage))
        p = n + m.match.ncodeunits + 1
    end
    li = p
    rs = collect(eachmatch(r"((?:[DL]+-)*[^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:end]))
    if !isempty(rs) && first(last(rs)) in ["L", "D"]
        rs = rs[begin:end - 1]
    end
    if !isempty(rs)
        ms = get_headgroup_compound.(rs)
        if !isnothing(linktocheck)
            dp = dehydrogenposition(first(ms))
            ls[linktocheck] = first(ls[linktocheck]) => Linkageposition(dp)
        end 
        dmp = Pair[dehydroxyposition(a) => dehydrogenposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
        push!(dmp, dehydroxyposition(last(ms)) => 0x00)
        trace = false
        for (x, c, d) in zip(rs, ms, dmp)
            trace = push_new_headgroup!(cs, ls, x, c, d)
            li = p + x.match.offset + x.match.ncodeunits
        end
    end
    linktocheck = trace ? lastindex(ls) : nothing
    if li < lastindex(s)
        push!(cs, parse_singleheadgroup(first(match(r"^\[?(.*?)\]?$", s[li:end]))))
        if !isnothing(linktocheck)
            ls[linktocheck] = first(ls[linktocheck]) => Linkageposition(dehydrogenposition(last(cs)))
        end 
    end
    if length(cs) == 1 && (isempty(ls) || first(first(ls)) in [Linkageposition(0x00), Linkageposition(nothing)])
        return first(cs)
    end
    all(x -> x isa Saccharide, cs) ? Glycan((cs..., ), convert(Vector{Pair{AbstractAnomerposition, Linkageposition}}, ls)) :  
    all(x -> x isa αAminoAcid, cs) ? Peptide((cs..., )) : DehydratedChemical((cs..., ), ls, nothing)
end
parse_headgroup(::Nothing) = nothing

function parse_singleheadgroup(s)
    s = startswith(s, r"\d+:\d+") ? string("FOH ", s) : s
    try 
        parse_lipid(s)
    catch
        try
            parse_monosaccharide(s)
        catch
            try
                parse_aa3(s)
            catch
                parse_lipidonlygroup(s)
            end
        end
    end
end

function parse_lipidonlygroup(s)
    if s == "SHex"
        Hexose([Sulfo() => 0x01])
    elseif s == "PI"
        makechemical(DehydratedChemical, Inositol(), PhosphoricAcid())
    elseif s == "PE"
        makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid())
    elseif s == "P"
        PhosphoricAcid()
    else
        throw(ArgumentError("Invalid abbreviation, \"$s\"."))
    end
end

function get_headgroup_compound(x)
    m, mp, mi, mj = x
    parse_singleheadgroup(first(match(r"^\[?(.*?)\]?$", m)))
end

function push_new_headgroup!(cs, ls, x, c, d)
    m, mp, mi, mj = x
    push!(cs, c)
    # deal with DehydratedChemical
    i = isempty(mi) ? first(d) : parse(UInt8, mi)
    lf = (isempty(mp) || mp == "(") ? (c isa Monosaccharide ? Anomerposition(i) : Linkageposition(i)) : (mp == "α" || mp == "(a") ? Alphaposition(i) : (mp == "β" || mp == "(b") ? Betaposition(i) : throw(ArgumentError("Invalid linkage, \"$(string(x))\""))
    rt = isempty(mj) ? Linkageposition(last(d)) : Linkageposition(parse(UInt8, mj))
    push!(ls, lf => rt)
    isempty(mj)
end

function split_tail_components(s)
    r = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?", s))
    i = firstindex(s)
    v = Tuple[]
    k = firstindex(s) - 1
    for m in r 
        j = prevind(s, m.match.offset + 1)
        ss = s[i:j]
        if isnothing(match(r"(?:\(?(\d*)-(\d*)([abαβ]?)\)?)?[DL]$", ss)) && isnothing(match(r"(>:\(?(\d*)-(\d*)([abαβ]?)\)?)?\([\dRS\,]*\)$", ss))
            ss = s[nextind(s, k):j]
            if !isempty(ss)
                if k < firstindex(s) 
                    push!(v, ("", "", "", ss)) 
                else
                    push!(v, (match(r"\(?(\d*)-(\d*)([abαβ]?)\)?", s[i:k]).captures..., ss))
                end
            end
            i = nextind(s, j)
            k = prevind(s, m.match.offset + m.match.ncodeunits + 1)
        end
    end
    m = match(r"\(?(\d*)-(\d*)([abαβ]?)\)?", s[i:k])
    push!(v, isnothing(m) ? ("", "", "", s[nextind(s, k):end]) : (m.captures..., s[nextind(s, k):end]))
    v
end

"""
    parse_tailgroup(s)

Parse string into chemical as tailgroup
"""
function parse_tailgroup(s::AbstractString)
    s = string(s)
    # linkage
    p = firstindex(s)
    cs = AbstractChemical[]
    ls = Pair{AbstractLinkageposition, AbstractLinkageposition}[]
    m = match(r"(\[\(?\d*-[^\]\[]*(?:(?1)[^\]\[]*)*\])", s)
    prep = 0x00
    if isnothing(m) # no branch
        # rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?((?:[DL]+-)*[^-\[]*[^\[])", s))
        rs = split_tail_components(s)
        # if isempty(rs) # single compound
        #     return (first(match(r"^\[?(.*?)\]?$", s)))
        # if first(rs)[begin:end - 1] == ("", "", "") # start with compound
        #     pushfirst!(cs, parse_singletailgroup(first(rs)[end]))
        #     prep = dehydrogenposition(first(cs))
        # else # start with -
        #     prep = 0x00
        # end
        ms = get_tailgroup_compound.(rs)
        dmp = Pair[dehydrogenposition(a) => dehydroxyposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
        pushfirst!(dmp, prep => dehydroxyposition(first(ms)))
        for (x, c, d) in zip(rs, ms, dmp)
            push_new_tailgroup!(cs, ls, x, c, d)
        end
        if length(ms) == 1 && first(first(ls)) in [Linkageposition(0x00), Linkageposition(nothing)]
            return first(ms)
        end
        p = lastindex(s)
    # elseif m.match.offset + 1 > firstindex(s) # start with compound
    #     p = m.match.offset + 1
    #     s2 = s[begin:m.match.offset]
    #     # m = match(r"\(?(\d*)-(\d*)([abαβ]?)\)?((?:[DL]+-)*[^-\[]*[^\[])", s2)
    #     rs = split_tail_components(s2)
    #     ms = get_tailgroup_compound.(rs)
    #     if isnothing(m) # single compound
    #         pushfirst!(cs, parse_singletailgroup(s2))
    #         prep = dehydrogenposition(first(cs))
    #     elseif m.match.offset + 1 > firstindex(s) # start with compound
    #         pushfirst!(cs, parse_singletailgroup(s2[begin:m.match.offset]))
    #         prep = dehydrogenposition(first(cs))
    #         p = m.match.offset + 1
    #     else # start with -
    #         prep = 0x00
    #         p = firstindex(s)
    #     end
    # else # start with -
    #     prep = 0x00
    end
    @inbounds for m in eachmatch(r"(\[\(?\d*-[^\]\[]*(?:(?1)[^\]\[]*)*\])", s)
        n = m.match.offset
        if n > p
            # rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?((?:[DL]+-)*[^-\[]*[^\[])", s[p:n]))
            rs = split_tail_components(s[p:n])
            ms = get_tailgroup_compound.(rs)
            dmp = [dehydrogenposition(a) => dehydroxyposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
            pushfirst!(dmp, prep => dehydroxyposition(first(ms)))
            for (x, c, d) in zip(rs, ms, dmp)
                push_new_tailgroup!(cs, ls, x, c, d)
            end
            prep = dehydrogenposition(last(ms))
        end
        pushfirst!(cs, parse_tailgroup(first(match(r"^\[?(.*?)\]?$", m.match))))
        pushfirst!(ls, last(first(cs).linkage))
        p = n + m.match.ncodeunits + 1
    end
    if p < lastindex(s) # tail main chain
        # rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?((?:[DL]+-)*[^-\[]*[^\[])", s[p:end]))
        rs = split_tail_components(s[p:end])
        ms = get_tailgroup_compound.(rs)
        dmp = [dehydrogenposition(a) => dehydroxyposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
        pushfirst!(dmp, prep => dehydroxyposition(first(ms)))
        for (x, c, d) in zip(rs, ms, dmp)
            push_new_tailgroup!(cs, ls, x, c, d)
        end
    end 
    all(x -> x isa Saccharide, cs) ? Glycan((cs..., ), convert(Vector{Pair{AbstractAnomerposition, Linkageposition}}, ls)) : 
    all(x -> x isa αAminoAcid, cs) ? Peptide((cs..., )) : DehydratedChemical((cs..., ), ls)
end
parse_tailgroup(::Nothing) = nothing

function parse_tailsubstituent(l::AbstractString, fg::AbstractString; onlinked = true)
    if isempty(l)
        # specieslevel
        # only carbon chain
        # parse_singletailgroup(fg)
        # throw(ArgumentError("Invalid species modification, \"$s\""))
        # not allowed ;5(16:1) => FOL 16:1 
        m = parse_singletailgroup(fg)
        return dehydroxygroup(m)
    end
    !onlinked && !isempty(l) && throw(ArgumentError("O-linked or N-linked modification not allowed, \"$s\""))
    # fg = startswith(fg, r"\(") ? replace(fg, r"^\(" => "", r"\)$" => "") : fg
    # fg = startswith(fg, r"\d+:\d+") ? string("FOH ", fg) : fg
    l2 = match(r"^\(?(\d*)-(\d*)([abαβ]?)\)?", string(fg))
    if !isnothing(l2)
        fg = fg[nextind(fg, l2.match.offset):end]
    end
    m = nothing
    if l == "O" || l == "N"
        for (k, v) in FG_CLINKED
            if k == fg 
                m = parentchemical(v)
                break 
            end
        end
        m = isnothing(m) ? parse_tailgroup(fg) : m
        mf = last(getchaincomponent(m))
        lg = Dehydroxy
        l = l == "O" ? OLinkage() : NLinkage()
    elseif l == "CO"
        for (k, v) in FG_CLINKED
            if k == fg 
                m = parentchemical(v)
                break 
            end
        end
        m = isnothing(m) ? parse_headgroup(fg) : m 
        mf = first(getchaincomponent(m))
        lg = Dehydrogen
        l = CarboxylicLinkage()
    else
        throw(ArgumentError("Invalid linkage, \"$l\""))
    end
    if isnothing(l2)  
        i = lg == Dehydroxy ? dehydroxyposition(mf) : lg == Dehydrogen ? dehydrogenposition(mf) : 0x00
        p = mf isa Monosaccharide ? Anomerposition(i) : Linkageposition(i)
    else
        mj, mi, mp = l2
        i = if isempty(mi) 
            lg == Dehydroxy ? dehydroxyposition(mf) : lg == Dehydrogen ? dehydrogenposition(mf) : 0x00 
        else 
            parse(UInt8, mi)
        end
        # rt = isempty(mj) ? Linkageposition(first(d)) : Linkageposition(parse(UInt8, mj))
        p = isempty(mp) ? (mf isa Monosaccharide ? Anomerposition(i) : Linkageposition(i)) : (mp == "α" || mp == "a") ? Alphaposition(i) : (mp == "β" || mp == "b") ? Betaposition(i) : throw(ArgumentError("Invalid linkage, \"$(string(l2))\""))
    end
    sub = lg == Dehydroxy ? dehydroxygroup(m; position = p) : 
            lg == Dehydrogen ? dehydrogengroup(m; position = p) : 
                lg == Odehydrogen ? dehydrogengroup(m; position = p) : 
                    lg == Ndehydrogen ? dehydrogengroup(m; position = p) : 
                        Substituent(lg, m, p)
    XLinkedFunctionalGroup(l, sub)
end

function parse_singletailgroup(s)
    s = startswith(s, r"\d+:\d+") ? string("FOH ", s) : s
    try 
        parse_lipid(s)
    catch
        try
            parse_monosaccharide(s)
        catch
            try
                parse_aa3(s)
            catch
                parse_lipidonlygroup(s)
            end
        end
    end
end

function get_tailgroup_compound(x)
    mj, mi, mp, m = x
    parse_singletailgroup(first(match(r"^\[?(.*?)\]?$", m)))
end

function push_new_tailgroup!(cs, ls, x, c, d)
    mj, mi, mp, m = x
    pushfirst!(cs, c)
    # deal with DehydratedChemical
    i = isempty(mi) ? last(d) : parse(UInt8, mi)
    lf = isempty(mp) ? (c isa Monosaccharide ? Anomerposition(i) : Linkageposition(i)) : (mp == "α" || mp == "a") ? Alphaposition(i) : (mp == "β" || mp == "b") ? Betaposition(i) : throw(ArgumentError("Invalid linkage, \"$(string(x))\""))
    rt = isempty(mj) ? Linkageposition(first(d)) : Linkageposition(parse(UInt8, mj))
    pushfirst!(ls, lf => rt)
    isempty(mj)
end
