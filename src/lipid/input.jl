parse_chemical(::Type{<: Lipid}, s) = parse_lipid(s)

function parse_lipid(s::AbstractString)
    any(==(s), SPECIFIC_LIPID) && return parse_spesific(s)
    headclspossil, schain = split_class_chain(s)
    headcls, pos, sil = match(r"^(.*?)[\s^\[]*(\([\d\,'\s]+\))?\s*(\[[^\]\[]+\])?$", headclspossil)
    result = nothing
    cls = ""
    for (c, o) in REGEX[:class]
        result = match(c, string(headcls))
        isnothing(result) || (cls = string(o); break)
    end
    isnothing(result) && throw(ArgumentError("Invalid lipid class"))
    head, pre, post = result
    head = isnothing(head) ? nothing : isempty(head) ? nothing : head
    pre = isnothing(pre) ? nothing : isempty(pre) ? nothing : pre
    post = isnothing(post) ? nothing : isempty(post) ? nothing : post
    parse_head = get(CLASS_STRUCT, cls, nothing)
    Con, bone, echain = parse_head(head, pre, cls, post, pos, sil)
    chain, sn = parse_carbonchain(Con, bone, echain, schain)
    make_lipid(Con, bone, pos, chain, sn)
end
# match(r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)
# cbdb, pos, fg = match(r"[/,_]?([d,t,e]?[P,O]?-?\d*:\d*)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)

"""
    parse_specific(s)

Parse specific lipid
"""
function parse_specific(s)
    throw(ArgumentError("`parse_specific` not implemented for \"$s\""))
end

"""
    make_lipid(Con, bone, pos, chain, sn)

Construct lipids from parsed information
"""
function make_lipid(Con::Type{<: FattyAcyl}, bone, pos, chain, sn)
    Con(bone, first(chain))
end

function make_lipid(Con::Type{<: NacylAmine}, bone, pos, chain, sn)
    if bone isa Tuple
        if length(chain) == 1
            Con(parse_lipid("FN 0:0"), first(chain))
        else length(chain) == nchainposition(Con)        
            bhead, bbone, bechain = bone
            Con(make_lipid(bhead, bbone, nothing, [first(chain)], 0x00), last(chain))
        end
    else
        Con(bone, first(chain))
    end
end

function make_lipid(Con::Type{<: FattyAcylEster}, bone, pos, chain, sn)
    bhead, bbone, bechain = bone
    if length(chain) == 1
        if bechain == (Acyl, )
            Con(parse_lipid("FA 0:0"), first(chain), sn)
        else
            Con(parse_lipid("FOH 0:0"), first(chain), sn)
        end
    elseif length(chain) == nchainposition(Con)        
        Con(make_lipid(bhead, bbone, nothing, [first(chain)], 0x00), last(chain), sn)
    else
        throw(ArgumentError("Invalid lipid name"))
    end
end

function make_lipid(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, pos, chain, sn)
    Con(bone, chain, sn)
end

function make_lipid(Con::Type{<: GlycerophosphoNacylethanolamine}, bone, pos, chain, sn)
    Con(bone, chain)
end

function make_lipid(Con::Type{<: SphingoBone}, bone, pos, chain, sn)
    if isnothing(pos) || bone isa GlyComp
        poss = [0x00]
    else
        poss = [parse(UInt8, x.match) for x in collect(eachmatch(r"\d+", pos))]
    end
    fc = first(chain)
    lv = annotationlevel(fc; partial = true, additional = true, pass = true)
    if isnothing(bone) || specieslevel in lv || molecularspecieslevel in lv
        all(iszero, poss) || throw(ArgumentError("Headgroup position cannot be specified at (molecular)specieslevel or without headgroup(s)"))
        nothing
    elseif any(>=(structurepositionpartiallevel), lv)
        if isnothing(fc.substituent)
            sub = Pair{UInt8, Hydroxy}[]
        elseif !isa(Hydroxy(), eltype(fc.substituent).parameters[end])
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
        else
            sub = fc.substituent
        end
        for i in poss
            push!(sub, UInt(i) => Hydroxy())
        end
        sort_chainmodification!(sub)
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, fc.isotopiclabel)
    elseif any(>=(structuredefinedpartiallevel), lv) # add o
        if isnothing(fc.substituent)
            sub = Pair{Hydroxy, UInt8}[]
        elseif !isa(Hydroxy(), eltype(fc.substituent).parameters[begin])
            sub = convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, fc.substituent)
        else
            sub = fc.substituent
        end
        for i in eachindex(poss)
            id = findfirst(x -> first(x) == Hydroxy(), fc.substituent)
            if isnothing(id)
                push!(sub, Hydroxy() => 0x01)
                sort_chainmodification!(sub)
            else
                n = last(sub[id])
                sub[id] = Hydroxy() => (n + 0x01)
            end
        end
        sort_chainmodification!(sub)
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, fc.isotopiclabel)
    elseif any(in(lv), [snpositionlevel, passsnpositionlevel, phosphatepositionlevel, passphosphatepositionlevel])
        all(iszero, poss) || throw(ArgumentError("Headgroup position cannot be specified at specieslevel or without headgroup(s)"))
        nothing
    else
        throw(ArgumentError("Unknown headgroup position when \"OH\"'s positions are known"))
    end
    if length(poss) == 1
        pos = first(poss)
    elseif length(poss) == 2
        pos = first(poss) * UInt8(32) + last(poss)
    else
        throw(ArgumentError("Invalid headgroup position, \"$pos\""))
    end
    Con(bone, length(chain) == 1 ? fc : (fc, chain[begin + 1:end]...), pos)
end

function make_lipid(Con::Type{<: MixSphingoBone}, bone, pos, chain, sn)
    pos = isnothing(bone) ? UInt8[] : isnothing(pos) ? zeros(UInt8, length(bone)) : split_class_position(pos)
    poss = vcat(collect.(pos)...)
    fc = first(chain)
    lv = annotationlevel(fc; partial = true, additional = true, pass = true)
    if isnothing(bone) || specieslevel in lv || molecularspecieslevel in lv
        all(iszero, poss) || throw(ArgumentError("Headgroup position cannot be specified at specieslevel or without headgroup(s)"))
        nothing
    elseif any(>=(structurepositionpartiallevel), lv)
        if isnothing(fc.substituent)
            sub = Pair{UInt8, Hydroxy}[]
        elseif !isa(Hydroxy(), eltype(sub).parameters[end])
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
        else
            sub = fc.substituent
        end
        for i in poss
            push!(sub, UInt(i) => Hydroxy())
        end
        sort_chainmodification!(sub)
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, fc.isotopiclabel)
    elseif any(>=(structuredefinedpartiallevel), lv) 
        if isnothing(fc.substituent)
            sub = Pair{Hydroxy, UInt8}[]
        elseif !isa(Hydroxy(), eltype(sub).parameters[begin])
            sub = convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, fc.substituent)
        else
            sub = fc.substituent
        end
        for i in eachindex(poss)
            id = findfirst(x -> first(x) == Hydroxy(), fc.substituent)
            if isnothing(id)
                push!(sub, Hydroxy() => 0x01)
                sort_chainmodification!(sub)
            else
                n = last(sub[id])
                sub[id] = Hydroxy() => (n + 0x01)
            end
        end
        sort_chainmodification!(sub)
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, fc.isotopiclabel)
    elseif any(in(lv), [snpositionlevel, passsnpositionlevel, phosphatepositionlevel, passphosphatepositionlevel])
        all(iszero, poss) || throw(ArgumentError("Headgroup position cannot be specified at specieslevel or without headgroup(s)"))
        nothing
    else
        throw(ArgumentError("Unknown head group position when \"OH\"'s positions are known"))
    end
    if any(x -> isa(x, Tuple), pos)
        pos = map(pos) do p 
            p isa Int ? UInt8(p) : length(p) == 2 ? (UInt8(first(p)) * UInt8(32) + UInt8(last(p))) : throw(ArgumentError("Invalid head group position, \"$p\""))
        end
    else
        pos = convert(Vector{UInt8}, pos)
    end
    Con(bone, length(chain) == 1 ? fc : (fc, chain[begin + 1:end]...), pos)
end

function make_lipid(Con::Type{<: Sterol}, bone, pos, chain, sn)
    throw(ArgumentError("`make_lipid` not implemented for `Sterol`"))
end

function make_lipid(Con::Type{<: Prenol}, bone, pos, chain, sn)
    throw(ArgumentError("`make_lipid` not implemented for `Prenol`"))
end

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
            push!(v, eval(Meta.parse(pos[s:e])))
            s = nextind(pos, i)
            continue
        elseif x == ')'
            dp -= 1
        end
        e = i
    end
    v
end

"""
    split_class_chain(s)
Split string into class and chain
"""
function split_class_chain(s)
    p = false
    e = 0
    dp = 0
    db = 0
    for i in eachindex(s)
        @inbounds x = s[i]
        if dp == 0 && db == 0 && x == ' '
            p = true
            e = i
            continue
        elseif dp == 0 && db == 0 && p && isnothing(match(r"(?:[dte]?[OP]-)?\d+:\d+", s[i:end]))
            break
        end
        p = false
        @match x begin
            '(' => (dp += 1)
            '[' => (db += 1)
            ')' => (dp -= 1)
            ']' => (db -= 1)
            _   => () 
        end
    end
    s[firstindex(s):e - 1], s[e:lastindex(s)]
end

"""
"""
function parse_glycomp(s::AbstractString)
    result = Pair{Monosaccharide, UInt8}[]
    GlyComp(sort!(iter_glycomp!(result, s); by = repr ∘ first))
end

function iter_glycomp!(result::Vector{Pair{Monosaccharide, UInt8}}, s::AbstractString)
    @inbounds for i in Iterators.reverse(eachindex(s))
        try
            mono = parse_monoglycomp(s[begin:i])
            if i == lastindex(s)
                push!(result, mono => 0x01)
                return result
            else
                s = string(s[nextind(s, i):end])
                n = match(r"^\d+", s)
                if isnothing(n)
                    push!(result, mono => 0x01)
                    return iter_glycomp!(result, s)
                else
                    i = n.match.offset + n.match.ncodeunits
                    n = parse(UInt8, n.match)
                    push!(result, mono => n)
                    return i == lastindex(s) ? result : iter_glycomp!(result, s[nextind(s, i):end])
                end
            end
        catch
            continue
        end
    end
    throw(ArgumentError("Invalid monosaccharide(s), \"$s\""))
end

function parse_monoglycomp(s::AbstractString)
    try
        parse_monosaccharide(s)
    catch
        s == "SHex" && return Hex([Sulfate() => 0x01])
        throw(ArgumentError("Invalid monosaccharide, \"$s\""))
    end
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
    @inbounds for m in eachmatch(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])", s)
        n = m.match.offset
        if n > p
            rs = collect(eachmatch(r"([^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:n]))
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
    rs = collect(eachmatch(r"([^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:end]))
    ms = get_headgroup_compound.(rs)
    if !isnothing(linktocheck)
        dp = dehydrogenposition(first(ms))
        ls[linktocheck] = first(ls[linktocheck]) => Linkageposition(dp)
    end 
    dmp = [dehydroxyposition(a) => dehydrogenposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
    push!(dmp, dehydroxyposition(last(ms)) => 0x00)
    trace = false
    for (x, c, d) in zip(rs, ms, dmp)
        trace = push_new_headgroup!(cs, ls, x, c, d)
        li = p + x.match.offset + x.match.ncodeunits
    end
    linktocheck = trace ? lastindex(ls) : nothing
    if li < lastindex(s)
        push!(cs, parse_singleheadgroup(first(match(r"^\[?(.*?)\]?$", s[li:end]))))
        if !isnothing(linktocheck)
            ls[linktocheck] = first(ls[linktocheck]) => dehydrogenposition(last(cs))
        end 
    end
    all(x -> x isa Saccharide, cs) ? Glycan((cs..., ), convert(Vector{Pair{AbstractAnomerposition, Linkageposition}}, ls)) :  
    all(x -> x isa αAminoAcid) ? Peptide((cs..., )) : DehydratedChemical((cs..., ), ls)
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
        Hex([Sulfate() => 0x01])
    elseif s == "PI"
        makemolecule(DehydratedChemical, Ino(), PhosphoricAcid())
    elseif s == "PE"
        makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid())
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
    if isnothing(m) # no branch
        rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?([^-\[]*[^\[])", s))
        if isempty(rs) # single compound
            return parse_singletailgroup(first(match(r"^\[?(.*?)\]?$", s)))
        elseif first(rs).match.offset + 1 > firstindex(s) # start with compound
            pushfirst!(cs, parse_singletailgroup(s[begin:first(rs).match.offset]))
            prep = dehydrogenposition(first(cs))
            p = first(rs).match.offset + 1
        else # start with -
            prep = 0x00
        end
        ms = get_tailgroup_compound.(rs)
        dmp = [dehydrogenposition(a) => dehydroxyposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
        pushfirst!(dmp, prep => dehydroxyposition(first(ms)))
        for (x, c, d) in zip(rs, ms, dmp)
            push_new_tailgroup!(cs, ls, x, c, d)
        end
        p = last(rs).match.offset + last(rs).match.ncodeunits + 1
    elseif m.match.offset + 1 > firstindex(s) # start with compound
        p = m.match.offset + 1
        s2 = s[begin:m.match.offset]
        m = match(r"\(?(\d*)-(\d*)([abαβ]?)\)?([^-\[]*[^\[])", s2)
        if isnothing(m) # single compound
            pushfirst!(cs, parse_singletailgroup(s2))
            prep = dehydrogenposition(first(cs))
        elseif m.match.offset + 1 > firstindex(s) # start with compound
            pushfirst!(cs, parse_singletailgroup(s2[begin:m.match.offset]))
            prep = dehydrogenposition(first(cs))
            p = m.match.offset + 1
        else # start with -
            prep = 0x00
            p = firstindex(s)
        end
    else # start with -
        prep = 0x00
    end
    @inbounds for m in eachmatch(r"(\[\(?\d*-[^\]\[]*(?:(?1)[^\]\[]*)*\])", s)
        n = m.match.offset
        if n > p
            rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?([^-\[]*[^\[])", s[p:n]))
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
    li = p
    if li < lastindex(s) # tail main chain
        rs = collect(eachmatch(r"\(?(\d*)-(\d*)([abαβ]?)\)?([^-\[]*[^\[])", s[p:end]))
        ms = get_tailgroup_compound.(rs)
        dmp = [dehydrogenposition(a) => dehydroxyposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
        pushfirst!(dmp, prep => dehydroxyposition(first(ms)))
        for (x, c, d) in zip(rs, ms, dmp)
            push_new_tailgroup!(cs, ls, x, c, d)
        end
    end 
    all(x -> x isa Saccharide, cs) ? Glycan((cs..., ), convert(Vector{Pair{AbstractAnomerposition, Linkageposition}}, ls)) : 
    all(x -> x isa αAminoAcid) ? Peptide((cs..., )) : DehydratedChemical((cs..., ), ls)
end
parse_tailgroup(::Nothing) = nothing

function parses_tailsubstituent(l::AbstractString, fg::AbstractString; onlinked = true)
    if isempty(l)
        # specieslevel
        # only carbon chain
        # parse_singletailgroup(fg)
        # throw(ArgumentError("Invalid species modification, \"$s\""))
        m = parse_singletailgroup(fg)
        return Substituent(Dehydroxy, m, lk(dehydroxyposition(m)))
    end
    !onlinked && !isempty(l) && throw(ArgumentError("O-linked or N-linked modification not allowed, \"$s\""))
    # fg = startswith(fg, r"\(") ? replace(fg, r"^\(" => "", r"\)$" => "") : fg
    # fg = startswith(fg, r"\d+:\d+") ? string("FOH ", fg) : fg
    l2 = match(r"^\(?(\d*)-(\d*)([abαβ]?)\)?", fg)
    if !isnothing(l2)
        fg = fg[l2.match.offset + 1:end]
    end
    if l == "O" || l == "N"
        m = parse_tailgroup(fg)
        mf = last(getchaincomponent(m))
        lg = Dehydroxy
        l = l == "O" ? OLinkage() : NLinkage()
    elseif l == "CO"
        m = parse_headgroup(fg)
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
        mj, mi, mp, _ = l2
        i = if isempty(mi) 
            lg == Dehydroxy ? dehydroxyposition(mf) : lg == Dehydrogen ? dehydrogenposition(mf) : 0x00 
        else 
            parse(UInt8, mi)
        end
        # rt = isempty(mj) ? Linkageposition(first(d)) : Linkageposition(parse(UInt8, mj))
        p = isempty(mp) ? (mf isa Monosaccharide ? Anomerposition(i) : Linkageposition(i)) : (mp == "α" || mp == "a") ? Alphaposition(i) : (mp == "β" || mp == "b") ? Betaposition(i) : throw(ArgumentError("Invalid linkage, \"$(string(l2))\""))
    end
    XLinkedFunctionalGroup(l, Substituent(lg, parse_tailgroup(fg), p))
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

function parse_gsl_isomer(cls, post)
    if endswith(cls, "?")
        iso = [replace(x, "alpha" => "α") for x in split(post, r"\s*,\s*")]
        sort!(iso)
        c = cls[begin:end - 1]
        eval(Meta.parse(cls))(ntuple(i -> eval((string(c, iso[i])))(), length(iso)))
    else
        eval(Meta.parse(replace(cls, "alpha" => "α")))() 
    end
end

"""
    parse_carbonchain(Con, bone, echain, schain)

Parse carbon chain string into (`CarbonChain`, snposition code) or ((`CarbonChain`..., ), snposition code)
"""
function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcyl}
    mchain = split_carbonchain(schain)
    length(mchain) > 1 && throw(ArgumentError("Maximal number of chains is 1, got $(length(mchain))"))
    [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: NacylAmine}
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 # sum or single
        [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
    elseif length(mchain) == nchainposition(Con) # NA xx:x/xx:x
        map(echain, mchain) do e, m
            p = parse_position_carbonchain(e, m; force = true)
            first(p) == :start || first(p) == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            last(p)
        end, 0x00
    else
        throw(ArgumentError("Maximal number of chains is $(nchainposition(Con)), got $(length(mchain))"))
    end
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcylEster}
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 
        [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
    elseif echain == (Acyl, Acyl)
        fa = last(parse_position_carbonchain(Acyl, first(mchain); force = true))
        m = match(r"/(\d+)O\(FA[^\s]*(.*)\)", last(mchain))
        if isnothing(m)
            foh = parse_position_carbonchain(Acyl, last(mchain); force = true)
            first(foh) == :inorder || throw(ArgumentError("Invalid chain position, \"$(last(mchain))\""))
            [last(foh), fa], 0x00
        else
            p, a = m
            cchain, pos, ox, sil, mod, sn = parse_singlecarbonchain(a)
            isnothing(ox) || throw(ArgumentError("This fatty acid should have position of functional group, $(last(mchain))"))
            isnothing(sn) || throw(ArgumentError("This fatty acid should not have sn position, $(last(mchain))"))
            # [last(parse_position_carbonchain(Acyl, string(a, ";", p, "OH"); force = true)), fa], parse(UInt8, p)
            [last(parse_position_carbonchain(Acyl, a; force = true)), fa], parse(UInt8, p)
        end
    elseif length(mchain) == nchainposition(Con)
        map(echain, mchain) do e, m
            p = parse_position_carbonchain(e, m; force = true)
            first(p) == :start || first(p) == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            last(p)
        end, 0x00
    else
        throw(ArgumentError("Maximal number of chains is $(nchainposition(Con)), got $(length(mchain))"))
    end
end

function parse_carbonchain(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, echain, schain)
    mchain = split_carbonchain(schain)
    maxsn = nchainposition(Con)
    length(mchain) > maxsn && throw(ArgumentError("Maximal number of chains is $maxsn, got $(length(mchain))"))
    infos = [parse_singlecarbonchain(m) for m in mchain]
    pos = [parse_chainposition(Radyl, info.cchain, info.sn, nothing) for info in infos]
    isn = findall(x -> !isnothing(match(r"^\(sn.*\)$", string(x))), pos)
    if !isempty(isn)
        id = collect(eachindex(pos))
        id = vcat(setdiff!(id, isn), isn)
        infos = infos[id]
        pos = pos[id]
        mchain = mchain[id]
    end
    if length(pos) > 1 && any(==(:inorder), @view pos[begin + 1:end])
        allequal(pos[begin + 1:end]) || throw(ArgumentError("Chain separater should be all \"/\" or \"\\\""))
        pos[begin] = :inorder
    end
    position = map(enumerate(pos)) do (i, c)
        @match c begin 
            :start => 0x00
            :inorder => UInt8(i)
            :noorder => 0x00
            r"\(sn.*\)" => UInt8(findfirst(==(replace(c, "(" => "", ")" => "")), chainposition(Con)))
            _  => throw(ArgumentError("Invalid chain position, \"$(mchain[i])\""))
        end
    end
    if !isempty(isn)
        id = sortperm(position; alg = MergeSort)
        infos = infos[id]
        position = position[id]
        mchain = mchain[id]
    end
    carbonchain = if length(mchain) == maxsn
        [make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for info in infos] # ordered/unordered, ex TG
    elseif length(mchain) == length(echain) # unordered
        all(==(:inorder), pos) && throw(ArgumentError("Expected number of chains in order is $maxsn, got $(length(mchain)); change chain separater to \"_\""))
        [make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for info in infos]
    else
        # remove sn
        # sort 
        nr = length(mchain) - length(isn)
        ne = length(echain) - length(isn)
        snx = ne ÷ nr
        sns = repeat([snx], nr)
        δ = ne - nr * snx
        i = firstindex(sns) # large => small ? 
        while δ > 0
            δ -= 1
            sns[i] += 1
            i += 1
        end
        sns = vcat(sns, repeat([1], length(isn)))
        [sn == 1 ? make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) : 
                        make_carbonchain(ntuple(i -> Radyl, sn), info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for (sn, info) in zip(sns, infos)]
    end

    allunique(filter(!=(0x00), position)) || throw(ArgumentError("Overlapped chain position, \"$schain\""))
    
    # check sn position distribution (sn)
    # :start => 0x01
    # :inorder => vector ID 
    # :noorder => 0x00
    # sn => chainposition ID 
    # base = maxsn + 1
    
    allunique(filter(>(0x00), position)) || throw(ArgumentError("sn position overlapped"))
    id = findall(x -> ncarbon(x) > 0, carbonchain)
    carbonchain = carbonchain[id]
    position = position[id]
    sn = maxsn > 3 ? 0x0000 : 0x00
    bs = convert(typeof(sn), maxsn + 1)
    for p in position
        sn *= bs
        sn += p
    end
    length(carbonchain) == 1 ? first(carbonchain) : ntuple(i -> carbonchain[i], length(carbonchain)), sn
    # prev = prev * base + next
end

function parse_carbonchain(Con::Type{<: Sphingolipid}, bone, echain, schain)
    mchain = split_carbonchain(schain)
    if length(echain) == 3 # ACer...
        # take chain from bone
        if bone isa FattyAcid && ncarbon(bone.chain) > 0 # FA x:y
            echain = first(echain, 2)
        elseif bone isa Tuple && ncarbon(first(bone).chain) > 0
            echain = first(echain, 2)
        else
            throw(ArgumentError("Invalid sphingolipid backbone, $bone"))
        end
        if length(mchain) == 1 # sum level
            pchain = [parse_position_carbonchain(echain, first(mchain))]
        elseif length(mchain) == 2 # species level
            pchain = [parse_position_carbonchain(C, m) for (C, m) in zip(first(echain, 2), mchain)]
        else
            throw(ArgumentError("Invalid number of chains, $(length(mchain))"))
        end
    elseif length(mchain) == length(echain)
        pchain = [parse_position_carbonchain(C, m) for (C, m) in zip(echain, mchain)]
    else # sum level
        pchain = [parse_position_carbonchain(echain, first(mchain))]
    end
    for (i, c) in enumerate(pchain)
        @match first(c) begin 
            :start => 0x00
            :inorder => UInt8(i)
            _ => throw(ArgumentError("Invalid chain position, \"$(mchain[i])\""))
        end
    end
    ntuple(i -> last(pchain[i]), length(pchain)), 0x00
end
function parse_carbonchain(Con::Type{<: Sterol}, bone, pos, chain, sn)
    throw(ArgumentError("`parse_carbonchain` not implemented for `Sterol`"))
end

function parse_carbonchain(Con::Type{<: Prenol}, bone, pos, chain, sn)
    throw(ArgumentError("`parse_carbonchain` not implemented for `Prenol`"))
end

"""
    split_carbonchain(schain)

Split carbon chain string by "/" and "_"; return a vector of carbon chain string
"""
function split_carbonchain(schain)
    s1 = split(schain, "/")
    fs = popfirst!(s1)
    chains = String[]
    s2 = split(fs, "_")
    push!(chains, popfirst!(s2))
    for s in s2
        push!(chains, string("_", s))
    end
    for s in s1
        s2 = split(string("/", s), "_")
        push!(chains, popfirst!(s2))
        for s in s2
            push!(chains, string("_", s))
        end
    end
    chains
end

"""
    parse_position_carbonchain(T, mchain; position = nothing, force = false)

Parse carbon chain string into snposition symbol => `CarbonChain`

`force`: force carbon chain to be `T`
"""
function parse_position_carbonchain(T, mchain; position = nothing, force = false)
    # pos => chain
    # pos:
    # :start start
    # :inorder ordered
    # :noorder unordered
    # number for ACer
    # sn-\d
    # cchain_sil = match(r"^([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?", mchain)
    # isnothing(cchain_sil) && throw(ArgumentError("Invalid fattyacyl chain, \"$mchain\""))
    # m = cchain_sil.match
    # cchain, pos, sil = cchain_sil
    # mchain = mchain[m.offset + m.ncodeunits + 1:end]
    # mod = collect(eachmatch(r"((?:;(([^)(;/_]*[\(\[][^)(]*+(?:(?3)[^)(]*)*+[\)\]])?[^)(;/_]*)))", mchain))
    #                           ((?:;[^)(;/_]*([\(\[][^)(]*+(?:(?2)[^)(]*)*+[\)\]])?[^)(;/_]*)*)
    # sn, = match(r"(\(sn-*\d*'*\))?$", mchain)
    # r"((?:;(([^)(\[\];/_]*\([^)(]*+(?:(?3)[^)(]*)*+\))?([^)(\[[^\[\]]*+(?:(?3)[^\[\]]*)*+\])?[^)(;/_]*)))"
    cchain, pos, ox, sil, mod, sn = parse_singlecarbonchain(mchain)
    parse_chainposition(T, cchain, sn, position) => make_carbonchain(T, cchain, pos, ox, sil, split_chainmodification(mod); force)
end

function parse_singlecarbonchain(mchain) 
    cchain, pos, ox, sil, mod, _, sn = (isnothing(x) ? nothing : isempty(x) ? nothing : x for x in match(REGEX[:chain], mchain))
    (; cchain, pos, ox, sil, mod, sn)
end


"""
    make_carbonchain(T, cchain::AbstractString, pos, ox, sil, mod; force = false)
    make_carbonchain(T, carbon::Number, doublebond, substituent, isotopiclabel)

Construct `CarbonChain` with parsed information
"""
function make_carbonchain(T::Tuple, cchain::AbstractString, pos, ox, sil, mod; force = false)
    if length(T) == 1
        make_carbonchain(first(T), cchain, pos, ox, sil, mod; force)
    elseif SPB in T
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(SPB, cchain, pos, ox, sil, mod)...)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(T[i], cchain, pos, ox, sil, mod)...)
    elseif force
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(Radyl, cchain, pos, ox, sil, mod)...)
    else
        n, rad, chain = match(r"([d,t,e]?)([P,O]?-?)(\d+:\d+.*)", cchain)
        if isnothing(rad) 
            Chain = Tuple{ntuple(i -> Acyl, length(T))...}
        else
            n = @match n begin
                ""  => 1
                "d" => 2
                "t" => 3
                "e" => 4
                _   => throw(ArgumentError("Invalid representation of numbers of alkyl or alkenyl chain"))
            end
            echain = [Acyl for i in eachindex(T)]
            rad = @match rad begin
                ""   => Acyl
                "O-" => Alkyl
                "P-" => Alkenyl
                _    => throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
            end
            echain[begin:begin + n - 1] .= rad
            Chain = Tuple{echain...}
        end
        CarbonChain{Chain}(parse_carbonchainbody(Radyl, chain, pos, ox, sil, mod)...)
    end
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: Radyl}
    force && return CarbonChain{T}(parse_carbonchainbody(T, cchain, pos, ox, sil, mod)...)
    n, rad, chain = match(r"([d,t,e]?)([P,O]?-?)(\d+:\d+.*)", cchain)
    isempty(n) || throw(ArgumentError("This should be a single fattyacyl chain, \"$cchain\""))
    Chain = @match rad begin
        ""   => Acyl
        "O-" => Alkyl
        "P-" => Alkenyl
        _    => throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
    end
    Chain <: T || throw(ArgumentError("Fattyacyl chain does not match to class"))
    CarbonChain{Chain}(parse_carbonchainbody(Chain, chain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: SPB}
    rad, chain = match(r"([d,t,e]?[P,O]?-?)(\d+:\d+.*)", cchain)
    isempty(rad) || throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
    CarbonChain{T}(parse_carbonchainbody(T, chain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchain` not implemented for `$T`"))
end

make_carbonchain(T, carbon::C, doublebond, substituent, isotopiclabel = nothing) where {C <: Number} = CarbonChain{T}(UInt8(carbon), uint8ize(doublebond), uint8ize(substituent), isotopiclabel)

"""
    parse_carbonchainbody(T, cchain, pos, ox, sil, mod; force = false) -> (carbon, doublebond, substituent, isotopiclabel)

Parse carbon chain info strings into correct type and format for `CarbonChain` construction
"""
function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: Radyl}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    # mod
    sub = parse_chainmodification(T, mod, ox)
    isnothing(pos) && return cb, db, sub, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid double bond configuration, \"$(string(ps))\""))
        @inbounds db[i] = x + e
    end
    return cb, db, sub, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: SPB}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    sub = parse_chainmodification(T, mod, ox)
    isnothing(pos) && return cb, db, sub, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid double bond configuration, \"$(string(ps))\""))
        @inbounds db[i] = x + e
    end
    return cb, db, sub, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchainbody` not implemented for `$T`"))
end

"""
    parse_chainposition(T, mc, sn, position)

Parse chain position into snposition symbol
"""
function parse_chainposition(::Type{<: Radyl}, mc, sn, position)
    !isnothing(position) ? position : 
    !isnothing(sn) ? sn : 
    startswith(mc, r"\s") ? :start : 
    startswith(mc, "/") ? :inorder :
    startswith(mc, "_") ? :noorder : throw(ArgumentError("Invalid fattyacyl chain, \"$mc\""))
end

function parse_chainposition(::Type{<: SPB}, mc, sn, position)
    isnothing(sn) || throw(ArgumentError("SPB does not have sn position"))
    !isnothing(position) ? position : 
    startswith(mc, r"\s") ? :start : 
    startswith(mc, "/") ? :inorder :
    startswith(mc, "_") ? :noorder : throw(ArgumentError("Invalid fattyacyl chain, \"$mc\""))
end

function parse_chainposition(::Type{T}, mc, sn, position) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_chainposition` not implemented for `$T"))
end

function parse_chainposition(T::Tuple, mc, sn, position)
    if length(T) == 1
        parse_chainposition(first(T), mc, sn, position)
    elseif SPB in T
        parse_chainposition(SPB, mc, sn, position)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        parse_chainposition(T[i], mc, sn, position)
    else
        parse_chainposition(Radyl, mc, sn, position)
    end
end

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
        c = parses_tailsubstituent(mm...; onlinked)
        return [parse(UInt8, first(x.captures)) => c for x in match.(r"^[;,](\d+)C?[ON]?", ms)]
    elseif length(ms) == 1
        for (k, v) in FG_SP
            m = match(Regex(string("^;\\(?", k, "\\)?(\\d*)\$")), mod)
            isnothing(m) ? continue : return [v => (isempty(first(m)) ? 0x01 : parse(UInt8, first(m)))]
        end
        m = startswith(first(ms), ";(") ? match(r"^;\((C?[ON]?)\((.*)\)\)(\d*)$", first(ms)) : match(r"^;(C?[ON]?)\((.*)\)()$", first(ms)) 
        isnothing(m) && throw(ArgumentError("Invalid chain modification, \"$mod\""))
        l, fg, n = m
        return [parses_tailsubstituent(l, fg) => (isempty(n) ? 0x01 : parse(UInt8, n))]
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

function parse_sil(s)
    s = replace(s, r"(\d+[A-Z][a-z]*)" => s"[\1]")
    s = replace(s, r"(\d*),([\[,A-Z])" => s"\1\2")
end # TO FORMULA

function distribute_sil(s)
    [parse(Int, a) => b for (a, b) in eachmatch(r"(\d)-(.*?[A-Z][a-z]*\d*)", s)]
end
# PC[1-D5] PC[3-1-(1,1,2,2)D4] PC[3-2-D9]
# PC[(1,1,2,3,3)D5] PC[(1',1',2',2')D4] PC[(3',3,'3',4',4',4',5',5',5')D9] PC[(1')15N]
# function hascarbonenumeration
