include("class_struct.jl")
function parse_lipid(s::AbstractString)
    any(==(s), SPECIFIC_LIPID) && return parse_specific(s)
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

include("parse_group.jl")
include("parse_carbonchain.jl")
include("parse_chainmodification.jl")

"""
    parse_specific(s)

Parse specific lipid
"""
function parse_specific(s)
    throw(ArgumentError("`parse_specific` not implemented for \"$s\""))
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
