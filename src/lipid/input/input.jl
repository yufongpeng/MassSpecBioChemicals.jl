include("class_struct.jl")

function parse_lipid(s::AbstractString; kwargs...)
    # try 
        if any(==(s), SPECIFIC_LIPID) 
            lipid =  parse_specific(s)
            # lipid isa lipidtype || throw(ArgumentError("The input string cannot be parsed as `$T`"))
            return lipid 
        end
        headclspossil, schain = split_class_chain(s)
        headcls, pos, sil = split_head_sil(headclspossil)
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
        # Con <: lipidtype || throw(ArgumentError("The input string cannot be parsed as `$lipidtype`"))
        chain, sn = parse_carbonchain(Con, bone, echain, schain)
        pos = parse_headposition(Con, bone, pos)
        make_lipid(Con, bone, pos, chain, sn)
    # catch e 
        # @warn "Construct generic `LipidChemical` instead" 
        # throw(e)
        # LipidChemical(parse_chemical(Chemical, s; kwargs...))
    # end
end
# match(r"([\s,/,_](?:\[[RS]\])*[d,t,e]?[P,O,E]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)
# cbdb, pos, fg = match(r"[/,_]?((?:\[[RS]\])*[d,t,e]?[P,O,E]?-?\d*:\d*)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)


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
    f = lastindex(s)
    dp = 0
    db = 0
    for i in eachindex(s)
        @inbounds x = s[i]
        if dp == 0 && db == 0 && x == ' '
            p = true
            e = i
            f = min(f, e)
            continue
        # elseif dp == 0 && db == 0 && p && isnothing(match(r"^(?:[dte]?[OPE]-)?\d+:\d+", s[i:end]))
        elseif dp == 0 && db == 0 && p
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
    e == 0 && return ("", string(' ', s))
    pcl = s[firstindex(s):prevind(s, f)]
    pch = s[nextind(s, e):lastindex(s)]
    ncl, nch = split_class_chain(pch)
    (isempty(ncl) ? pcl : string(pcl, " ", ncl), nch)
end

"""
"""
split_head_sil(s) = match(r"^(.*?)[\s^\[]*(\([DLRS\d\,'\s]+\))?\s*(\[[^\]\[]+\])?$", s)

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

include("parse_backbone.jl")
include("parse_group.jl")
include("parse_carbonchain.jl")
include("parse_chainmodification.jl")
include("parse_position.jl")
include("const.jl")
