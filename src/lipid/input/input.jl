include("class_struct.jl")

parse_chemical(::Type{T}, s; kwargs...) where {T <: AbstractLipid} = parse_lipid(s; lipidtype = T, kwargs...)

function parse_lipid(s::AbstractString; lipidtype = AbstractLipid, kwargs...)
    if lipidtype <: LipidChemical 
        return LipidChemical(parse_chemical(Chemical, s; kwargs...))
    end
    try 
        if any(==(s), SPECIFIC_LIPID) 
            lipid =  parse_specific(s)
            lipid isa lipidtype || throw(ArgumentError("The input string cannot be parsed as `$T`"))
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
        Con <: lipidtype || throw(ArgumentError("The input string cannot be parsed as `$lipidtype`"))
        chain, sn = parse_carbonchain(Con, bone, echain, schain)
        make_lipid(Con, bone, pos, chain, sn)
    catch e 
        if lipidtype <: AbstractLipid
            throw(e)
            @warn "Construct generic `LipidChemical` instead" 
            LipidChemical(parse_chemical(Chemical, s; kwargs...))
        else
            throw(e)
        end
    end
end
# match(r"([\s,/,_](?:\[[RS]\])*[d,t,e]?[P,O,E]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)
# cbdb, pos, fg = match(r"[/,_]?((?:\[[RS]\])*[d,t,e]?[P,O,E]?-?\d*:\d*)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)

include("parse_group.jl")
include("parse_carbonchain.jl")
include("parse_chainmodification.jl")
include("check_configuration.jl")

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

"""
    make_lipid(Con, bone, pos, chain, sn)

Construct lipids from parsed information
"""
function make_lipid(Con::Type{<: FattyAcyl}, bone, pos, chain, sn)
    Con(bone, check_configuration!(bone, first(chain)))
end

function make_lipid(Con::Type{<: NacylAmine}, bone, pos, chain, sn)
    if bone isa Tuple
        if length(chain) == 1
            hbone = parse_lipid("FN 0:0")
            Con(hbone, check_configuration!(hbone, first(chain)))
        else length(chain) == nchainposition(Con)        
            bhead, bbone, bechain = bone
            hbone = make_lipid(bhead, bbone, nothing, [first(chain)], 0x00)
            Con(hbone, check_configuration!(hbone, last(chain)))
        end
    else
        Con(bone, check_configuration!(bone, first(chain)))
    end
end

function make_lipid(Con::Type{<: FattyAcylEster}, bone, pos, chain, sn)
    bhead, bbone, bechain = bone
    if length(chain) == 1
        if bechain == (Acyl, )
            hbone = parse_lipid("FA 0:0")
        else
            hbone = parse_lipid("FOH 0:0")
        end
        Con(hbone, check_configuration!(hbone, first(chain)), sn)
    elseif length(chain) == nchainposition(Con)     
        hbone = make_lipid(bhead, bbone, nothing, [first(chain)], 0x00)
        Con(hbone, check_configuration!(hbone, last(chain)), sn)
    else
        throw(ArgumentError("Invalid lipid name"))
    end
end

# Radylglycerophosphoglycerol Single -> Tuple
function make_lipid(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, pos, chain, sn)
    ch = class_rs(bone)
    if length(ch) == 1
        ch = first(ch)
    else
        ch = (reverse!(ch)..., )
    end
    for c in tuplize(chain)
        check_configuration!.(Ref(bone), c)
        if ispropertyposition(c.doublebond) && any(x -> first(x) == 0x00, c.doublebond)
            throw(ArgumentError("Position `1` cannot contain double bond."))
        end
    end
    Con(bone, chain, sn, ch)
end

function make_lipid(Con::Type{<: GlycerophosphoNacylethanolamine}, bone, pos, chain, sn)
    ch = class_rs(bone)
    if length(ch) == 1
        ch = first(ch)
    else
        ch = (reverse!(ch)..., )
    end
    Con(bone, check_configuration!.(Ref(bone), chain), ch)
end

function parse_headposition(x)
    p, c = x 
    parse(UInt8, p) => isempty(c) ? RSChirality() : c == "R" ? RChirality() : c == "S" ? SChirality() : throw(ArgumentError("Invalid chirality."))
end

function make_lipid(Con::Type{<: SphingoBone}, bone, pos, chain, sn)
    bone_species = false
    if isnothing(pos) || isnothing(bone) || bone isa GlyComp
        poss = Pair{UInt8, RSSystem}[0x00 => RSChirality()]
        bone_species = true
    else
        poss = Pair{UInt8, RSSystem}[parse_headposition(x) for x in eachmatch(r"(\d+)([RS]*)", pos)]
    end
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if bone_species && chain_species || isnothing(bone)
        check_configuration!(bone, fc) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    elseif chain_species
        all(iszero ∘ first, poss) || throw(ArgumentError("Headgroup position cannot be specified at (molecular)specieslevel or without headgroup(s)"))
        check_configuration!(bone, fc) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    elseif ispropertyposition(fc.substituent) || (isnothing(fc.substituent) && any(!iszero ∘ first, poss))
        if isnothing(fc.substituent)
            sub = Pair{UInt8, AbstractFunctionalGroup}[]
        elseif XLinkedFunctionalGroup <: eltype(fc.substituent).parameters[end]
            sub = fc.substituent
        else
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
        end
        for p in poss
            push!(sub, first(p) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone)))
        end
        sort_chainmodification!(sub) 
        if (0x02 => Hydrogen()) in sub 
            ch = vcat(fc.chirality, poss) 
        else
            push!(sub, 0x02 => Hydrogen())
            ch = vcat(fc.chirality, poss, [0x02 => RSChirality()]) 
        end
        # check poss chirality
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, ch, fc.isotopiclabel)
        check_configuration!(bone, fc) 
        chi = [last(fc.chirality[findfirst(x -> first(x) == p, fc.chirality)]) for (p, _) in poss]
        chi = length(chi) == 1 ? first(chi) : (chi..., )
    elseif ispropertynumber(fc.substituent) || isnothing(fc.substituent)
        if isnothing(fc.substituent)
            sub = Pair{XLinkedFunctionalGroup, UInt8}[]
        else
            sub = convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, fc.substituent)
        end
        m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone))
        id = findfirst(x -> first(x) == m, fc.substituent)
        if isnothing(id)
            push!(sub, m => UInt8(length(poss)))
        else
            n = last(sub[id])
            sub[id] = m => (n + UInt8(length(poss)))
        end
        sort_chainmodification!(sub)
        headposition = [first(p) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone)) for p in poss]
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, poss, fc.isotopiclabel)
        check_configuration!(bone, fc; headposition) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    else
        throw(ArgumentError("Unknown headgroup position when \"OH\"'s positions are known"))
    end
    if length(poss) == 1
        enpos = first(first(poss))
    elseif length(poss) == 2
        enpos = first(first(poss)) * UInt8(32) + first(last(poss))
    else
        throw(ArgumentError("Invalid headgroup position, \"$pos\""))
    end
    Con(bone, length(chain) == 1 ? fc : (fc, check_configuration!.(nothing, chain[begin + 1:end]...)), enpos, chi)
end

function make_lipid(Con::Type{<: MixSphingoBone}, bone, pos, chain, sn)
    bone_species = false
    if isnothing(bone)
        poss = UInt8[]
        bone_species = true
    elseif isnothing(pos)
        poss = [[0x00 => RSChirality()] for i in eachindex(bone)]
        bone_species = true
    else
        poss = [[parse_headposition(x) for x in eachmatch(r"(\d+)([RS]*)", p)] for p in split_class_position(pos)]
    end
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if bone_species && chain_species || isnothing(bone)
        check_configuration!(bone, fc) 
        chi = last.(first.(poss))
    elseif chain_species
        all(iszero ∘ first, vcat(poss...)) || throw(ArgumentError("Headgroup position cannot be specified at specieslevel or without headgroup(s)"))
        check_configuration!(bone, fc) 
        chi = last.(first.(poss))
    elseif ispropertyposition(fc.substituent) || (isnothing(fc.substituent) && any(x -> any(!iszero ∘ first, x), poss))
        if isnothing(fc.substituent)
            sub = Pair{UInt8, AbstractFunctionalGroup}[]
        else
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
        end
        length(poss) == length(bone) || throw(ArgumentError("Number of headgroups does not equal number of positions."))
        for (b, p) in zip(bone, poss)
            m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b))
            for x in p
                push!(sub, first(x) => m)
            end
        end
        push!(sub, 0x02 => Hydrogen())
        sort_chainmodification!(sub)
        ch = vcat(fc.chirality, poss..., [0x02 => RSChirality()]) 
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, ch, fc.isotopiclabel)
        check_configuration!(bone, fc) 
        chi = map(poss) do pos 
            if length(pos) == 1
                last(fc.chirality[findfirst(x -> first(x) == first(first(pos)), fc.chirality)])
            else 
                [last(fc.chirality[findfirst(x -> first(x) == first(p), fc.chirality)]) for p in pos]
            end
        end
    elseif ispropertynumber(fc.substituent) || isnothing(fc.substituent)
        if isnothing(fc.substituent)
            sub = Pair{AbstractFunctionalGroup, UInt8}[]
        else
            sub = fc.substituent
        end
        length(poss) == length(bone) || throw(ArgumentError("Number of headgroups does not equal number of positions."))
        for (b, p) in zip(bone, poss)
            m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b))
            id = findfirst(x -> first(x) == m, fc.substituent)
            if isnothing(id)
                push!(sub, m => UInt8(length(p)))
            else
                n = last(sub[id])
                sub[id] = m => (n + UInt8(length(p)))
            end
        end
        sort_chainmodification!(sub)
        headposition = mapreduce(vcat, zip(bone, poss)) do (b, p) 
            m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b))
            [first(pp) => m for pp in p]
        end
        fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, vcat(poss...), fc.isotopiclabel)
        check_configuration!(bone, fc; headposition) 
        chi = last.(first.(poss))
    else
        throw(ArgumentError("Unknown head group position when \"OH\"'s positions are known"))
    end
    if any(x -> isa(x, Vector), poss)
        enpos = map(poss) do p 
            p isa Int ? UInt8(p) : length(p) == 1 ? UInt8(first(first(p))) : length(p) == 2 ? (UInt8(first(first(p))) * UInt8(32) + UInt8(first(last(p)))) : throw(ArgumentError("Invalid head group position, \"$p\""))
        end
    else
        enpos = convert(Vector{UInt8}, pos)
    end
    if eltype(chi) <: RSSystem
        chi = convert(Vector{RSSystem}, chi)
    end
    Con(bone, length(chain) == 1 ? fc : (fc, check_configuration!.(nothing, chain[begin + 1:end]...)), enpos, chi)
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

function parse_singlechirality(x)
    if isnothing(x) || isempty(x)
        RSChirality
    elseif x == "(R)-"
        RChirality
    elseif x == "(S)-"
        SChirality
    elseif x == "L-"
        LForm
    elseif x == "D-"
        DForm
    else
        @warn "Invalid chirality"
        RSChirality
    end
end
# specific for glycerol-phospho-aa
function parse_glycerol_rs(pre)
    if isnothing(pre)
        Glycerol{RSChirality}()
    elseif pre == "(R)-"
        Glycerol{RChirality}()
    elseif pre == "(S)-"
        Glycerol{SChirality}()
    else
        @warn "Invalid chirality"
        Glycerol{RSChirality}()
    end
end

function parse_glycerol2_rs(pre)
    isnothing(pre) && return (Glycerol{RSChirality}(), Glycerol{RSChirality}())
    rss = split(replace(pre, "(" => "", ")" => "", " " => "", "-" => ""), ",")
    gly1, gly2 = (Glycerol{RSChirality}(), Glycerol{RSChirality}())
    for rs in rss 
        if rs == "R"
            gly1 = Glycerol{RChirality}()
        elseif rs == "S"
            gly1 = Glycerol{SChirality}()
        elseif rs == "R'"
            gly2 = Glycerol{RChirality}()
        elseif rs == "S'"
            gly2 = Glycerol{SChirality}()
        else
            @warn "Invalid chirality"
        end
    end
    (gly2, gly1)
end

function parse_glycerol3_rs(pre)
    isnothing(pre) && return (Glycerol{RSChirality}(), Glycerol{RSChirality}(), Glycerol{RSChirality}())
    rss = split(replace(pre, "(" => "", ")" => "", " " => "", "-" => ""), ",")
    gly1, gly2, gly3 = (Glycerol{RSChirality}(), Glycerol{RSChirality}(), Glycerol{RSChirality}())
    for rs in rss 
        if rs == "R"
            gly1 = Glycerol{RChirality}()
        elseif rs == "S"
            gly1 = Glycerol{SChirality}()
        elseif rs == "R'"
            gly2 = Glycerol{RChirality}()
        elseif rs == "S'"
            gly2 = Glycerol{SChirality}()
        elseif rs == "R''"
            gly3 = Glycerol{RChirality}()
        elseif rs == "S''"
            gly3 = Glycerol{SChirality}()
        else
            @warn "Invalid chirality"
        end
    end
    (gly3, gly2, gly1)
end

function parse_aa_glycerol_rs(pre, aa = "Ser")
    isnothing(pre) && return (parse_aa(aa), Glycerol{RSChirality}())
    rs, dl = match(r"(?:\(([RS\,']+)\)-)*([DL]-)*", pre)
    if isnothing(dl)
        aaa = parse_aa(aa)
    else
        aaa = parse_aa(string(dl, aa))
    end
    if isnothing(rs)
        gly = Glycerol{RSChirality}()
    else
        rss = split(replace(rs, " " => ""), ",")
        gly = Glycerol{RSChirality}()
        for rs in rss 
            if rs == "R"
                gly = Glycerol{RChirality}()
            elseif rs == "S"
                gly = Glycerol{SChirality}()
            elseif rs == "R'"
                aaa = parse_aa(string("(R)-", aa))
            elseif rs == "S'"
                aaa = parse_aa(string("(S)-", aa))
            else
                @warn "Invalid chirality"
            end
        end
    end
    (aaa, gly)
end

function class_rs(backbone)
    s = getchaincomponent(backbone)
    rs = collect(map(glycerol_aa_rs, s))
    filter!(!isnothing, rs)
    all(==(AChirality()), rs) ? nothing : rs
end

glycerol_aa_rs(x) = nothing
glycerol_aa_rs(::Glycerol{T}) where T = T()
glycerol_aa_rs(aa::αAminoAcid) = alpha_rs(aa)

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
