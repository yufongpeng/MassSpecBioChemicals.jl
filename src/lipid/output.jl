
function chemicalname(lipid::T) where {T <: Glycerolipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(class_abbr(lipid), " ", join(rp, "/"))
end
chemicalname(lipid::GlycerophosphoNacylethanolamine) = string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
function chemicalname(lipid::T) where {T <: Glycerophospholipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(class_abbr(lipid), " ", join(rp, "/"))
end

# function chemicalname(lipid::T) where {T <: FattyAcyl{B, <: Tuple}} where B
#     string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
# end

function chemicalname(lipid::T) where {T <: FattyAcyl}
    string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
end

function chemicalname(lipid::T) where {T <: Union{NacylAlkylAmine, FattyAcylEster}}
    if ncarbon(lipid.backbone.chain) == 0
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    else
        string(replace(chemicalname(lipid.backbone), class_abbr(lipid.backbone) => class_abbr(lipid)), "/", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    end
end

function chemicalname(lipid::T) where {T <: FattyAcylEstolid}
    if ncarbon(lipid.backbone.chain) == 0
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    elseif iszero(lipid.position)
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", replace(chemicalname(lipid.backbone), string(class_abbr(lipid.backbone), " ") => ""))
    else
        string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", Int(lipid.position), "O(", chemicalname(lipid.backbone), ")")
    end
end

function chemicalname(lipid::T) where {T <: Sphingolipid{H, <: Tuple}} where H
    isnothing(lipid.headgroup) && return string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
    lv = annotationlevel(first(getlipidchain(lipid)); partial = true, additional = true, pass = true)
    if specieslevel in lv || molecularspecieslevel in lv
        return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    end
    position = decode_position(lipid)
    pos = any(==(0), position) ? "" : replace(string(position), " " => "", ",)" => ")")
    fc = deepcopy(first(getlipidchain(lipid)))
    del = Int[]
    position = vcat(collect.(position)...)
    if any(>=(structurepositionpartiallevel), lv)
        for i in position
            id = findfirst(==(i => Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            while id in del
                id = findnext(==(i => Hydroxy()), fc.substituent)
                isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
            push!(del, id)
        end
    elseif any(>=(structuredefinedpartiallevel), lv)
        for i in position
            id = findfirst(x -> ==(first(x) == Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            n = last(fc.substituent[id])
            if n > 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
            elseif n == 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
                push!(del, id)
            else
                throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
        end
    end
    deleteat!(fc.substituent, unique!(del))
    c = collect(getlipidchain(lipid))
    c[begin] = fc
    string(class_abbr(lipid), pos, " ", join(repr_singlechain.(c), "/"))
end

function chemicalname(lipid::T) where {T <: Sphingolipid{H, <: CarbonChain}} where H
    isnothing(lipid.headgroup) && return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    lv = annotationlevel(first(getlipidchain(lipid)); partial = true, additional = true, pass = true)
    if specieslevel in lv || molecularspecieslevel in lv
        return string(class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    end
    position = decode_position(lipid)
    pos = any(iszero, position) ? "" : replace(string(position), " " => "", ",)" => ")")
    fc = deepcopy(first(getlipidchain(lipid)))
    del = Int[]
    position = vcat(collect.(position)...)
    if any(>=(structurepositionpartiallevel), lv)
        for i in position
            id = findfirst(==(i => Hydroxy()), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            while id in del
                id = findnext(==(i => Hydroxy()), fc.substituent)
                isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
            push!(del, id)
        end
    elseif any(>=(structuredefinedpartiallevel), lv)
        for i in position
            id = findfirst(x -> first(x) == Hydroxy(), fc.substituent)
            isnothing(id) && throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            n = last(fc.substituent[id])
            if n > 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
            elseif n == 1
                fc.substituent[id] = Hydroxy() => (n - 0x01)
                push!(del, id)
            else
                throw(ArgumentError("No corresponding hydroxy group for headgroup at $(Int(i))"))
            end
        end
    end
    deleteat!(fc.substituent, unique!(del))
    string(class_abbr(lipid), pos, " ", repr_singlechain(fc))
end

"""
    decode_sn(lipid)

Decode snpositioncode into vector of snposition
"""
function decode_sn(lipid::T) where {T <: Lipid}
    n = length(getlipidchain(lipid))
    position = zeros(Int, n)
    sn = lipid.sn
    ep = n - 1
    bs = nchainposition(T) + 1
    for i in eachindex(position)
        p, sn = divrem(sn, bs ^ ep)
        position[i] = p
        ep -= 1
    end
    position
end

"""
    decode_position(lipid)

Decode headgroup position code into vector of headgroup position
"""
function decode_position(lipid::T) where {T <: SphingoBone}
    p = Int(lipid.position)
    if hascycloheadgroup(lipid)
        p1, p2 = divrem(p, 32)
        iszero(p1) ? (p2, ) : (p1, p2)
    else
        (p, )
    end
end

decode_position(lipid::T) where {T <: SphingoBone{Nothing}} = ()
"""
    repr_db(chain)
    repr_db(dbs)

Decode double bond code into double bond representation
"""
function repr_db(c::CarbonChain) 
    p = repr_db_position(c)
    isempty(p) ? string(ndoublebond(c)) : string(ndoublebond(c), "(", p, ")")
end
"""
    repr_db_position(chain)
    repr_db_position(dbs)

Decode double bond code into double bond position representation
"""
repr_db_position(chain::CarbonChain) = repr_db_position(chain.doublebond)
function repr_db_position(dbs)
    join(decode_db_position(dbs), ",")
end
"""
    decode_db_position(db)

Decode double bond code into double bond position
"""
decode_db_position(chain::CarbonChain) = decode_db_position(chain.doublebond)
function decode_db_position(dbs::Vector)
    filter!(!isnothing, [_decode_db_position(x) for x in dbs])
end
decode_db_position(dbs) = String[]
function _decode_db_position(db::UInt8)
    a, b = divrem(db, 3)
    a == 0x00 && return nothing
    string(a, b == 0 ? "" : b == 1 ? "Z" : b == 2 ? "E" : throw(ArgumentError("Invalid fattyacyl chain")))
end

"""
    repr_sub(sub)

Convert chain substituent(s) into readble representation
"""
function repr_sub(sub::UInt8)
    sub == 0 ? "" : string(";O", sub > 1 ? Int(sub) : "")
end 
repr_sub(sub::Vector{<: Pair{OxygenAtom, UInt8}}) = repr_sub(only(sub))
function repr_sub(sub::Pair{OxygenAtom, UInt8})
    last(sub) == 0 ? "" : string(";O", last(sub) > 1 ? Int(last(sub)) : "")
end

function repr_sub(sub::Vector{<: Pair{F, UInt8} where F})
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ first)
    subs = String[""]
    ns = Int[0]
    for x in sub
        next, p = x
        next = sub_abbr(next)
        if next == last(subs)
            subs[end] += Int(p)
        else
            push!(subs, next)
            push!(ns, Int(p))
        end
    end
    popfirst!(subs)
    popfirst!(ns)
    s = ""
    for (sub, n) in zip(subs, ns)
        s *= n > 1 ? (length(sub) > 1 ? string(";(", sub, ")", n) : string(";", sub, n)) : endswith(sub, r"\d") ? string(";(", sub, ")") : string(";", sub)
    end
    s
end 

function repr_sub(sub::Vector{<: Pair{UInt8, F} where F})
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ last)
    s = ""
    prev = ""
    for x in sub
        p, next = x
        next = sub_abbr(next)
        if next == prev
            s *= string(",", Int(p), next)
        else
            s *= string(";", Int(p), next)
            prev = next
        end
    end
    s
end 
repr_sub(::Nothing) = ""

"""
    class_abbr(lipid)

Class abbreviation
"""
class_abbr(x) = chemicalabbr(x)
class_abbr(::PhosphoricAcid) = "P"
class_abbr(::Hydrocarbon) = "HC"
class_abbr(::FattyAcid) = "FA"
class_abbr(::FattyAldehyde) = "FAL"
class_abbr(::FattyAlcohol) = "FOH"
class_abbr(::WaxEster) = "WE"
class_abbr(::FattyAmide) = "FAM"
class_abbr(::FattyAmine) = "FN"
class_abbr(::NacylAmine) = "NA"
class_abbr(::NacylAmine{<: includeSIL(Ethanolamine)}) = "NAE"
class_abbr(::NacylAmine{<: includeSIL(Taurine)}) = "NAT"
for x in NAAA
    T = typeof(parse_aa(x))
    @eval head_abbr(aa::$T) = letter3_abbr(aa)
    @eval class_abbr(c::NacylAmine{<: includeSIL($T)}) = string("NA", head_abbr(c.backbone))
end
class_abbr(::FattyAcylCarnitine) = "CAR"
class_abbr(::FattyAcylCoA) = "CoA"
class_abbr(::FattyAcylEstolid) = "FAHFA"

class_abbr(::Monoradylglycerol) = "MG"
class_abbr(::Diradylglycerol) = "DG"
class_abbr(::Triradylglycerol) = "TG"
function class_abbr(c::Omodifiedmonoradylglycerol)
    b = head_abbr(c.backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "MG", r)
end
function class_abbr(c::Omodifieddiradylglycerol)    
    b = head_abbr(c.backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "DG", r)
end
class_abbr(::Sulfoquinovosylmonoradylglycerol) = "SQMG"
class_abbr(::Sulfoquinovosyldiradylglycerol) = "SQDG"
class_abbr(::Monogalactosylmonoradylglycerol) = "MGMG"
class_abbr(::Monogalactosyldiradylglycerol) = "MGDG"
class_abbr(::Digalactosylmonoradylglycerol) = "DGMG"
class_abbr(::Digalactosyldiradylglycerol) = "DGDG"
class_abbr(::Glucuronosylmonoradylglycerol) = "GlcAMG"
class_abbr(::Glucuronosyldiradylglycerol) = "GlcADG"

class_abbr(::Phosphatidicacid) = "PA"
class_abbr(::Phosphatidylcholine) = "PC"
class_abbr(::Phosphatidylethanolamine) = "PE"
function class_abbr(c::PhosphatidylNmodifiedethanolamine)
    mod = head_abbr(first(getchaincomponent(c.backbone)))
    (occursin(" ", mod) || occursin("-", mod) || startswith(mod, r"\d")) ? string("PE-N(", mod, ")") : string("PE-N", mod)
end
class_abbr(::PhosphatidylNmethylethanolamine) = "PE-NMe"
class_abbr(::PhosphatidylNNdimethylethanolamine) = "PE-NMe2"
class_abbr(::Phosphatidylserine) = "PS"
function class_abbr(c::PhosphatidylNmodifiedserine)
    mod = head_abbr(first(getchaincomponent(c.backbone)))
    (occursin(" ", mod) || occursin("-", mod) || startswith(mod, r"\d")) ? string("PS-N(", mod, ")") : string("PS-N", mod)
end
function class_abbr(c::Phosphatidylinositol)
    pi = first(getchaincomponent(c.backbone))
    if isnothing(pi.substituent)
        "PI"
    else
        n = 0
        p = String[]
        for i in first(getchaincomponent(c.backbone)).substituent
            if i isa Phosphate
                n += 1
            elseif i isa Pair && last(i) isa Phosphate
                n += 1
                push!(p, string(Int(first(i)), "'"))
            end
        end
        if n < 1
            "PI"
        else
            string("PIP", n > 1 ? n : "", isempty(p) ? "" : "(" * join(p, ",") * ")")
        end
    end
end
class_abbr(::Phosphatidylglycerol) = "PG"
class_abbr(::Phosphatidylglycerolphosphate) = "PGP"
class_abbr(::Phosphatidylmethanol) = "PMeOH"
class_abbr(::Phosphatidylethanol) = "PEtOH"
class_abbr(c::Diradylglycerophosphate) = replace(head_abbr(c.backbone), r"P-\d*Glycerol[^-]*-*$" => "GP")

class_abbr(::Lysophosphatidicacid) = "LPA"
class_abbr(::Lysophosphatidylcholine) = "LPC"
class_abbr(::Lysophosphatidylethanolamine) = "LPE"
function class_abbr(c::LysophosphatidylNmodifiedethanolamine)
    mod = head_abbr(first(getchaincomponent(c.backbone)))
    (occursin(" ", mod) || occursin("-", mod) || startswith(mod, r"\d")) ? string("LPE-N(", mod, ")") : string("LPE-N", mod)
end
class_abbr(::LysophosphatidylNmethylethanolamine) = "LPE-NMe"
class_abbr(::LysophosphatidylNNdimethylethanolamine) = "LPE-NMe2"
class_abbr(::Lysophosphatidylserine) = "LPS"
function class_abbr(c::LysophosphatidylNmodifiedserine)
    mod = head_abbr(first(getchaincomponent(c.backbone)))
    (occursin(" ", mod) || occursin("-", mod) || startswith(mod, r"\d")) ? string("LPS-N(", mod, ")") : string("LPS-N", mod)
end
function class_abbr(c::Lysophosphatidylinositol)
    pi = first(getchaincomponent(c.backbone))
    if isnothing(pi.substituent)
        "LPI"
    else
        n = 0
        p = String[]
        for i in first(getchaincomponent(c.backbone)).substituent
            if i isa Phosphate
                n += 1
            elseif i isa Pair && last(i) isa Phosphate
                n += 1
                push!(p, string(Int(first(i)), "'"))
            end
        end
        if n < 1
            "LPI"
        else
            string("LPIP", n > 1 ? n : "", isempty(p) ? "" : "(" * join(p, ",") * ")")
        end
    end
end
class_abbr(::Lysophosphatidylglycerol) = "LPG"
class_abbr(::Lysophosphatidylglycerolphosphate) = "LPGP"
class_abbr(::Lysophosphatidylmethanol) = "LPMeOH"
class_abbr(::Lysophosphatidylethanol) = "LPEtOH"
class_abbr(c::Monoradylglycerophosphate) = replace(head_abbr(c.backbone), r"P-\d*Glycerol[^-]*-*$" => "LGP")

class_abbr(::Bisphosphatidicacid) = "BPA"
class_abbr(::Semilysobisphosphatidicacid) = "SLBPA"
class_abbr(::Lysobisphosphatidicacid) = "LBPA"
class_abbr(::Cardiolipin) = "CL"
class_abbr(::Monolysocardiolipin) = "MLCL"
class_abbr(::Dilysocardiolipin) = "DLCL"
class_abbr(::GlycerophosphoNacylethanolamine) = "GP-NAE"

class_abbr(::Ceramide) = "Cer"
class_abbr(c::CeramideBone) = string(head_abbr(c.headgroup), "-Cer")
class_abbr(c::Glycosylceramide{<: AbstractGlycan}) = head_abbr(c.headgroup)
class_abbr(c::Glycosylceramide{<: Lac}) = string(head_abbr(c.headgroup), "Cer")
class_abbr(c::Glycosylceramide{<: Glycan}) = string(head_abbr(c.headgroup), "Cer")
class_abbr(c::Glycosylceramide{<: GlyComp}) = string(head_abbr(c.headgroup), "Cer")
class_abbr(::CeramidePhosphate) = "CerP"
class_abbr(::Inositolphosphorylceramide) = "IPC"
class_abbr(c::Glycosylinositolphosphorylceramide) = string(head_abbr(first(getchaincomponent(c.headgroup))), "IPC")
class_abbr(::Ethanolaminephosphorylceramide) = "EPC"
class_abbr(::Mannosylinositolphosphorylceramide) = "MIPC"
class_abbr(::Mannosyldiinositolphosphorylceramide) = "M(IP)2C"
class_abbr(::Sphingomyelin) = "SM"
class_abbr(::Sulfonolipid) = "SL"

class_abbr(::SphingoidBase) = "SPB"
class_abbr(c::SphingoidBaseBone) = string(head_abbr(c.headgroup), "-SPB")
class_abbr(c::Glycosylsphingoidbase{<: AbstractGlycan}) = string("Lyso", head_abbr(c.headgroup))
class_abbr(c::Glycosylsphingoidbase{<: Lac}) = string(head_abbr(c.headgroup), "SPB")
class_abbr(c::Glycosylsphingoidbase{<: Glycan}) = string(head_abbr(c.headgroup), "-SPB")
class_abbr(c::Glycosylsphingoidbase{<: GlyComp}) = string(head_abbr(c.headgroup), "SPB")
class_abbr(::SphingoidBasePhosphate) = "SPBP"
class_abbr(::Lysoinositolphosphorylceramide) = "LIPC"
class_abbr(c::Lysoglycosylinositolphosphorylceramide) = string(head_abbr(first(getchaincomponent(c.backbone))), "LIPC")
class_abbr(::Lysoethanolaminephosphorylceramide) = "LEPC"
class_abbr(::Lysomannosylinositolphosphorylceramide) = "LMIPC"
class_abbr(::Lysomannosyldiinositolphosphorylceramide) = "LM(IP)2C"
class_abbr(::Lysosphingomyelin) = "LSM"
class_abbr(::Lysosulfonolipid) = "LSL"
function class_abbr(c::Acylceramide)
    mod = head_abbr(c.headgroup)
    isempty(mod) ? "ACer" : string(mod, "-ACer")
end
function class_abbr(c::Acylhexosylceramide)
    mod1, mod2 = head_abbr.(c.headgroup)
    isempty(mod1) ? string("A", mod2, "Cer") : string(mod1, "-A", mod2, "Cer")
end
function class_abbr(c::Acylsphingomyelin)
    mod1, mod2 = head_abbr.(c.headgroup)
    isempty(mod1) ? "ASM" : string(mod1, "-ASM")
end

head_abbr(x) = chemicalabbr(x)
head_abbr(::PhosphoricAcid) = "P"
head_abbr(glycan::GlyComp) = join(map(glycan.comp) do (x, n)
    c = head_abbr(x)
    if n > 1
        if endswith(c, r"\d")
            c = string("(", c, ")", Int(n))
        else
            c = string(c, Int(n))
        end
    elseif endswith(c, r"\d")
        c = string("(", c, ")")
    end
    c
end, "")

function head_abbr(glycan::Glycan)
    s = ""
    for (i, x) in enumerate(glycan.linkage)
        c = glycan.chain[i]
        xs = repr_linkage(x)
        cs = head_abbr(c)
        if c isa Glycan
            m = match(r"\(?((?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?$", cs)
            sp, si, sj = m
            mp, mi, mj = match(r"([αβab]*)(\d*)-(\d*)$", xs)
            if isempty(sp)
                sp = mp
            end
            if isempty(si)
                si = mi
            end
            if isempty(mj)
                mj = sj
            end
            x = string(sp, si, "-", mj)
            if startswith(x, "-") || startswith(x, "α") || startswith(x, "β")
                s = string(s, "[", cs[begin:m.match.offset], sp, si, "-", mj, "]")
            else
                s = string(s, "[", cs[begin:m.match.offset], "(", sp, si, "-", mj, ")]")
            end
        else
            if startswith(xs, "-") || startswith(xs, "α") || startswith(xs, "β")
                s = string(s, cs, xs)
            else
                s = string(s, cs, "(", xs, ")")
            end
        end
    end
    if length(glycan.linkage) == length(glycan.chain) - 1
        s = string(s, head_abbr(last(glycan.chain)))
    end
    s
end

function head_abbr(x::Hex)
    if x.substituent == [Sulfate() => 0x01]
       "SHex" 
    else
        chemicalname(x)
    end
end
function head_abbr(lipid::T) where {T <: FattyAcyl}
    ncarbon(lipid.chain) > 0 ? string(class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => "")) : 
        string("(", class_abbr(lipid), ")")
end
function head_abbr(lipid::T) where {T <: FattyAlcohol}
    ncarbon(lipid.chain) > 0 ? string(replace(repr_singlechain(lipid.chain), r"[OP]-" => "")) : "(Alk)"
end
function head_abbr(dc::DehydratedChemical)
    s = ""
    for (i, x) in enumerate(getchainlinkage(dc))
        c = getchaincomponent(dc)[i]
        xs = repr_linkage(x)
        cs = head_abbr(c)
        if ischainedchemical(c)
            m = match(r"\(((?:\(a)*(?:\(b)*[αβ]*)(\d*)-(\d*)\)?$", cs)
            sp, si, sj = m
            mp, mi, mj = match(r"([αβab]*)(\d*)-(\d*)$", xs)
            if isempty(sp)
                sp = mp
            end
            if isempty(si)
                si = mi
            end
            if isempty(mj)
                mj = sj
            end
            x = string(sp, si, "-", mj)
            if startswith(x, "-") || startswith(x, "α") || startswith(x, "β")
                s = string(s, "[", cs[begin:m.match.offset], x, "]")
            else
                s = string(s, "[", cs[begin:m.match.offset], "(", x, ")]")
            end
        else
            if startswith(xs, "-") || startswith(xs, "α") || startswith(xs, "β")
                s = string(s, cs, xs)
            else
                s = string(s, cs, "(", xs, ")")
            end
        end
    end
    if length(getchainlinkage(dc)) == length(getchaincomponent(dc)) - 1
        s = string(s, head_abbr(last(getchaincomponent(dc))))
    end
    s
end

"""
    sub_abbr(lipid)

Substituent abbreviation
"""
sub_abbr(x) = chemicalabbr(x)
sub_abbr(::Tauryl) = "T"
sub_abbr(::PhosphoricAcid) = "P"
function sub_abbr(x::Hex)
    if x == Hex([Sulfate() => 0x01])
       "SHex" 
    else
        chemicalname(x)
    end
end
function sub_abbr(dc::Union{Glycan, DehydratedChemical})
    s = ""
    for (i, x) in enumerate(getchainlinkage(dc))
        c = getchaincomponent(dc)[i]
        xs = repr_linkage(x)
        cs = sub_abbr(c)
        mp, mi, mj = match(r"([αβab]*)(\d*)-(\d*)$", xs)
        if ischainedchemical(c)
            m = match(r"^\(?(\d*)-(\d*)([abαβ]?)\)?", cs)
            sj, si, sp = m
            if isempty(sp)
                sp = mp
            end
            if isempty(si)
                si = mi
            end
            if isempty(mj)
                mj = sj
            end
            x = string(mj, "-", si, sp)
            if endswith(x, "-") || endswith(x, "α") || endswith(x, "β")
                s = string("[", x, cs[m.match.offset + m.match.ncodeunits + 1:end], "]", s)
            else
                s = string("[","(", x, ")", cs[m.match.offset + m.match.ncodeunits + 1:end], "]", s)
            end
        else
            x = string(mj, "-", mi, mp)
            if endswith(x, "-") || endswith(x, "α") || endswith(x, "β")
                s = string(x, cs, s)
            else
                s = string("(", x, ")", cs, s)
            end
        end
    end
    if length(getchainlinkage(dc)) == length(getchaincomponent(dc)) - 1
        s = string(sub_abbr(last(getchaincomponent(dc))), s)
    end
    s
end

function sub_abbr(lf::XLinkedFunctionalGroup)
    l = sub_abbr(lf.xlinkage)
    f = sub_abbr(lf.functionalgroup) # depends on linkage, no internal linkage
    if occursin(" ", f) || startswith(f, r"\d") || occursin("-", f)
        string(l, "(", f, ")")
    else
        string(l, f)
    end
end

function sub_abbr(lipid::FattyAlcohol)
    replace(chemicalname(lipid), r"^FOH[^\s]* " => "")
end

function sub_abbr(dc::Substituent)
    s = sub_abbr(dc.molecule)
    # if leavinggroup(dc) == Dehydroxy()
    #     dp = dehydroxyposition(last(dc.molecule))
    # elseif leavinggroup(dc) == Dehydrogen()
    #     dp = dehydroxyposition(last(dc.molecule))
    # else
    #     return s
    # end
    s
end
"""
    repr_singlechain(c::CarbonChain)

Representation of single carbon chain
"""
repr_singlechain(c::CarbonChain{<: SPB}) = string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Acyl}) = string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Alkyl}) = string("O-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Alkenyl}) = string("P-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
function repr_singlechain(c::CarbonChain{<: T}) where {T <: Tuple}
    if SPB in T.parameters
        return string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
    elseif Alkenyl in T.parameters # assume only P and A
        n = count(==(Alkenyl), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkenyl chain"))
        pre = string(pre, "P-")
    elseif Alkyl in T.parameters # assume only O and A
        n = count(==(Alkyl), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkyl chain"))
        pre = string(pre, "O-")
    else
        pre = ""
    end
    string(pre, ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
end

# to full structure level
# delete headgroup
# delete charged chain components

"""
    repr_smiles_carbonchain(lipid)

SMILES representation of lipid without headgroup and charged functional group
"""
repr_smiles_carbonchain(back::Glycerol) = repr_smiles(back)
repr_smiles_carbonchain(lipid::Diradylglycerol) = repr_smiles(lipid)
repr_smiles_carbonchain(lipid::Monoradylglycerol) = repr_smiles(lipid)
function repr_smiles(lipid::Monoradylglycerol)
    # check
    p = only(decode_sn(lipid))
    c = repr_smiles_carbonchain(lipid.chain)
    smi = repr_smiles(lipid.backbone)
    r = collect(eachmatch(r"\(O\)", smi))
    next = r[lastindex(r) - p + 1].match.offset
    string(smi[begin:next + 2], c, smi[next + 3:end])
end
function repr_smiles(lipid::Diradylglycerol)
    # check
    position = decode_sn(lipid)
    p = sortperm(position; rev = true)
    position = position[p]
    cs = [repr_smiles_carbonchain(c) for c in lipid.chain[p]]
    smi = repr_smiles(lipid.backbone)
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[[length(r) - t + 1 for t in position]], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

repr_smiles_carbonchain(lipid::Glycerophospholipid) = repr_smiles_carbonchain(dissociate_headgroup(lipid))

function repr_smiles_carbonchain(lipid::Ceramide)
    # check
    cs = [repr_smiles_carbonchain(c) for c in lipid.chain]
    id = findfirst("(N)", first(cs))
    string(first(cs)[begin:id[2]], cs[begin + 1], first(cs)[last(id):end])
end

repr_smiles_carbonchain(lipid::SphingoidBase) = repr_smiles_carbonchain(lipid.chain)
repr_smiles_carbonchain(lipid::CeramideBone) = repr_smiles_carbonchain(dissociate_headgroup(lipid))
repr_smiles_carbonchain(lipid::SphingoidBaseBone) = repr_smiles_carbonchain(dissociate_headgroup(lipid))

initial_smiles(::CarbonChain{Acyl}) = (["C(=O)"], [0])
initial_smiles(::CarbonChain{Alkyl}) = (["C"], [0])
initial_smiles(::CarbonChain{Alkenyl}) = (["\\C", "=C", "/C"], [0, 0, 0])
initial_smiles(::CarbonChain{SPB}) = (["C", "C(N)"], [0, 0])

function push_smiles_db!(quec, posc, slash)
    if length(quec) > 1
        quec[begin] = string("=", quec[begin])
        quec[begin + 1] = string(slash, quec[begin + 1])
    else
        if isempty(quec)
            push!(quec, "=X")
            push!(posc, 0)
        else
            quec[begin] = string("=", quec[begin])
        end
        push!(quec, string(slash, "X"))
        push!(posc, 0)
    end
end

function repr_smiles_carbonchain(chain::CarbonChain{<: CarbonChainType})
    i = j = k = 1
    pos = iszero(chain.doublebond) ? [0x00] : sort!([divrem(db, 3) for db in chain.doublebond]; by = first)
    subs = isnothing(chain.substituent) ? Pair[] : sort(chain.substituent; by = first)
    s = ""
    cn = 0
    quec, posc = initial_smiles(chain)
    while i <= ncarbon(chain)
        if isempty(quec)
            prec = "C"
            numc = 0
        else
            prec = replace(popfirst!(quec), r"[XC]" => "C")
            numc = popfirst!(posc)
        end
        if numc > 0
            prec *= string(numc) 
        end
        if !startswith(prec, r"[\[\\/=#]*C") # -O-, -OO-
            s *= prec
            continue
        end
        if j <= lastindex(pos) && i == first(pos[j])
            if last(pos[j]) == 0
                if isempty(quec)
                    push!(quec, "=X")
                    push!(posc, 0)
                else
                    quec[begin] = string("=", quec[begin])
                end
            elseif last(pos[j]) == 1
                if startswith(prec, '\\')
                    push_smiles_db!(quec, posc, "/")
                elseif startswith(prec, '/')
                    push_smiles_db!(quec, posc, "\\")
                else
                    prec = replace(prec, r"[XC]" => "\\C") # replace ? 
                    push_smiles_db!(quec, posc, "/")
                end
            elseif last(pos[j]) == 2
                if startswith(prec, '\\')
                    push_smiles_db!(quec, posc, "\\")
                elseif startswith(prec, '/')
                    push_smiles_db!(quec, posc, "/")
                else
                    prec = replace(prec, r"[XC]" => "\\C")
                    push_smiles_db!(quec, posc, "\\")
                end
            end
            j += 1
        end
        while k <= lastindex(subs) && first(subs[k]) == i
            smi = repr_smiles(last(subs[k]))
            if !endswith(smi, ")")
                bsmi, msmi = match(r"^((?:\([^\)]*\))*)(.*)", smi)
                if isempty(bsmi) || !occursin(r"\d+", bsmi)
                    add = [e for (e, ) in eachmatch(r"(=*#*\[*[A-Z][a-z]*@*H*\]*)\d*", msmi)] # NO OTHER BRANCH
                    for i in eachindex(add)
                        if i > lastindex(posc)
                            push!(quec, add[i])
                            push!(posc, 0)
                        else
                            quec[i] = replace(quec[i], "X" => add[i])
                        end
                    end
                else
                    add = [string(E) => isempty(n) ? 0 : parse(Int, string(n)) for (E, n) in eachmatch(r"(=*#*\[*[A-Z][a-z]*@*H*\]*)(\d*)", msmi)] # NO OTHER BRANCH
                    new = last.(add)
                    org = new .- cn
                    for i in eachindex(new)
                        if i > lastindex(posc)
                            push!(quec, first(add[i]))
                            push!(posc, new[i])
                        elseif posc[i] > 0
                            if new[i] - posc[i] > 0
                                new[i:end] .-= 1
                                new[i] = posc[i]
                            end
                            quec[i] = replace(quec[i], "X" => first(add[i]))
                        else
                            posc[i] = new[i]
                            quec[i] = replace(quec[i], "X" => first(add[i]))
                        end
                    end
                    id = findall(>(0), org)
                    bsmi = replace(bsmi, (string.(org[id]) .=> string.(new[id]))...)
                    cn = maximum(posc)
                end
                prec = string(prec, bsmi)
            else
                prec = string(prec, smi)
            end
            k += 1
        end
        i += 1
        s *= prec
    end
    s
end

function Base.show(io::IO, level::AnnotationLevel)
    print(io, @match level begin
        &specieslevel                   => "specieslevel"
        &molecularspecieslevel          => "molecularspecieslevel"
        &phosphatepositionlevel         => "phosphatepositionlevel"
        &snpositionlevel                => "snpositionlevel"
        &dbpositionpartiallevel         => "dbpositionpartiallevel"
        &dbpositionlevel                => "dbpositionlevel"
        &dbconfigpartiallevel           => "dbconfigpartiallevel"
        &dbconfiglevel                  => "dbconfiglevel"
        &structuredefinedpartiallevel   => "structuredefinedpartiallevel"
        &structuredefinedlevel          => "structuredefinedlevel"
        &structurepositionpartiallevel  => "structurepositionpartiallevel"
        &structurepositionlevel         => "structurepositionlevel"
        &structureconfigpartiallevel    => "structureconfigpartiallevel"
        &structureconfiglevel           => "structureconfiglevel"
        &fullstructurelevel             => "fullstructurelevel"
        &completestructurelevel         => "completestructurelevel"
    end)
end


function Base.show(io::IO, level::PassLevel)
    print(io, @match level begin
        &passphosphatepositionlevel     => "passphosphatepositionlevel"
        &passsnpositionlevel            => "passsnpositionlevel"
    end)
end