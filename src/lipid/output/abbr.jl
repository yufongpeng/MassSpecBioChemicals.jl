# Lipid
# MassSpecChemicals
getchemicalattr(lipid::T, ::Val{:abbreviation}; kwargs...) where {T <: Lipid} = getchemicalattr(lipid, Val(:name); kwargs...)

chiral_class_abbr(x) = class_abbr(x)
chiral_class_abbr(x::Union{Glycerolipid, Glycerophospholipid}) = string(repr_class_rs(x), class_abbr(x))
function chiral_class_abbr(c::Omodifiedmonoradylglycerol)
    bone = getchaincomponent(c.backbone)
    if length(bone) == 2 
        b = parentchemical(first(bone))
        b isa Quinovose{DForm, PyranoseForm} && b.substituent == [0x06 => Sulfo()] && getchainlinkage(c.backbone) == [α(0x01) => lk(0x03)] && return string(repr_class_rs(c), "SQMG")
        b isa Galactose{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(c.backbone) == [β(0x01) => lk(0x03)] && return string(repr_class_rs(c), "MGMG")
        b isa GlucuronicAcid{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(c.backbone) == [α(0x01) => lk(0x03)] && return string(repr_class_rs(c), "GlcAMG")
    elseif length(bone) == 3 
        b = parentchemical(first(bone))
        b2 = parentchemical(bone[begin + 1])
        b isa Galactose{DForm, PyranoseForm} && b2 isa Galactose{DForm, PyranoseForm} && getchainlinkage(c.backbone) == [α(0x01) => lk(0x06), β(0x01) => lk(0x03)] && return string(repr_class_rs(c), "DGMG")
    end
    b = head_abbr(c.backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "MG", r)
end
function chiral_class_abbr(c::Omodifieddiradylglycerol) 
    bone = getchaincomponent(c.backbone)
    if length(bone) == 2 
        b = parentchemical(first(bone))
        b isa Quinovose{DForm, PyranoseForm} && b.substituent == [0x06 => Sulfo()] && getchainlinkage(c.backbone) == [α(0x01) => lk(0x03)] && return string(repr_class_rs(c), "SQDG")
        b isa Galactose{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(c.backbone) == [β(0x01) => lk(0x03)] && return string(repr_class_rs(c), "MGDG")
        b isa GlucuronicAcid{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(c.backbone) == [α(0x01) => lk(0x03)] && return string(repr_class_rs(c), "GlcADG")
    elseif length(bone) == 3 
        b = parentchemical(first(bone))
        b2 = parentchemical(bone[begin + 1])
        b isa Galactose{DForm, PyranoseForm} && b2 isa Galactose{DForm, PyranoseForm} && getchainlinkage(c.backbone) == [α(0x01) => lk(0x06), β(0x01) => lk(0x03)] && return string(repr_class_rs(c), "DGDG")
    end
    b = head_abbr(c.backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "DG", r)
end

chiral_class_abbr(x::FattyAcylCarnitine) = string(repr_class_rs(x), class_abbr(x))
chiral_class_abbr(x::Union{Phosphatidicacid, Phosphatidylcholine, Phosphatidylethanol, Phosphatidylethanolamine, Phosphatidylinositol, Phosphatidylmethanol, PhosphatidylNmethylethanolamine, PhosphatidylNmodifiedethanolamine, PhosphatidylNNdimethylethanolamine, Lysophosphatidicacid, Lysophosphatidylcholine, Lysophosphatidylethanol, Lysophosphatidylethanolamine, Lysophosphatidylinositol, Lysophosphatidylmethanol, LysophosphatidylNmethylethanolamine, LysophosphatidylNmodifiedethanolamine, LysophosphatidylNNdimethylethanolamine}) = 
    string(repr_class_rs(x), class_abbr(x))
chiral_class_abbr(c::Monoradylglycerophosphate) = replace(head_abbr(c.backbone), r"P-\d*\(*[RS]*\)*-*Glycerol[^-]*-*$" => string(repr_class_rs(c), "LGP"))
chiral_class_abbr(c::Diradylglycerophosphate) = replace(head_abbr(c.backbone), r"P-\d*\(*[RS]*\)*-*Glycerol[^-]*-*$" => string(repr_class_rs(c), "GP"))

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
    T = PROTEIN_3LETTER_AA[string(x)]
    @eval class_abbr(c::NacylAmine{<: includeSIL($T)}) = string("NA", letter3_abbr(c.backbone))
    @eval chiral_class_abbr(c::NacylAmine{<: includeSIL($T)}) = string(repr_class_rs(c), class_abbr(c))
end
class_abbr(::FattyAcylCarnitine) = "CAR"
class_abbr(::FattyAcylCoA) = "CoA"
class_abbr(::FattyAcylEstolid) = "FAHFA"

class_abbr(::Monoradylglycerol) = "MG"
class_abbr(::Diradylglycerol) = "DG"
class_abbr(::Triradylglycerol) = "TG"

# class_abbr(::Sulfoquinovosylmonoradylglycerol) = "SQMG"
# class_abbr(::Sulfoquinovosyldiradylglycerol) = "SQDG"
# class_abbr(::Monogalactosylmonoradylglycerol) = "MGMG"
# class_abbr(::Monogalactosyldiradylglycerol) = "MGDG"
# class_abbr(::Digalactosylmonoradylglycerol) = "DGMG"
# class_abbr(::Digalactosyldiradylglycerol) = "DGDG"
# class_abbr(::Glucuronosylmonoradylglycerol) = "GlcAMG"
# class_abbr(::Glucuronosyldiradylglycerol) = "GlcADG"

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
            if i isa Pair && first(i) isa Phosphoryl
                n += last(i)
            elseif i isa Pair && last(i) isa Phosphoryl
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
            if i isa Pair && first(i) isa Phosphoryl
                n += last(i)
            elseif i isa Pair && last(i) isa Phosphoryl
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
    string(head_abbr(c.headgroup), "-ACer")
end
# function class_abbr(c::Acylhexosylceramide)
#     mod1, mod2 = head_abbr.(c.headgroup)
#     isempty(mod1) ? string("A", mod2, "Cer") : string(mod1, "-A", mod2, "Cer")
# end
function class_abbr(c::Acylsphingomyelin)
    mod1, mod2 = head_abbr.(c.headgroup)
    isempty(mod1) ? "ASM" : string(mod1, "-ASM")
end

"""
    head_abbr(x)

Head group abbreviation
"""
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

function head_abbr(x::Hexose)
    if x.substituent == [Sulfo() => 0x01]
       "SHex" 
    else
        chemicalabbr(x)
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

function head_abbr(dc::Substituent)
    head_abbr(dc.chemical)
end

"""
    sub_abbr(lipid)

Substituent abbreviation
"""
sub_abbr(x) = chemicalabbr(x)
sub_abbr(::Tauryl) = "T"
sub_abbr(::PhosphoricAcid) = "P"
function sub_abbr(x::Hexose)
    if x == Hexose([Sulfo() => 0x01])
       "SHex" 
    else
        chemicalabbr(x)
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
    f = lf.xlinkage == CarboxylicLinkage() ? head_abbr(lf.functionalgroup) : sub_abbr(lf.functionalgroup) # depends on linkage, no internal linkage
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
    s = sub_abbr(dc.chemical)
    # if leavinggroup(dc) == Dehydroxy()
    #     dp = dehydroxyposition(last(dc.chemical))
    # elseif leavinggroup(dc) == Dehydrogen()
    #     dp = dehydroxyposition(last(dc.chemical))
    # else
    #     return s
    # end
    s
end