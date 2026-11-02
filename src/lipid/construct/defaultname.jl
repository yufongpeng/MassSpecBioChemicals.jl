"""
    lipiddefaultname(::Type{<:Lipid}, components...)
    lipiddefaultname(::Type{<:FattyAcyl}, backbone, chain, position = 0x00)
    lipiddefaultname(::Type{<:Glycerolipid}, backbone, chain, sn, chirality)
    lipiddefaultname(::Type{<:GlycerophosphoNacylethanolamine}, backbone, chain, chirality)
    lipiddefaultname(::Type{<:Glycerophospholipid}, backbone, chain, sn, chirality)
    lipiddefaultname(::Type{<:Sphingolipid}, backbone, chain, position, chirality)

Generate the default name for a lipid based on its type and structural components.
"""
lipiddefaultname(::Type{T}, chain) where {T <: Lipid} =
    string(nameof(Base.typename(T).wrapper), " ", repr_singlechain(chain))
lipiddefaultname(::Type{T}, backbone, chain, args...) where {T <: Lipid} =
    string(default_lipid_classabbr(T, backbone), " ", repr_singlechain(chain))
lipiddefaultname(::Type{<:LipidChemical}, chemical) = chemicalname(chemical)
lipiddefaultname(::Type{T}, backbone, chain, position = 0x00) where {T <: FattyAcyl} =
    string(default_lipid_classabbr(T, backbone), " ", replace(repr_singlechain(chain), r"[OP]-" => ""))
function lipiddefaultname(::Type{T}, backbone, chain, position = 0x00) where {T <: Union{NacylAlkylAmine, FattyAcylEster}}
    if ncarbon(backbone.chain) == 0
        string(default_lipid_classabbr(T, backbone), " ", replace(repr_singlechain(chain), r"[OP]-" => ""))
    else
        string(replace(chemicalname(backbone), default_lipid_classabbr(typeof(backbone), nothing) => default_lipid_classabbr(T, backbone)), "/", replace(repr_singlechain(chain), r"[OP]-" => ""))
    end
end

function lipiddefaultname(::Type{T}, backbone, chain, position) where {T <: FattyAcylEstolid}
    if ncarbon(backbone.chain) == 0
        string(default_lipid_classabbr(T, backbone), " ", replace(repr_singlechain(chain), r"[OP]-" => ""))
    elseif iszero(position)
        string(default_lipid_classabbr(T, backbone), " ", replace(repr_singlechain(chain), r"[OP]-" => ""), "/", replace(chemicalname(backbone), string(default_lipid_classabbr(typeof(backbone), nothing), " ") => ""))
    else
        i = Int(position)
        id = findfirst(x -> (Int(first(x)) == i) && (last(x) == Hydroxy()), backbone.chain.substituent)
        sub = backbone.chain.substituent[setdiff(eachindex(backbone.chain.substituent), [id])]
        id = findfirst(x -> Int(first(x)) == i, backbone.chain.chirality)
        rs = last(backbone.chain.chirality[id])
        l = MonoFattyAcyl(backbone.backbone, CarbonChain{Acyl}(backbone.chain.carbon, backbone.chain.doublebond, sub, backbone.chain.chirality, backbone.chain.isotopiclabel))
        string(default_lipid_classabbr(T, backbone), " ", replace(repr_singlechain(chain), r"[OP]-" => ""), "/", i, "O(", l, ")", repr_chirality(rs))
    end
end

lipiddefaultname(::Type{<:GlycerophosphoNacylethanolamine}, backbone, chain, chirality) =
    string(default_lipid_chirality_prefix(chirality), "GP-NAE ", repr_singlechain(chain))

function lipiddefaultname(::Type{T}, backbone, chain, sn, chirality) where {T <: Union{Glycerolipid, Glycerophospholipid}}
    chains = tuplize(chain)
    class = default_lipid_classabbr(T, backbone, chirality)
    positions = decode_sn(T, length(chains), sn)
    if any(iszero, positions) && any(>(0), positions)
        pos = chainposition(T)
        return string(class, " ", join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(positions, chains)], "_"))
    elseif any(iszero, positions)
        return string(class, " ", join(repr_singlechain.(chains), "_"))
    end
    names = fill("0:0", nchainposition(T))
    for (c, p) in zip(chains, positions)
        names[p] = repr_singlechain(c)
    end
    string(class, " ", join(names, "/"))
end

function lipiddefaultname(::Type{T}, backbone, chain, position, chirality) where {T <: Sphingolipid}
    string(default_lipid_classabbr(T, backbone), decode_headgroup_position(T, position, chirality),
        " ", join(repr_singlechain.(tuplize(chain)), "/"))
end

lipiddefaultname(::Type{<:SubstitutedSterol}, backbone) = chemicalname(backbone)

"""
    default_lipid_classabbr(::Type{<:Lipid}, backbone)
    default_lipid_classabbr(::Type{<:Lipid}, headgroup)

Default lipid class abbreviation.
"""
default_lipid_classabbr(::Type{T}, backbone) where {T <: Lipid} =
    string(nameof(Base.typename(T).wrapper))

default_lipid_classabbr(::Type{<:Hydrocarbon}, backbone) = "HC"
default_lipid_classabbr(::Type{<:FattyAcid}, backbone) = "FA"
default_lipid_classabbr(::Type{<:FattyAldehyde}, backbone) = "FAL"
default_lipid_classabbr(::Type{<:FattyAlcohol}, backbone) = "FOH"
default_lipid_classabbr(::Type{<:FattyAmide}, backbone) = "FAM"
default_lipid_classabbr(::Type{<:FattyAmine}, backbone) = "FN"
default_lipid_classabbr(::Type{<:FattyAcylCoA}, backbone) = "CoA"
default_lipid_classabbr(::Type{<:FattyAcylEstolid}, backbone) = "FAHFA"
default_lipid_classabbr(::Type{<:WaxEster}, backbone) = "WE"
# default_lipid_classabbr(::Type{<:FattyAcylEster}, backbone) = "Lipid"
default_lipid_classabbr(::Type{<:NacylAlkylAmine}, backbone) = "NA"
default_lipid_classabbr(::Type{<:NacylAmine{<:Ethanolamine}}, backbone) = "NAE"
default_lipid_classabbr(::Type{<:NacylAmine{<:Taurine}}, backbone) = "NAT"
default_lipid_classabbr(::Type{<:NacylAmine}, backbone) =
    applicable(letter3_abbr, backbone) ? string("NA", letter3_abbr(backbone)) : "NA"
for x in NAAA
    T = PROTEIN_3LETTER_AA[string(x)]
    @eval default_lipid_classabbr(::Type{<:NacylAmine{<: $T}}, backbone::S) where S =
        string(default_backbone_chirality_prefix(S), "NA", letter3_abbr(backbone))
end
default_lipid_classabbr(::Type{<:FattyAcylCarnitine{B}}, backbone) where B =
    string(default_backbone_chirality_prefix(B), "CAR")

function default_lipid_classabbr(::Type{T}, backbone, chirality) where {T <: Glycerolipid}
    string(default_lipid_chirality_prefix(chirality), default_achiral_lipid_classabbr(T, backbone))
end
default_lipid_classabbr(::Type{<:Radyldiglycerol}, backbone, chirality) = head_abbr(backbone)
default_lipid_classabbr(::Type{<:Radyltriglycerol}, backbone, chirality) = head_abbr(backbone)
function default_lipid_classabbr(::Type{T}, backbone, chirality) where {T <: Glycerophospholipid}
    string(default_lipid_chirality_prefix(chirality), default_achiral_lipid_classabbr(T, backbone))
end

"""
    default_lipid_chirality_prefix(chirality)

Return the default lipid chirality prefix for the given chirality.
"""
function default_lipid_chirality_prefix(chirality)
    rs = map(x -> repr_chirality(x; bracket = false), tuplize(chirality))
    labels = [string(x, repeat("'", i - 1)) for (i, x) in enumerate(rs) if !isempty(x)]
    isempty(labels) ? "" : string("(", join(labels, ","), ")-")
end

"""
    default_backbone_chirality_prefix(::Type)

Return the default backbone chirality prefix for the given backbone type.
"""
function default_backbone_chirality_prefix(::Type{B}) where B
    B isa UnionAll && return ""
    rs = repr_chirality(first(B.parameters)())
    isempty(rs) ? "" : string(rs, "-")
end

"""
    default_achiral_lipid_classabbr(::Type{<:Lipid}, backbone)

Return the default achiral lipid class abbreviation for the given lipid type and backbone.
"""
default_achiral_lipid_classabbr(::Type{<:Monoradylglycerol}, backbone) = "MG"
default_achiral_lipid_classabbr(::Type{<:Diradylglycerol}, backbone) = "DG"
default_achiral_lipid_classabbr(::Type{<:Triradylglycerol}, backbone) = "TG"
default_achiral_lipid_classabbr(::Type{<:Radyldiglycerol}, backbone) = head_abbr(backbone)
default_achiral_lipid_classabbr(::Type{<:Radyltriglycerol}, backbone) = head_abbr(backbone)
# default_achiral_lipid_classabbr(::Type{<:Omodifiedmonoradylglycerol}, backbone) =
    # replace(head_abbr(backbone), r"(.*-)\d*Glycerol((?:\[[^-]\])?)" => s"\1MG\2")
# default_achiral_lipid_classabbr(::Type{<:Omodifieddiradylglycerol}, backbone) =
    # replace(head_abbr(backbone), r"(.*-)\d*Glycerol((?:\[[^-]\])?)" => s"\1DG\2")
# default_achiral_lipid_classabbr(::Type{<:Omodifiedradylglycerol}, backbone) =
#     replace(head_abbr(backbone), r"(.*-)\d*Glycerol((?:\[[^-]\])?)" => s"\1G\2")
default_achiral_lipid_classabbr(::Type{T}, backbone) where {T <: Glycerolipid} =
    string(nameof(Base.typename(T).wrapper))

function default_achiral_lipid_classabbr(::Type{<:Omodifiedmonoradylglycerol}, backbone)
    bone = getchaincomponent(backbone)
    if length(bone) == 2 
        b = parentchemical(first(bone))
        b isa Quinovose{DForm, PyranoseForm} && b.substituent == [0x06 => Sulfo()] && getchainlinkage(backbone) == [α(0x01) => lk(0x03)] && return "SQMG"
        b isa Galactose{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(backbone) == [β(0x01) => lk(0x03)] && return "MGMG"
        b isa GlucuronicAcid{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(backbone) == [α(0x01) => lk(0x03)] && return "GlcAMG"
    elseif length(bone) == 3 
        b = parentchemical(first(bone))
        b2 = parentchemical(bone[begin + 1])
        b isa Galactose{DForm, PyranoseForm} && b2 isa Galactose{DForm, PyranoseForm} && getchainlinkage(backbone) == [α(0x01) => lk(0x06), β(0x01) => lk(0x03)] && return "DGMG"
    end
    b = head_abbr(backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "MG", r)
end

function default_achiral_lipid_classabbr(::Type{<:Omodifieddiradylglycerol}, backbone)
    bone = getchaincomponent(backbone)
    if length(bone) == 2 
        b = parentchemical(first(bone))
        b isa Quinovose{DForm, PyranoseForm} && b.substituent == [0x06 => Sulfo()] && getchainlinkage(backbone) == [α(0x01) => lk(0x03)] && return "SQDG"
        b isa Galactose{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(backbone) == [β(0x01) => lk(0x03)] && return "MGDG"
        b isa GlucuronicAcid{DForm, PyranoseForm} && isnothing(b.substituent) && getchainlinkage(backbone) == [α(0x01) => lk(0x03)] && return "GlcADG"
    elseif length(bone) == 3 
        b = parentchemical(first(bone))
        b2 = parentchemical(bone[begin + 1])
        b isa Galactose{DForm, PyranoseForm} && b2 isa Galactose{DForm, PyranoseForm} && getchainlinkage(backbone) == [α(0x01) => lk(0x06), β(0x01) => lk(0x03)] && return "DGDG"
    end
    b = head_abbr(backbone)
    l, r = match(r"(.*-)\d*Glycerol((?:\[[^-]\])?)", b)
    string(l, "DG", r)
end

default_achiral_lipid_classabbr(::Type{<:Lysophosphatidicacid}, backbone) = "LPA"
default_achiral_lipid_classabbr(::Type{<:Phosphatidicacid}, backbone) = "PA"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylcholine}, backbone) = "LPC"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylcholine}, backbone) = "PC"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylethanolamine}, backbone) = "LPE"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylethanolamine}, backbone) = "PE"
default_achiral_lipid_classabbr(::Type{<:LysophosphatidylNmethylethanolamine}, backbone) = "LPE-NMe"
default_achiral_lipid_classabbr(::Type{<:PhosphatidylNmethylethanolamine}, backbone) = "PE-NMe"
default_achiral_lipid_classabbr(::Type{<:LysophosphatidylNNdimethylethanolamine}, backbone) = "LPE-NMe2"
default_achiral_lipid_classabbr(::Type{<:PhosphatidylNNdimethylethanolamine}, backbone) = "PE-NMe2"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylserine}, backbone) = "LPS"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylserine}, backbone) = "PS"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylglycerol}, backbone) = "LPG"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylglycerol}, backbone) = "PG"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylglycerolphosphate}, backbone) = "LPGP"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylglycerolphosphate}, backbone) = "PGP"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylmethanol}, backbone) = "LPMeOH"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylmethanol}, backbone) = "PMeOH"
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylethanol}, backbone) = "LPEtOH"
default_achiral_lipid_classabbr(::Type{<:Phosphatidylethanol}, backbone) = "PEtOH"
default_achiral_lipid_classabbr(::Type{<:Bisphosphatidicacid}, backbone) = "BPA"
default_achiral_lipid_classabbr(::Type{<:Semilysobisphosphatidicacid}, backbone) = "SLBPA"
default_achiral_lipid_classabbr(::Type{<:Lysobisphosphatidicacid}, backbone) = "LBPA"
default_achiral_lipid_classabbr(::Type{<:Cardiolipin}, backbone) = "CL"
default_achiral_lipid_classabbr(::Type{<:Monolysocardiolipin}, backbone) = "MLCL"
default_achiral_lipid_classabbr(::Type{<:Dilysocardiolipin}, backbone) = "DLCL"
default_achiral_lipid_classabbr(::Type{<:LysophosphatidylNmodifiedethanolamine}, backbone) =
    _modified_headgroup_classabbr("LPE-N", backbone)
default_achiral_lipid_classabbr(::Type{<:PhosphatidylNmodifiedethanolamine}, backbone) =
    _modified_headgroup_classabbr("PE-N", backbone)
default_achiral_lipid_classabbr(::Type{<:LysophosphatidylNmodifiedserine}, backbone) =
    _modified_headgroup_classabbr("LPS-N", backbone)
default_achiral_lipid_classabbr(::Type{<:PhosphatidylNmodifiedserine}, backbone) =
    _modified_headgroup_classabbr("PS-N", backbone)
default_achiral_lipid_classabbr(::Type{<:Lysophosphatidylinositol}, backbone) =
    _inositolphospholipid_classabbr("LPI", "LPIP", backbone)
default_achiral_lipid_classabbr(::Type{<:Phosphatidylinositol}, backbone) =
    _inositolphospholipid_classabbr("PI", "PIP", backbone)
default_achiral_lipid_classabbr(::Type{<:Monoradylglycerophosphate}, backbone) =
    replace(head_abbr(backbone), r"P-\d*\(*[RS]*\)*-*Glycerol[^-]*-*$" => "LGP")
default_achiral_lipid_classabbr(::Type{<:Diradylglycerophosphate}, backbone) =
    replace(head_abbr(backbone), r"P-\d*\(*[RS]*\)*-*Glycerol[^-]*-*$" => "GP")
default_achiral_lipid_classabbr(::Type{<:GlycerophosphoNacylethanolamine}, backbone) = "GP-NAE"
default_achiral_lipid_classabbr(::Type{T}, backbone) where {T <: Glycerophospholipid} =
    string(nameof(Base.typename(T).wrapper))

function _modified_headgroup_classabbr(prefix, backbone)
    mod = head_abbr(first(getchaincomponent(backbone)))
    (occursin(" ", mod) || occursin("-", mod) || startswith(mod, r"\d")) ?
        string(prefix, "(", mod, ")") : string(prefix, mod)
end

function _inositolphospholipid_classabbr(base, phosphorylated, backbone)
    inositol = first(getchaincomponent(backbone))
    isnothing(inositol.substituent) && return base
    n = 0
    positions = String[]
    for sub in inositol.substituent
        if sub isa Pair && first(sub) isa Phosphoryl
            n += last(sub)
        elseif sub isa Pair && last(sub) isa Phosphoryl
            n += 1
            push!(positions, string(Int(first(sub)), "'"))
        end
    end
    n < 1 && return base
    string(phosphorylated, n > 1 ? n : "",
        isempty(positions) ? "" : "(" * join(positions, ",") * ")")
end

default_lipid_classabbr(::Type{<:Ceramide}, headgroup) = "Cer"
default_lipid_classabbr(::Type{<:CeramideBone}, headgroup) = string(head_abbr(headgroup), "-Cer")
default_lipid_classabbr(::Type{<:Glycosylceramide{<: AbstractGlycan}}, headgroup) = head_abbr(headgroup)
default_lipid_classabbr(::Type{<:Glycosylceramide{<: Lac}}, headgroup) = string(head_abbr(headgroup), "Cer")
default_lipid_classabbr(::Type{<:Glycosylceramide{<: Glycan}}, headgroup) = string(head_abbr(headgroup), "Cer")
default_lipid_classabbr(::Type{<:Glycosylceramide{<: GlyComp}}, headgroup) = string(head_abbr(headgroup), "Cer")
default_lipid_classabbr(::Type{<:CeramidePhosphate}, headgroup) = "CerP"
default_lipid_classabbr(::Type{<:Inositolphosphorylceramide}, headgroup) = "IPC"
default_lipid_classabbr(::Type{<:Glycosylinositolphosphorylceramide}, headgroup) = string(head_abbr(first(getchaincomponent(headgroup))), "IPC")
default_lipid_classabbr(::Type{<:Ethanolaminephosphorylceramide}, headgroup) = "EPC"
default_lipid_classabbr(::Type{<:Mannosylinositolphosphorylceramide}, headgroup) = "MIPC"
default_lipid_classabbr(::Type{<:Mannosyldiinositolphosphorylceramide}, headgroup) = "M(IP)2C"
default_lipid_classabbr(::Type{<:Sphingomyelin}, headgroup) = "SM"
default_lipid_classabbr(::Type{<:Sulfonolipid}, headgroup) = "SL"

default_lipid_classabbr(::Type{<:SphingoidBase}, headgroup) = "SPB"
default_lipid_classabbr(::Type{<:SphingoidBaseBone}, headgroup) = string(head_abbr(headgroup), "-SPB")
default_lipid_classabbr(::Type{<:Glycosylsphingoidbase{<: AbstractGlycan}}, headgroup) = string("Lyso-", head_abbr(headgroup))
default_lipid_classabbr(::Type{<:Glycosylsphingoidbase{<: Lac}}, headgroup) = string(head_abbr(headgroup), "SPB")
default_lipid_classabbr(::Type{<:Glycosylsphingoidbase{<: Glycan}}, headgroup) = string(head_abbr(headgroup), "SPB")
default_lipid_classabbr(::Type{<:Glycosylsphingoidbase{<: GlyComp}}, headgroup) = string(head_abbr(headgroup), "SPB")
default_lipid_classabbr(::Type{<:SphingoidBasePhosphate}, headgroup) = "SPBP"
default_lipid_classabbr(::Type{<:Lysoinositolphosphorylceramide}, headgroup) = "LIPC"
default_lipid_classabbr(::Type{<:Lysoglycosylinositolphosphorylceramide}, headgroup) = string(head_abbr(headgroup), "LIPC")
default_lipid_classabbr(::Type{<:Lysoethanolaminephosphorylceramide}, headgroup) = "LEPC"
default_lipid_classabbr(::Type{<:Lysomannosylinositolphosphorylceramide}, headgroup) = "LMIPC"
default_lipid_classabbr(::Type{<:Lysomannosyldiinositolphosphorylceramide}, headgroup) = "LM(IP)2C"
default_lipid_classabbr(::Type{<:Lysosphingomyelin}, headgroup) = "LSM"
default_lipid_classabbr(::Type{<:Lysosulfonolipid}, headgroup) = "LSL"
function default_lipid_classabbr(::Type{<:Acylceramide}, headgroup)
    string(head_abbr(headgroup), "-ACer")
end
# function default_lipid_classabbr(c::Type{<:Acylhexosylceramide)
#     mod1, mod2 = head_abbr.(c.headgroup)
#     isempty(mod1) ? string("A", mod2, "Cer") : string(mod1, "-A", mod2, "Cer")
# end
function default_lipid_classabbr(::Type{<:Acylsphingomyelin}, headgroup)
    head = head_abbr(first(headgroup))
    isempty(head) ? "ASM" : string(head, "-ASM")
end

"""
    unknownlipidchirality(::Type{<:Lipid})

Return the unknown lipid chirality for the given lipid type.
"""
unknownlipidchirality(::Type{<:Lipid}) = RSChirality()
unknownlipidchirality(::Type{<:Radyldiglycerol}) = (RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:Radylglycerophosphoaminoacid}) = (RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:Radylglycerophosphoglycerol}) = (RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:Bisradylglycerophosphate}) = (RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:Radyltriglycerol}) = (RSChirality(), RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:Bisradylglycerophosphoglycerol}) = (RSChirality(), RSChirality(), RSChirality())
unknownlipidchirality(::Type{<:MixSphingoBone}) = RSSystem[]