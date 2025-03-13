const NAAA = [    
    "Ala", 
    "Arg", 
    "Asn", 
    "Asp", 
    "Cys", 
    "Gln", 
    "Glu", 
    "Gly", 
    "His", 
    "Ile", 
    "Leu", 
    "Lys", 
    "Met", 
    "Phe", 
    "Pro", 
    "Ser", 
    "Thr", 
    "Trp", 
    "Tyr", 
    "Val", 
    "Orn"
]
const GLYCAN_CER = ["GlcCer", "GalCer", "LacCer"]
const GLYCAN_AS_LIPID = ["GM3", "GM2", "GM1a", "GD1a", "GD1aα", "GD1aa", "GT1a", "GT1aα", "GT1aa", "GD3", "GD2", "GD1b", "GT1b", "GT1bα", "GT1ba", "GQ1b", "GQ1bα", "GQ1ba",
"GT3", "GT2", "GT1c", "GQ1c", "GQ1cα", "GQ1ca", "GP1c", "GP1cα", "GP1ca", "GA2", "GA1", "GM1b", "GM1α", "GM1alpha", "GD1c", "GD1α", "GD1alpha", "GD1e", "GM1", "GD1", "GT1", "GQ1", "GP1", 
"Gb3", "Gb4", "iGb3", "iGb4", "Lc3", "LM1", "GM4", "SM4", "SM3", "SM2", "SM1", "SM1a", "SM1b", "SB1", "SB1a", "GM1?", "GD1?", "GT1?", "GQ1?"]
const FREESTEROL = ["C", "DSMS", "DC", "CAS", "BRS", "EGS", "DEGS", "SIS", "STS"]
const CLASS = [
    "HC", "FA", "FAL", "FOH", "WE", "FAM", "FattyAmide", "FN", "FattyAmine", "NA", "NAE", "NAT", [string("NA", x) for x in NAAA]..., "CAR", "CoA", "FAHFA",
    "MG", "DG", "TG", "SQMG", "SQDG", "MGMG", "MGDG", "DGMG", "DGDG", "GlcAMG", "GlcADG",
    "PA", "LPA", "PC", "LPC", "PE", "LPE", "PE-NMe", "LPE-NMe", "PE-NMe2", "LPE-NMe2", "PE-N", "LPE-N", "PS", "LPS", "PS-N", "LPS-N", "PI", "LPI",
    "PG", "LPG", "PMeOH", "LPMeOH", "PEtOH", "LPEtOH", "GP", "LGP", "PIP", "LPIP", "PIP2", "LPIP2", "PIP3", "LPIP3", "PGP", "LPGP", "BPA", "SLBPA", "BMP", "LBPA", "CL", "MLCL", "DLCL", "GP-NAE", 
    "Cer", "SPB", "CerP", "SPBP", "EPC", "LEPC", "PE-Cer", "PE-SPB", "IPC", "LIPC", "PI-Cer", "PI-SPB", "MIPC", "LMIPC", "M(IP)2C", "LM(IP)2C", "SM", "LSM",
    GLYCAN_AS_LIPID..., [endswith(x, "Cer") ? replace(x, "Cer" => "SPB") : string("Lyso", x) for x in GLYCAN_AS_LIPID]..., "SL", "LSL", "ACer", "ASM", 
    "ST", "SE", "BA", "SG", "ASG", (FREESTEROL .* "E")..., "FC", FREESTEROL[begin + 1:end]...
]

const POST_MODIFIER = Dict{String, String}(
    "LPE-N" => "(\\([^)(]*+(?:(?3)[^)(]*)*+\\))",
    "LPS-N" => "(\\([^)(]*+(?:(?3)[^)(]*)*+\\))",
    "PE-N"  => "(\\([^)(]*+(?:(?3)[^)(]*)*+\\))",
    "PS-N"  => "(\\([^)(]*+(?:(?3)[^)(]*)*+\\))",
    "GM1?"   => "\\((.*)\\)",
    "GD1?"   => "\\((.*)\\)",
    "GT1?"   => "\\((.*)\\)",
    "GQ1?"   => "\\((.*)\\)",
    # "GP1?"   => "(\\(.*\\))",
    # "SM1?"   => "(\\(.*\\))",
)

const PRE_MODIFIER = Dict{String, String}(
    # "ASM"       => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)-",
    # "ACer"      => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)-",
    # "AHexCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)-",
    # "AGlcCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)-",
    # "AGalCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)-",
    "ASM"       => "(FA\\s+\\d+:\\d+[^-]*)-",
    "ACer"      => "(FA\\s+\\d+:\\d+[^-]*)-",
    "AHexCer"   => "(FA\\s+\\d+:\\d+[^-]*)-",
    "AGlcCer"   => "(FA\\s+\\d+:\\d+[^-]*)-",
    "AGalCer"   => "(FA\\s+\\d+:\\d+[^-]*)-",
)

const CLASS_STRUCT = Dict{String, Function}(
    "HC"    => class_struct_HC,
    "FA"    => class_struct_FA,
    "FAL"   => class_struct_FAL, 
    "FOH"   => class_struct_FOH, 
    "FAM"   => class_struct_FAM, 
    "FattyAmide" => class_struct_FAM,
    "FN"   => class_struct_FN, 
    "FattyAmine" => class_struct_FN,
    "CAR"   => class_struct_CAR,
    "CoA"   => class_struct_CoA, 
    "NAE"   => class_struct_NAE,
    "NAT"   => class_struct_NAT,
    [string("NA", x) => class_struct_NAAA for x in NAAA]...,
    "WE"    =>class_struct_WE, 
    "NA"    =>class_struct_NA, 
    "FAHFA"    =>class_struct_FAHFA, 

    "MG"    => class_struct_MG,
    "DG"    => class_struct_DG,
    "TG"    => class_struct_TG,
    "SQMG"  => class_struct_SQMG,
    "SQDG"  => class_struct_SQDG,
    "MGMG"  => class_struct_MGMG,
    "MGDG"  => class_struct_MGDG,
    "DGMG"  => class_struct_DGMG,
    "DGDG"  => class_struct_DGDG,
    "GlcAMG"    => class_struct_GlcAMG,
    "GlcADG"    => class_struct_GlcADG,

    "PA"    => class_struct_PA, 
    "LPA"   => class_struct_LPA, 
    "PC"    => class_struct_PC, 
    "LPC"   => class_struct_LPC, 
    "PE"    => class_struct_PE, 
    "LPE"   => class_struct_LPE, 
    "PE-NMe"    => class_struct_PE_NMe, 
    "LPE-NMe"   => class_struct_LPE_NMe, 
    "PE-NMe2"   => class_struct_PE_NMe2, 
    "LPE-NMe2"  => class_struct_LPE_NMe2, 
    "PE-N"  => class_struct_PE_N, 
    "LPE-N" => class_struct_LPE_N, 
    "PS"    => class_struct_PS, 
    "LPS"   => class_struct_LPS, 
    "PS-N"  => class_struct_PS_N, 
    "LPS-N" => class_struct_LPS_N, 
    "PI"    => class_struct_PI, 
    "LPI"   => class_struct_LPI, 
    "PG"    => class_struct_PG, 
    "LPG"   => class_struct_LPG, 
    "PMeOH" => class_struct_PMeOH, 
    "LPMeOH"    => class_struct_LPMeOH, 
    "PEtOH" => class_struct_PEtOH, 
    "LPEtOH"    => class_struct_LPEtOH, 
    "GP"    => class_struct_GP, 
    "LGP"   => class_struct_LGP, 
    "PIP"   => class_struct_PIP, 
    "LPIP"  => class_struct_LPIP, 
    "PIP2"  => class_struct_PIP2, 
    "LPIP2" => class_struct_LPIP2, 
    "PIP3"  => class_struct_PIP3,
    "LPIP3" => class_struct_LPIP3,
    "PGP"   => class_struct_PGP, 
    "LPGP"  => class_struct_LPGP, 
    "BPA"   => class_struct_BPA, 
    "SLBPA" => class_struct_SLBPA, 
    "BMP"   => class_struct_LBPA, 
    "LBPA"  => class_struct_LBPA, 
    "CL"    => class_struct_CL, 
    "MLCL"  => class_struct_MLCL, 
    "DLCL"  => class_struct_DLCL, 
    "GP-NAE"    => class_struct_GP_NAE, 
    
    "Cer"   => class_struct_Cer,
    "SPB"   => class_struct_SPB,
    "CerP"  => class_struct_CerP,
    "SPBP"  => class_struct_SPBP,
    "EPC"   => class_struct_EPC,
    "LEPC"  => class_struct_LEPC,
    "PE-Cer"    => class_struct_EPC,
    "PE-SPB"    => class_struct_LEPC,
    "IPC"   => class_struct_IPC,
    "LIPC"  => class_struct_LIPC,
    "PI-Cer"   => class_struct_IPC,
    "PI-SPB"  => class_struct_LIPC,
    "MIPC"  => class_struct_MIPC,
    "LMIPC" => class_struct_LMIPC,
    "M(IP)2C"   => class_struct_MIPIPC,
    "LM(IP)2C"  => class_struct_LMIPIPC,
    "SM"    => class_struct_SM,
    "LSM"   => class_struct_LSM,
    "GlcCer"    => class_struct_GlcCer,
    "GalCer"    => class_struct_GalCer,
    "LacCer"    => class_struct_LacCer,
    "GlcSPB"    => class_struct_GlcSPB,
    "GalSPB"    => class_struct_GalSPB,
    "LacSPB"    => class_struct_LacSPB,
    [sugar => class_struct_GSL for sugar in GLYCAN_AS_LIPID]...,
    [string("Lyso", sugar) => class_struct_GSL for sugar in GLYCAN_AS_LIPID]...,
    "SL"    => class_struct_SL,
    "LSL"   => class_struct_LSL,
    "ACer"  => class_struct_ACer,
    "ASM"   => class_struct_ASM,
    # AHexCer ?
    "ST"    => class_struct_ST,
    "SE"    => class_struct_SE,
    "SG"    => class_struct_SG,
    "ASG"   => class_struct_ASG,
    "FC"    => class_struct_FC,
    [cls => class_struct_SST for cls in FREESTEROL[begin + 1:end]]..., 
    [string(cls, "E") => class_struct_SSE for cls in FREESTEROL]..., 
    [string(cls, "G") => class_struct_SSG for cls in FREESTEROL]..., 
    [string("A", cls, "G") => class_struct_SASG for cls in FREESTEROL]..., 
    "BA"    => class_struct_BA
)

# mod has O or no mod/db position
# const FG_SPECIES_ONLY = Dict{String, FunctionalGroup}{

# }
const FG_SPECIES_CLINKED = Dict{String, AbstractFunctionalGroup}(
    "O"     => OxygenAtom(),
    "Br"    => Bromo(),
    "Cl"    => Chloro(),
    "F"     => Fluoro(),
    "I"     => Iodo()
)

const FG_SPECIES = Dict{String, AbstractFunctionalGroup}(
    "O"     => OxygenAtom(),
    "Br"    => Bromo(),
    "Cl"    => Chloro(),
    "F"     => Fluoro(),
    "I"     => Iodo(),
    # ASG/BA
    # FA 
    # O CHAIN  (Hydroxy, CHAIN(Dehydroxy))
    # N CHAIN  (Amino, CHAIN(Dehydroxy))
    # isotope ?
    "T"     => Tauryl(),
    [x => conjugation(parse_aa(x)) for x in NAAA]...,
    [k => Substituent(Dehydroxy, v(), Anomerposition(0x00)) for (k, v) in MONO_STRUCT]...,

    "Me"    => Methyl(),
    "Et"    => Ethyl(),
    "Fo"    => Formyl(),
    "Ac"    => Acetyl(),
    "S"     => Sulfate(),
)


const FG_CLINKED = Dict{String, AbstractFunctionalGroup}(
    "H"     => Hydrogen(),
    "Me"    => Methyl(),
    "Et"    => Ethyl(),
    "Fo"    => Formyl(),
    "Ac"    => Acetyl(),
    "CN"    => Cyano(),
    "Br"    => Bromo(),
    "Cl"    => Chloro(),
    "F"     => Fluoro(),
    "I"     => Iodo(),
    "OH"    => Hydroxy(),
    "COOH"  => CarboxylicAcidGroup(),
    "oxo"   => Oxo(),
    "oxy"   => Alkoxy(),
    "OOH"   => Hydroperoxyl(),
    "Ep"    => Epoxy(),
    "OO"    => Peroxy(),
    "NH2"   => Amino(),
    "SH"    => Sulfanyl(),
    "NO2"   => Nitro(),
)

const FG = Dict{String, AbstractFunctionalGroup}(
    "H"     => Hydrogen(),
    "Me"    => Methyl(),
    "OMe"   => Methoxy(),
    "Et"    => Ethyl(),
    "OEt"   => Ethoxy(),
    "Fo"    => Formyl(),
    "OFo"   => Oformyl(),
    "NFo"   => Nformyl(),
    "Ac"    => Acetyl(),
    "OAc"   => Oacetyl(),
    "NAc"   => Nacetyl(),
    "CN"    => Cyano(),
    "Br"    => Bromo(),
    "Cl"    => Chloro(),
    "F"     => Fluoro(),
    "I"     => Iodo(),
    "OH"    => Hydroxy(),
    "COOH"  => CarboxylicAcidGroup(),
    "oxo"   => Oxo(),
    "oxy"   => Alkoxy(),
    "OOH"   => Hydroperoxyl(),
    "Ep"    => Epoxy(),
    "OO"    => Peroxy(),
    "NH2"   => Amino(),
    "SH"    => Sulfanyl(),
    "P"     => Phosphate(),
    "MP"    => Phosphate(),
    "DP"    => Diphosphate(),
    "TP"    => Triphosphate(),
    "S"     => Sulfate(),
    "SO3"   => Sulfite(),
    "NO2"   => Nitro(),
    "COT"   => XLinkedFunctionalGroup(CarboxylicAcidGroup(), Tauryl()),
    [string("CO", x) => XLinkedFunctionalGroup(CarboxylicAcidGroup(), conjugation(parse_aa(x))) for x in NAAA]...,
    [string("O", k) => XLinkedFunctionalGroup(Hydroxy(), Substituent(Dehydroxy, v(), Anomerposition(dehydrogenposition(v())))) for (k, v) in MONO_STRUCT]...,
    # O..., on chain (Hydroxy, ...(Dehydroxy))
    # N..., on chain (Amino, ...(Dehydroxy))
    # isotope ?
)

# const MOD_STRUCT = Dict{String, }(
#     r"^;O(\d*)$"        => ,
#     r"^;((\d+)OH,?)+$"  => ,
#     r"^;\(+OH\)+(\d*)$" => ,
# )

function class_regex()
    C = collect(keys(CLASS_STRUCT))
    RC = map(C) do x                                                                                                                                   
        "^(.*?)" * get(PRE_MODIFIER, x, "()") * replace(x, "(" => "\\(", ")" => "\\)") * get(POST_MODIFIER, x, "()") * "\$"                                                                                                                                             
    end
    id = sortperm(RC; by = x -> length(x), rev = true)
    map((x, y) -> (Regex(x), y), RC[id], C[id])
    # sort!(C; by = x -> length(x), rev = true)
    # KC = join(RC, "|")
    # DC = join(OC, ")(?:")
    # Regex("^(.*?)(" * KC * ")[^(?:" * DC * ")]*(\\([\\d,\\,,',\\s]+\\))?\\s*(\\[[^\\]\\[]+\\])?[\\s,;]")
    # Regex("(" * KC * ")(\\([\\d\\,'\\s]+\\))?\\s*(\\[[^\\]\\[]+\\])?[\\s,;]")
end

function chain_regex()
    # X O,P-
    # r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^(?:sn\-)][^)(]*+(?:(?1)[^)(]*)*+\))?(?:(;O\d*)?(\[[^\[]+\])?)?((?:;[^)(;(/_]*(\([^(?:sn\-)][^)(]*+(?:(?2)[^)(]*)*+\))?[^)(;(/_]*)*)(\(sn-*\d*'*\))?$"
    r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^(?:sn\-)][^)(]*+(?:(?1)[^)(]*)*+\))?(?:(;O\d*)?(\[[^\[]+\])?)?((?:;[^)(;/_]*([\(\[][^(?:s)][^)(]*+(?:(?2)[^)(]*)*+[\)\]])?[^)(;/_]*)*)(\(sn-*\d*'*\))?$"
    #r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^(?:sn\-)][^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^)(;(/_]*(\([^(?:sn\-)][^)(]*+(?:(?2)[^)(]*)*+\))?[^)(;(/_]*)*)(\(sn-*\d*'*\))?+"
    #r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^)(;(/_]*(\([^)(]*+(?:(?2)[^)(]*)*+\))?[^)(;(/_]*)*)(\(sn-*\d*'*\))?"
end
# O/P-18:1(5Z);.....[(1,2,3,4,5)13C5](sn-1)
# r"((?:;(([^)(;/_]*[\(\[][^)(\]\[]*+(?:(?3)[^)(\]\[]*)*+[\)\]])?[^)(;(/_]*)))"
# map(x -> collect(eachmatch(r"((?:;(([^)(;/_]*[\(\[][^)(\]\[]*+(?:(?3)[^)(\]\[]*)*+[\)\]])?[^)(;(/_]*)))", x)), c)
# map(x -> match(r"^([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?", x), c)
# map(x -> match(r"(\(sn-*\d*'*\))?$", x), c)
# function clsmod_regex()
#     pre = collect(keys(PRE_MODIFIER))
#     post = collect(keys(POST_MODIFIER))
#     C = [pre..., post...]
#     unique!(C)
#     RC = map(C) do x
#         replace(x, "(" => "\\(", ")" => "\\)")
#     end
#     sort!(RC; by = x -> length(x), rev = true)
#     KC = join(RC, "|")
#     DC = join(RC, ")(?:")
#     Regex("^(.*?)(" * KC * ")[^(?:" * DC * ")]*")
# end

const REGEX = Dict{Symbol, Any}(
    :class => class_regex(),
    :chain => chain_regex()
    # :clsmod => clsmod_regex()
)

function register_class!(cls::String)
    push!(CLASS, cls)
    REGEX[:class] = class_regex()
end

const SPECIFIC_LIPID = ["Cholesterol", "Desmosterol", "Dihydrocholesterol", "Campesterol", "Brassicasterol", "Ergosterol", "Dehydroergosterol", "Sitosterol", "Stigmasterol"]
