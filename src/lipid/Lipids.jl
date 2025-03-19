module Lipids
using Reexport, IterTools, MLStyle, Dates, Downloads, ZipArchives
using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical, tuplize
@reexport using ..MassSpecBioChemicals.BasicCompounds, ..MassSpecBioChemicals.Metabolites, ..MassSpecBioChemicals.Proteins, ..MassSpecBioChemicals.Glycans
import ..MassSpecBioChemicals: parentchemical, leavinggroup, conjugation, repr_linkage, dehydroxyposition, dehydrogenposition
import MassSpecChemicals: parse_chemical, getchemicalattr
import Base: isless, getindex, length, size, eltype, getproperty, propertynames, setproperty!, keys, values, pairs
using ..MassSpecBioChemicals: lk, makechemical, makelinkage, concatchemical, AbstractFunctionalGroup, UnknownGroup, dehydrogenposition, dehydroxyposition
using ..MassSpecBioChemicals.Glycans: ap, α, β, parse_monosaccharide, MONO_STRUCT
using ..MassSpecBioChemicals.Proteins: parse_aa, parse_aa_fg, parse_aa3, letter3_abbr
export AbstractCarbonChain, CarbonChain, IsoprenoidChain, Acyl, Alkyl, Alkenyl, SPB, AbstractSTRing, STRing, SRing, DSMSRing, DCRing, CASRing, BRSRing, EGSRing, DEGSRing, SISRing, STSRing,
        Lipid, 
        FattyAcyl, MonoFattyAcyl, Hydrocarbon, FattyAcid, FattyAlcohol, FattyAldehyde, FattyAmide, FattyAmine, FattyAcylCarnitine, FattyAcylCoA, NacylAmine, FattyAcylEster, WaxEster, NacylAlkylAmine, FattyAcylEstolide,
        Glycerolipid, Monoradylglycerol, Diradylglycerol, Triradylglycerol, Estolide, 
        Omodifiedmonoradylglycerol, Sulfoquinovosylmonoradylglycerol, Monogalactosylmonoradylglycerol, Digalactosylmonoradylglycerol,
        Omodifieddiradylglycerol, Sulfoquinovosyldiradylglycerol, Monogalactosyldiradylglycerol, Digalactosyldiradylglycerol,

        Glycerophospholipid, GlycerophosphoNacylethanolamine, 
        Phosphatidicacid, Phosphatidylcholine, Phosphatidylethanolamine, PhosphatidylNmethylethanolamine, PhosphatidylNNdimethylethanolamine, Phosphatidylserine, Phosphatidylinositol, Phosphatidylglycerol, Phosphatidylmethanol, Phosphatidylethanol,
        AbstractPhosphatidylinositolphosphate, Phosphatidylinositolphosphate, Phosphatidylinositolbiphosphate, Phosphatidylinositoltriphosphate, Phosphatidylglycerolphosphate, PhosphatidylNmodifiedethanolamine, PhosphatidylNmodifiedserine,
        Lysophosphatidicacid, Lysophosphatidylcholine, Lysophosphatidylethanolamine, LysophosphatidylNmethylethanolamine, LysophosphatidylNNdimethylethanolamine, Lysophosphatidylserine, Lysophosphatidylinositol, Lysophosphatidylglycerol, Lysophosphatidylmethanol, Lysophosphatidylethanol, 
        AbstractLysophosphatidylinositolphosphate, Lysophosphatidylinositolphosphate, Lysophosphatidylinositolbiphosphate, Lysophosphatidylinositoltriphosphate, Lysophosphatidylglycerolphosphate, LysophosphatidylNmodifiedethanolamine, LysophosphatidylNmodifiedserine,
        Lysobisphosphatidicacid, Semilysobisphosphatidicacid, Bisphosphatidicacid, Dilysocardiolipin, Monolysocardiolipin, Cardiolipin,

        Sphingolipid, Ceramide, SphingoidBase, Glycosylceramide, Glycosylsphingoidbase,
        CeramidePhosphate, SphingoidBasePhosphate, Sphingomyelin, Lysosphingomyelin, Inositolphosphorylceramide, Ethanolaminephosphorylceramide, Glycosylinositolphosphorylceramide, Mannosylinositolphosphorylceramide, Mannosyldiinositolphosphorylceramide,
        Lysoinositolphosphorylceramide, Lysoethanolaminephosphorylceramide, Lysoglycosylinositolphosphorylceramide, Lysomannosylinositolphosphorylceramide, Lysomannosyldiinositolphosphorylceramide,
        Sulfonolipid, Lysosulfonolipid, Acylceramide, MixSphingoBone, Acylhexosylceramide, Acylsphingomyelin, 

        Sterol, FreeSterol, Sterylester, SubstitutedSterol,
        Prenol, Retinylester, CoenzymeQ,

        parse_lipid, chemicalsmiles_carbonchain, 

        specieslevel, molecularspecieslevel, phosphatepositionlevel, snpositionlevel, 
        passphosphatepositionlevel, passsnpositionlevel, 
        dbpositionpartiallevel, dbpositionlevel, dbconfigpartiallevel, dbconfiglevel, 
        structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel,
        structureconfigpartiallevel, structureconfiglevel,
        fullstructurelevel, completestructurelevel,

        LMSD

include("utils.jl")
abstract type AbstractCarbonChain <: FunctionalGroup{Nothing, Nothing} end
struct CarbonChain{T, D, S, I <: Union{Nothing, String}} <: AbstractCarbonChain
    carbon::UInt8
    doublebond::D
    substituent::S
    isotopiclabel::I
end
function CarbonChain{T}(carbon::UInt8, 
        doublebond::D, 
        substituent::S, 
        isotopiclabel::I = nothing
    ) where {D <: Union{UInt8, Vector{UInt8}}, S <: Union{Nothing, UInt8, <: Vector}, I, T}
    CarbonChain{T, D, S, I}(carbon, doublebond, substituent, isotopiclabel)
end

struct IsoprenoidChain{N, I <: Union{Nothing, String}} <: AbstractCarbonChain
    isotopiclabel::I
end

abstract type CarbonChainType end
abstract type Radyl <: CarbonChainType end
struct Acyl <: Radyl end
struct Alkyl <: Radyl end
struct Alkenyl <: Radyl end
struct SPB <: CarbonChainType end
abstract type AbstractSTRing <: CarbonChainType end
struct STRing <: AbstractSTRing end
struct CRing <: AbstractSTRing end # Cholesterol
struct DSMSRing <: AbstractSTRing end # Desmosterol
struct DCRing <: AbstractSTRing end # Dihydrocholesterol
struct CASRing <: AbstractSTRing end # Campesterol
struct BRSRing <: AbstractSTRing end # Brassicasterol
struct EGSRing <: AbstractSTRing end # Ergosterol
struct DEGSRing <: AbstractSTRing end # Dehydroergosterol
struct SISRing <: AbstractSTRing end # Sitosterol
struct STSRing <: AbstractSTRing end # Stigmasterol
struct BAing <: AbstractSTRing end # Bile acid

#= 
T
CarbonChainType
Tuple: order STRing, SPB, Alkenyl Alkyl Acyl, more -> few
D
UInt8: ndoublebond
Vector{UInt8}: doublebond position divrem(x, 3) 0 unknown 1 Z 2 E

F
UInt8: noxygen
Vector{Pair{FunctionalGroup, UInt8}}: FunctionalGroup => number
Vector{Pair{AbstractLinkageposition, FunctionalGroup}}: Position => FunctionalGroup

C 
CarbonChain, Tuple{<: CarbonChain}
=#
const AlkylAcylChain = Union{<: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, <: Acyl}}, <: CarbonChain{Acyl}}
const Acyl2Chain = Union{<: CarbonChain{<: Tuple{<: Acyl, <: Acyl}}, <: CarbonChain{Acyl}}
const Alkyl2AcylChain = Union{<: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, <: Union{Alkyl, Alkenyl}, <: Acyl}}, <: Tuple{<: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, <: Union{Alkyl, Alkenyl}}}, <: CarbonChain{Acyl}}, <: Tuple{<: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, Acyl}}}, <: Tuple{<: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{Acyl}}}
abstract type Lipid{B, C} <: AbstractChemical end
abstract type FattyAcyl{B, C} <: Lipid{B, C} end
struct MonoFattyAcyl{B, C} <: FattyAcyl{B, C}
    backbone::B
    chain::C
end

struct NacylAmine{B, C} <: FattyAcyl{B, C}
    backbone::B
    chain::C
end
# NAE, NAGly, NASer, NAT...
struct FattyAcylEster{B <: MonoFattyAcyl, C} <: FattyAcyl{B, C}
    backbone::B
    chain::C
    position::UInt8
end

const Hydrocarbon{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(Dihydrogen), C <: Union{<: CarbonChain{Alkyl}, <: CarbonChain{Alkenyl}}}
const FattyAcid{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: CarbonChain{Acyl}}
const FattyAldehyde{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(Dihydrogen), C <: CarbonChain{Acyl}} # FAL
const FattyAlcohol{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: Union{<: CarbonChain{Alkyl}, <: CarbonChain{Alkenyl}}} # FOH
const FattyAmide{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(Ammonia), C <: CarbonChain{Acyl}}
const FattyAmine{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(Ammonia), C <: Union{<: CarbonChain{Alkyl}, <: CarbonChain{Alkenyl}}}
const FattyAcylCarnitine{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(Carnitine), C <: CarbonChain{Acyl}}  
const FattyAcylCoA{B, C} = MonoFattyAcyl{B, C} where {B <: includeSIL(CoA), C <: CarbonChain{Acyl}}
const WaxEster{B, C} = FattyAcylEster{B, C} where {B <: FattyAlcohol, C <: AlkylAcylChain}
const NacylAlkylAmine{B, C} = NacylAmine{B, C} where {B <: FattyAmine, C <: AlkylAcylChain}
const FattyAcylEstolid{B, C} = FattyAcylEster{B, C} where {B <: FattyAcid, C <: Acyl2Chain}

const MonoradylChain = CarbonChain{<: Radyl}
const HeadMonoradylChain = Union{<: MonoradylChain, <: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}}
const DiradylChain = Union{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}
# divrem 3
const MonoDiradylChain = Union{<: MonoradylChain, <: DiradylChain}
const HeadDiradylChain = Union{<: DiradylChain, <: CarbonChain{<: Tuple{<: Radyl, <: Radyl, <: Radyl}}}
const TriradylChain = Union{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}, <: CarbonChain{<: Radyl}}, <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}
const TetraradylChain = Union{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl, <: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl, <: Radyl}}, <: CarbonChain{<: Radyl}}, 
                            <: Tuple{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}, <: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}}, 
                            <: Tuple{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}, 
                            <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}

abstract type Glycerolipid{B, C} <: Lipid{B, C} end
struct Monoradylglycerol{B <: includeSIL(Glycerol), C <: MonoradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # MG [0] [1] [2] [3]
struct Diradylglycerol{B <: includeSIL(Glycerol), C <: DiradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # DG [0] [0, 0] [1, 0] [1, 2] [1, 3]
struct Triradylglycerol{B <: includeSIL(Glycerol), C <: TriradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # TG [0] [0, 0] [0, 0, 0] [1/2/3, 0] [1/3, 0, 0] [1, 2, 3] # divrem 16 divrem 4
# struct Estolide{B <: includeSIL(Glycerol), C <: TriradylChain} <: Glycerolipid{B, C}
#     backbone::B
#     chain::C
#     sn::UInt8
# end
struct Omodifiedradylglycerol{B <: DehydratedChemical, C <: MonoDiradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # DehydratedChemical((Glycerol(), ?), [[Odehydrogen(), ?], [1, ?]]) # Omodeified @1
const Omodifiedmonoradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B, C <: MonoradylChain}
const Omodifieddiradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B, C <: DiradylChain}
const Sulfoquinovosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Sulfoquinovose), <: includeSIL(Glycerol)}}, C}
const Sulfoquinovosylmonoradylglycerol{B, C} = Sulfoquinovosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Sulfoquinovose), <: includeSIL(Glycerol)}}, C <: MonoradylChain}
const Sulfoquinovosyldiradylglycerol{B, C} = Sulfoquinovosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Sulfoquinovose), <: includeSIL(Glycerol)}}, C <: DiradylChain}
const Monogalactosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C}
const Monogalactosylmonoradylglycerol{B, C} = Monogalactosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C <: MonoradylChain}
const Monogalactosyldiradylglycerol{B, C} = Monogalactosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C <: DiradylChain}
const Digalactosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C}
const Digalactosylmonoradylglycerol{B, C} = Digalactosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C <: MonoradylChain}
const Digalactosyldiradylglycerol{B, C} = Digalactosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Galactose), <: includeSIL(Galactose), <: includeSIL(Glycerol)}}, C <: DiradylChain}
const Glucuronosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(GlucuronicAcid), <: includeSIL(Glycerol)}}, C}
const Glucuronosylmonoradylglycerol{B, C} = Glucuronosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(GlucuronicAcid), <: includeSIL(Glycerol)}}, C <: MonoradylChain}
const Glucuronosyldiradylglycerol{B, C} = Glucuronosylradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(GlucuronicAcid), <: includeSIL(Glycerol)}}, C <: DiradylChain}
# SMGDG, CDPDAG
abstract type Glycerophospholipid{B, C} <: Lipid{B, C} end # PX [0] [0, 0] [1, 2], LPX [0] [1/2]
struct Radylglycerophosphate{B <: DehydratedChemical, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end

const PPA = includeSIL(PhosphoricAcid)
glycerophospho(T) = DehydratedChemical{<: Tuple{<: Union{<: T, <: IsotopiclabeledChemical{<: T}}, <: PPA, <: includeSIL(Glycerol)}}
glycerophospho(T, S) = DehydratedChemical{<: Tuple{<: Union{<: T, <: IsotopiclabeledChemical{<: T}}, <: Union{<: S, <: IsotopiclabeledChemical{<: S}}, <: PPA, <: includeSIL(Glycerol)}}
glycerophosphophosphate(T) = DehydratedChemical{<: Tuple{<: PPA, <: Union{<: T, <: IsotopiclabeledChemical{<: T}}, <: PPA, <: includeSIL(Glycerol)}}
const Monoradylglycerophosphate{B, C} = Radylglycerophosphate{B, C} where {B, C <: HeadMonoradylChain}
const Diradylglycerophosphate{B, C} = Radylglycerophosphate{B, C} where {B, C <: HeadDiradylChain}
const Lysophosphatidicacid{B, C} = Monoradylglycerophosphate{B, C} where {B <: DehydratedChemical{<: Tuple{<: PPA, <: includeSIL(Glycerol)}}, C <: MonoradylChain}
const Phosphatidicacid{B, C} = Diradylglycerophosphate{B, C} where {B <: DehydratedChemical{<: Tuple{<: PPA, <: includeSIL(Glycerol)}}, C <: DiradylChain}
const Lysophosphatidylcholine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Choline), C <: MonoradylChain}
const Phosphatidylcholine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Choline), C <: DiradylChain}
const Lysophosphatidylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanolamine), C <: MonoradylChain}
const Phosphatidylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanolamine), C <: DiradylChain}
const LysophosphatidylNmodifiedethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {D, B <: glycerophospho(D, Ethanolamine), C <: HeadMonoradylChain}
const PhosphatidylNmodifiedethanolamine{B, C} = Diradylglycerophosphate{B, C} where {D, B <: glycerophospho(D, Ethanolamine), C <: HeadDiradylChain}
const LysophosphatidylNmethylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Nmethylethanolamine), C <: MonoradylChain}
const PhosphatidylNmethylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Nmethylethanolamine), C <: DiradylChain}
const LysophosphatidylNNdimethylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(NNdimethylethanolamine), C <: MonoradylChain}
const PhosphatidylNNdimethylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(NNdimethylethanolamine), C <: DiradylChain}
const Lysophosphatidylserine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Serine), C <: MonoradylChain}
const Phosphatidylserine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Serine), C <: DiradylChain}
const LysophosphatidylNmodifiedserine{B, C} = Monoradylglycerophosphate{B, C} where {D, B <: glycerophospho(D, Serine), C <: HeadMonoradylChain}
const PhosphatidylNmodifiedserine{B, C} = Diradylglycerophosphate{B, C} where {D, B <: glycerophospho(D, Serine), C <: HeadDiradylChain}
const Lysophosphatidylinositol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Inositol), C <: MonoradylChain}
const Phosphatidylinositol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Inositol), C <: DiradylChain} # include PIP
const Lysophosphatidylglycerol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Glycerol), C <: MonoradylChain}
const Phosphatidylglycerol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Glycerol), C <: DiradylChain}
const Lysophosphatidylglycerolphosphate{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophosphophosphate(Glycerol), C <: MonoradylChain}
const Phosphatidylglycerolphosphate{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophosphophosphate(Glycerol), C <: DiradylChain}
const Lysophosphatidylmethanol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Methanol), C <: MonoradylChain}
const Phosphatidylmethanol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Methanol), C <: DiradylChain}
const Lysophosphatidylethanol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanol), C <: MonoradylChain}
const Phosphatidylethanol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanol), C <: DiradylChain}

struct Bisradylglycerophosphate{B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol)}}, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt16
end
const Bisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: TetraradylChain}
# BPA [0] [0, 0] [0, 0, 0] [0, 0, 0, 0] [1/2, 0] [1, PXP, 3] [1, PXP, 3, 3] [1, PXP, 4, 3] [1, 1, PXP, 3, 3] [1, 1, PXP, 4, 3] [1, 2, PXP, 4, 3]
const Semilysobisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: TriradylChain}
# SLBPA [0] [0, 0, 0] [1, 0] [1, P, 3/4] [1/2, 0, P, 0] [1, 1, P, 0] [1, 1, P, 3/4] [1, 2, P, 0] [1, 2, P, 3/4]  
const Lysobisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: DiradylChain} 
# BMP/LBPA [0] [0, 0] [1/2, P, 0] [0, P, 1/2] [1/2, P, 1/2]

struct Bisradylglycerophosphoglycerol{B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol)}}, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt16
end 
const Cardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: TetraradylChain}
# CL [0] [0, 0] [0, 0, 0] [0, 0, 0, 0] [1/2, 0] [1, PXP, 3] [1, PXP, 3, 3] [1, PXP, 4, 3] [1, 1, PXP, 3, 3] [1, 1, PXP, 4, 3] [1, 2, PXP, 4, 3]
const Monolysocardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: TriradylChain} 
# MLCL [0] [0, 0, 0] [1, 0] [1, PXP, 3/4] [1/2, 0, PXP, 0] [1, 1, PXP, 0] [1, 1, PXP, 3/4] [1, 2, PXP, 0] [1, 2, PXP, 3/4]
const Dilysocardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: DiradylChain} 
# DLCL [0] [0, 0] [1/2, PXP, 0] [0, PXP, 1/2] [1/2, PXP, 1/2]

struct GlycerophosphoNacylethanolamine{B <: glycerophospho(Ethanolamine), C <: CarbonChain{Acyl}} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
end # GP-NAE
# PnC, PnE, PPA

abstract type Sphingolipid{B, C} <: Lipid{B, C} end
const SUM_CER = Union{<: CarbonChain{<: Tuple{SPB, Acyl}}, <: Tuple{<: CarbonChain{SPB}, <: CarbonChain{Acyl}}}
const SUM_ACYLCER = Union{<: CarbonChain{<: Tuple{SPB, Acyl, Acyl}}, <: Tuple{<: CarbonChain{SPB}, <: CarbonChain{<: Tuple{Acyl, Acyl}}}, <: Tuple{<: CarbonChain{<: Tuple{SPB, Acyl}}, <: CarbonChain{Acyl}}}
struct SphingoBone{H, C <: Union{<: CarbonChain{SPB}, <: SUM_CER, <: SUM_ACYLCER}} <: Sphingolipid{H, C}
    headgroup::H
    chain::C
    position::UInt8
end
# HexCer/IPC/SM... headgroup position
# CerP/SPBP.. headgroup position < 32, divrem(x, 32) 35 => [1, 3]
# ACer [0] [0(SPB), 0(ACYL-ACYL)] [0(SPB-ACYL), 0(ACYL)] [2(SPB-ACYL), 0(ACYL)] [1~(SPB-ACYL), 2(ACYL)] [0(SPB), 2(ACYL), 1~(ACYL)] 
# divrem(x, 3) = 
# (0, 2) = [2(SPB-ACYL), 0(ACYL)]
# (a, 2) = [a(SPB-ACYL), 2(ACYL)]/[0(SPB), 2(ACYL), a(ACYL)] 

struct MixSphingoBone{H, C <: Union{<: CarbonChain{SPB}, <: SUM_CER, <: SUM_ACYLCER}} <: Sphingolipid{H, C}
    headgroup::H
    chain::C
    position::Vector{UInt8}
end
# > 1 head group
# pos = (hg1, ..., hg2)
# Hex-(FA....-)ACer/FA....-ASM 
# P-[Hex-]Cer(4, 1) 18:0;3OH/18:0 = Hex-[P-]Cer(1, 4) 18:0;3OH/18:0
# P-[P-][Hex-]Cer(3, 4, 1) 18:0/18:0 = Hex-[P-][P-]Cer(1, 3, 4) 18:0/18:0 = P-[Hex-][P-]Cer(3, 1, 4) 18:0/18:0 = 
# Hex-[P-][FA 18:0-]ACer(1, 4, 3) 18:0/18:0

const CeramideBone{H, C} = SphingoBone{H, C} where {H, C <: SUM_CER}
const SphingoidBaseBone{H, C} = SphingoBone{H, C} where {H, C <: CarbonChain{SPB}}
const Ceramide{C} = CeramideBone{Nothing, C} where {C <: SUM_CER}
const SphingoidBase{C} = SphingoidBaseBone{Nothing, C} where {C <: CarbonChain{SPB}}
const Glycosylceramide{H, C} = CeramideBone{H, C} where {H <: AbstractGlycan, C <: SUM_CER}
const Glycosylsphingoidbase{H, C} = SphingoidBaseBone{H, C} where {H <: AbstractGlycan, C <: CarbonChain{SPB}}

const Hexlike = includeSIL(Glycan{<: Tuple{<: AbstractHexose}})

const CeramidePhosphate{H, C} = CeramideBone{H, C} where {H <: PPA, C <: SUM_CER}
const SphingoidBasePhosphate{H, C} = SphingoidBaseBone{H, C} where {H <: PPA, C <: CarbonChain{SPB}}
const Inositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Inositol), <: PPA}}, C <: SUM_CER}
const Lysoinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Inositol), <: PPA}}, C <: CarbonChain{SPB}}
const Ethanolaminephosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Ethanolamine), <: PPA}}, C <: SUM_CER}
const Lysoethanolaminephosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Ethanolamine), <: PPA}}, C <: CarbonChain{SPB}}
const Mannosylinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Mannose), <: includeSIL(Inositol), <: PPA}}, C <: SUM_CER}
const Lysomannosylinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Mannose), <: includeSIL(Inositol), <: PPA}}, C <: CarbonChain{SPB}}
const Mannosyldiinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Inositol), <: PPA, <: includeSIL(Mannose), <: includeSIL(Inositol), <: PPA}}, C <: SUM_CER}
const Lysomannosyldiinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Inositol), <: PPA, <: includeSIL(Mannose), <: includeSIL(Inositol), <: PPA}}, C <: CarbonChain{SPB}}
const Glycosylinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: AbstractGlycan, <: includeSIL(Inositol), <: PPA}}, C <: SUM_CER}
const Lysoglycosylinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: AbstractGlycan, <: includeSIL(Inositol), <: PPA}}, C <: CarbonChain{SPB}}
const Sphingomyelin{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Choline), <: PPA}}, C <: SUM_CER}
const Lysosphingomyelin{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: includeSIL(Choline), <: PPA}}, C <: CarbonChain{SPB}}

const Sulfonolipid{H, C} = CeramideBone{H, C} where {H <: includeSIL(SulfurousAcid), C <: SUM_CER}
const Lysosulfonolipid{H, C} = SphingoidBaseBone{H, C} where {H <: includeSIL(SulfurousAcid), C <: CarbonChain{SPB}}

const Acylceramide{C} = SphingoBone{H, C} where {H <: FattyAcid, C <: Union{<: SUM_CER, SUM_ACYLCER}}
const Acylhexosylceramide{H, C} = MixSphingoBone{H, C} where {H <: Tuple{<: FattyAcid, <: Hexlike}, C <: Union{<: SUM_CER, SUM_ACYLCER}}
const Acylsphingomyelin{H, C} = MixSphingoBone{H, C} where {H <: Tuple{<: FattyAcid, DehydratedChemical{<: Tuple{<: includeSIL(Choline), <: PPA}}}, C <: Union{<: SUM_CER, SUM_ACYLCER}}

abstract type Sterol{C} <: Lipid{Nothing, C} end
struct SterolBone{C} <: Sterol{C}
    chain::C
end
const FreeSterol{C} = SterolBone{C} where {C <: CarbonChain{<: AbstractSTRing}}
const Sterylester{C} = SterolBone{C} where {C <: Tuple{<: CarbonChain{<: AbstractSTRing}, <: CarbonChain{Acyl}}}

struct SubstitutedSterol{C, S, T} <: Sterol{C}
    backbone::SubstitutedChemical{FreeSterol{C}, S, T}
end
# ST;..., BA;..., SG, ASG

abstract type Prenol{C} <: Lipid{Nothing, C} end
struct Retinylester{C <: CarbonChain{<: Acyl}} <: Prenol{C}
    chain::C
end
struct CoenzymeQ{C <: IsoprenoidChain} <: Prenol{C} end

include("database.jl")
include("interface.jl")
include("annotationlevel.jl")
include("transform.jl")
include(joinpath("input", "input.jl"))
include("const.jl")
include("abbr.jl")
include("output.jl")
 
end