module BasicCompounds
using Reexport
using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
import MassSpecChemicals: getchemicalattr
import ..MassSpecBioChemicals: dehydroxyposition, dehydrogenposition, leavinggroup, leavinggroupelements, dehydrogengroup, dehydroxygroup, dehydroxyfunctionalgroup, isdissociated, nlinkage, nbridge
export BasicCompound, 
       isdissociated,
       Dihydrogen,
       Hydrogen,
       Ethane, 
       Ethyl, # Et
       Ethanol, 
       Ethoxy, # OEt
       Methane, 
       Methyl, # Me
       Methanol, 
       Methoxy, # OMe
       FormicAcid, 
       Formyl, # Fo
       Oformyl, # OFo
       Formamide, 
       Nformyl, # NFo
       AceticAcid, 
       Acetyl, # Ac
       Oacetyl, # OAc
       Acetamide, 
       Nacetyl, # NAc

       HydrogenBromide, 
       Bromo, # Br
       HydrogenChloride, 
       Chloro, # Cl
       HydrogenFluoride, 
       Fluoro, # F
       HydrogenIodide, 
       Iodo, # I
       HydrogenOxide, 
       OxygenAtom, 
       Hydroxy, # OH
       OLinkage, # O...
       CarboxylicAcidGroup, # COOH
       CarboxylicLinkage, # CO...
       Oxo, # oxo 2
       Ether, # oxy 2
       HydrogenPeroxide, 
       Hydroperoxyl, # OOH
       Epoxy, # Ep 2
       Peroxy, # OO 2

       Ammonia, 
       Amino, # NH2/N
       MethylAmine, 
       MethylAmino,
       DimethylAmine,
       DimethylAmino,
       TrimethylAmine,
       TrimethylAmino,
       NLinkage, # N...

       HydrogenSulfide, 
       Sulfanyl, # SH
       SLinkgae, 
	   Methanethiol,
	   MethanethiolGroup,
       HydrogenCyanide, 
       Cyano, # CN
       PhosphoricAcid, # P
       Phosphoryl, # P/MP
       DiphosphoricAcid, 
       Diphosphoryl, # DP
       TriphosphoricAcid, 
       Triphosphoryl, # TP
       SulfuricAcid, # S
       Sulfo, # S
       # Sulfate, 
       SulfonicAcid, 
       NitricAcid, 
       Nitro, # NO2
       Phenolhydroxy,
       Amide,
       Ester,
       Thioester




abstract type BasicCompound <: AbstractChemical end 

struct Dihydrogen <: BasicCompound end 
struct Hydrogen <: FunctionalGroup{Dihydrogen, Dehydrogen} end # H
struct Ethane <: BasicCompound end 
struct Ethyl <: FunctionalGroup{Ethane, Dehydrogen} end # Et
struct Ethanol <: BasicCompound end 
struct Ethoxy <: FunctionalGroup{Ethanol, Dehydrogen} end # OEt
struct Methane <: BasicCompound end 
struct Methyl <: FunctionalGroup{Methane, Dehydrogen} end # Me
struct Methanol <: BasicCompound end 
struct Methoxy <: FunctionalGroup{Methanol, Dehydrogen} end # OMe
struct FormicAcid <: BasicCompound end 
struct Formyl <: FunctionalGroup{FormicAcid, Dehydroxy} end # Fo
struct Oformyl <: FunctionalGroup{FormicAcid, Dehydrogen} end # OFo
struct Formamide <: BasicCompound end 
struct Nformyl <: FunctionalGroup{Formamide, Dehydrogen} end # NFo
struct AceticAcid <: BasicCompound end 
struct Acetyl <: FunctionalGroup{AceticAcid, Dehydroxy} end # Ac
struct Oacetyl <: FunctionalGroup{AceticAcid, Dehydrogen} end # OAc
struct Acetamide <: BasicCompound end 
struct Nacetyl <: FunctionalGroup{Acetamide, Dehydrogen} end # NAc
struct HydrogenBromide <: BasicCompound end
struct Bromo <: FunctionalGroup{HydrogenBromide, Dehydrogen} end # Br
struct HydrogenChloride <: BasicCompound end
struct Chloro <: FunctionalGroup{HydrogenChloride, Dehydrogen} end # Cl
struct HydrogenFluoride <: BasicCompound end
struct Fluoro <: FunctionalGroup{HydrogenFluoride, Dehydrogen} end # F
struct HydrogenIodide <: BasicCompound end
struct Iodo <: FunctionalGroup{HydrogenIodide, Dehydrogen} end # I
struct HydrogenOxide <: BasicCompound end
struct OxygenAtom <: UnknownGroup{HydrogenOxide, Dehydrogen} end 
struct Hydroxy <: FunctionalGroup{HydrogenOxide, Dehydrogen} end # OH
struct OLinkage <: FunctionalGroup{Hydroxy, Dehydrogen} end # O...
struct CarboxylicAcidGroup <: FunctionalGroup{FormicAcid, Demethine} end # COOH
struct CarboxylicLinkage <: FunctionalGroup{CarboxylicAcidGroup, Dehydroxy} end # CO...
struct Oxo <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxo 2
struct Ether <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxy 2
struct HydrogenPeroxide <: BasicCompound end
struct Hydroperoxyl <: FunctionalGroup{HydrogenPeroxide, Dehydrogen} end # OOH
struct Epoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # Ep 2
struct Peroxy <: FunctionalGroup{HydrogenPeroxide, Didehydrogen} end # OO 2
struct Ammonia <: BasicCompound end
struct Amino <: FunctionalGroup{Ammonia, Dehydrogen} end # NH2/N
struct NLinkage <: FunctionalGroup{Amino, Dehydrogen} end # N...
struct MethylAmine <: BasicCompound end
struct MethylAmino <: FunctionalGroup{MethylAmine, Dehydrogen} end 
struct DimethylAmine <: BasicCompound end
struct DimethylAmino <: FunctionalGroup{DimethylAmine, Dehydrogen} end 
struct TrimethylAmine <: BasicCompound end
struct TrimethylAmino <: FunctionalGroup{TrimethylAmine, Dehydrogen} end 

struct HydrogenSulfide <: BasicCompound end
struct Sulfanyl <: FunctionalGroup{HydrogenSulfide, Dehydrogen} end # SH
struct Methanethiol <: BasicCompound end 
struct MethanethiolGroup <: FunctionalGroup{Methanethiol, Dehydrogen} end
struct HydrogenCyanide <: BasicCompound end
struct Cyano <: FunctionalGroup{HydrogenCyanide, Dehydrogen} end # CN
struct PhosphoricAcid <: BasicCompound end # P
struct Phosphoryl <: FunctionalGroup{PhosphoricAcid, Dehydroxy} end # P/MP
struct DiphosphoricAcid <: BasicCompound end
struct Diphosphoryl <: FunctionalGroup{DiphosphoricAcid, Dehydroxy} end # DP
struct TriphosphoricAcid <: BasicCompound end
struct Triphosphoryl <: FunctionalGroup{TriphosphoricAcid, Dehydroxy} end # TP
struct SulfuricAcid <: BasicCompound end # S
struct Sulfo <: FunctionalGroup{SulfuricAcid, Dehydroxy} end # S
struct SLinkage <: FunctionalGroup{Sulfo, Dehydrogen} end # S...

# struct Sulfate <: FunctionalGroup{SulfuricAcid, Dehydrogen} end 
struct NitricAcid <: BasicCompound end
struct Nitro <: FunctionalGroup{NitricAcid, Dehydroxy} end # NO2

struct Phenolhydroxy <: AbstractFunctionalGroup end
struct Amide <: AbstractFunctionalGroup end
struct Ester <: AbstractFunctionalGroup end
struct Thioester <: AbstractFunctionalGroup end

getchemicalattr(::Dihydrogen, ::Val{:name}; kwargs...) = "H2"
getchemicalattr(::Ethane, ::Val{:name}; kwargs...) = "Ethane"
getchemicalattr(::Ethanol, ::Val{:name}; kwargs...) = "EtOH"
getchemicalattr(::Methane, ::Val{:name}; kwargs...) = "Methane" 
getchemicalattr(::Methanol, ::Val{:name}; kwargs...) = "MeOH"
getchemicalattr(::FormicAcid, ::Val{:name}; kwargs...) = "Formic acid" 
getchemicalattr(::Formamide, ::Val{:name}; kwargs...) = "Formamide"
getchemicalattr(::AceticAcid, ::Val{:name}; kwargs...) = "Acetic acid"
getchemicalattr(::Acetamide, ::Val{:name}; kwargs...) = "Acetamide"
getchemicalattr(::HydrogenBromide, ::Val{:name}; kwargs...) = "Hydrogen bromide"
getchemicalattr(::HydrogenChloride, ::Val{:name}; kwargs...) = "Hydrogen chloride"
getchemicalattr(::HydrogenFluoride, ::Val{:name}; kwargs...) = "Hydrogen fluoride"
getchemicalattr(::HydrogenIodide, ::Val{:name}; kwargs...) = "Hydrogen iodide"
getchemicalattr(::HydrogenOxide, ::Val{:name}; kwargs...) = "Hydrogen oxide"
getchemicalattr(::HydrogenPeroxide, ::Val{:name}; kwargs...) = "Hydrogen peroxide"
getchemicalattr(::Ammonia, ::Val{:name}; kwargs...) = "Ammonia"
getchemicalattr(::MethylAmine, ::Val{:name}; kwargs...) = "Methylamine"
getchemicalattr(::DimethylAmine, ::Val{:name}; kwargs...) = "Dimethylamine"
getchemicalattr(::TrimethylAmine, ::Val{:name}; kwargs...) = "Trimethylamine"
getchemicalattr(::HydrogenSulfide, ::Val{:name}; kwargs...) = "Hydrogen sulfide"
getchemicalattr(::Methanethiol, ::Val{:name}; kwargs...) = "Methanethiol"
getchemicalattr(::HydrogenCyanide, ::Val{:name}; kwargs...) = "Hydrogen cyanide"
getchemicalattr(::PhosphoricAcid, ::Val{:name}; kwargs...) = "Phosphoric acid"
getchemicalattr(::DiphosphoricAcid, ::Val{:name}; kwargs...) = "Diphosphoric acid"
getchemicalattr(::TriphosphoricAcid, ::Val{:name}; kwargs...) = "Triphosphoric acid" 
getchemicalattr(::SulfuricAcid, ::Val{:name}; kwargs...) = "Sulfuric acid"
getchemicalattr(::NitricAcid, ::Val{:name}; kwargs...) = "Nitric acid"

getchemicalattr(::Dihydrogen, ::Val{:elements}; kwargs...) = ["H" => 2]
getchemicalattr(::Ethane, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3]
getchemicalattr(::Ethanol, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "H" => 2, "O" => 1, "H" => 1]
getchemicalattr(::Methane, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 4] 
getchemicalattr(::Methanol, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "O" => 1, "H" => 1]
getchemicalattr(::FormicAcid, ::Val{:elements}; kwargs...) = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1] 
getchemicalattr(::Formamide, ::Val{:elements}; kwargs...) = ["H" => 1, "C" => 1, "O" => 1, "N" => 1, "H" => 2]
getchemicalattr(::AceticAcid, ::Val{:elements}; kwargs...) =["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Acetamide, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "N" => 1, "H" => 2]
getchemicalattr(::HydrogenBromide, ::Val{:elements}; kwargs...) = ["H" => 1, "Br" => 1]
getchemicalattr(::HydrogenChloride, ::Val{:elements}; kwargs...) = ["H" => 1, "Cl" => 1]
getchemicalattr(::HydrogenFluoride, ::Val{:elements}; kwargs...) = ["H" => 1, "F" => 1]
getchemicalattr(::HydrogenIodide, ::Val{:elements}; kwargs...) = ["H" => 1, "I" => 1]
getchemicalattr(::HydrogenOxide, ::Val{:elements}; kwargs...) = ["H" => 1, "O" => 1, "H" => 1]
getchemicalattr(::HydrogenPeroxide, ::Val{:elements}; kwargs...) = ["H" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Ammonia, ::Val{:elements}; kwargs...) = ["N" => 1, "H" => 3]
getchemicalattr(::MethylAmine, ::Val{:elements}; kwargs...) = ["N" => 1, "H" => 2, "C" => 1, "H" => 3]
getchemicalattr(::DimethylAmine, ::Val{:elements}; kwargs...) = ["N" => 1, "H" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 3]
getchemicalattr(::TrimethylAmine, ::Val{:elements}; kwargs...) = ["N" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3]
getchemicalattr(::Ammonia, ::Val{:elements}; kwargs...) = ["N" => 1, "H" => 3]
getchemicalattr(::HydrogenSulfide, ::Val{:elements}; kwargs...) = ["H" => 2, "S" => 1]
getchemicalattr(::Methanethiol, ::Val{:elements}; kwargs...) = ["H" => 1, "S" => 1, "C" => 1, "H" => 3]
getchemicalattr(::HydrogenCyanide, ::Val{:elements}; kwargs...) = ["H" => 1, "C" => 1, "N" => 1]
getchemicalattr(::PhosphoricAcid, ::Val{:elements}; kwargs...) = ["H" => 3, "P" => 1, "O" => 4]
getchemicalattr(::DiphosphoricAcid, ::Val{:elements}; kwargs...) = ["H" => 4, "P" => 2, "O" => 7]
getchemicalattr(::TriphosphoricAcid, ::Val{:elements}; kwargs...) = ["H" => 5, "P" => 3, "O" => 10]
getchemicalattr(::SulfuricAcid, ::Val{:elements}; kwargs...) = ["H" => 2, "S" => 1, "O" => 4]
getchemicalattr(::NitricAcid, ::Val{:elements}; kwargs...) = ["H" => 1, "N" => 1, "O" => 3]
getchemicalattr(m::BasicCompound, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)

getchemicalattr(::Dihydrogen, ::Val{:SMILES}; kwargs...) = "[H2]"
getchemicalattr(::Ethane, ::Val{:SMILES}; kwargs...) = "CC"
getchemicalattr(::Ethanol, ::Val{:SMILES}; kwargs...) = "OCC"
getchemicalattr(::Methane, ::Val{:SMILES}; kwargs...) = "C" 
getchemicalattr(::Methanol, ::Val{:SMILES}; kwargs...) = "OC"
getchemicalattr(::FormicAcid, ::Val{:SMILES}; kwargs...) = "OC(=O)" 
getchemicalattr(::Formamide, ::Val{:SMILES}; kwargs...) = "NC(=O)"
getchemicalattr(::AceticAcid, ::Val{:SMILES}; kwargs...) = "OC(=O)C"
getchemicalattr(::Acetamide, ::Val{:SMILES}; kwargs...) = "NC(=O)C"
getchemicalattr(::HydrogenBromide, ::Val{:SMILES}; kwargs...) = "Br"
getchemicalattr(::HydrogenChloride, ::Val{:SMILES}; kwargs...) = "Cl"
getchemicalattr(::HydrogenFluoride, ::Val{:SMILES}; kwargs...) = "F"
getchemicalattr(::HydrogenIodide, ::Val{:SMILES}; kwargs...) = "I"
getchemicalattr(::HydrogenOxide, ::Val{:SMILES}; kwargs...) = "O"
getchemicalattr(::HydrogenPeroxide, ::Val{:SMILES}; kwargs...) = "OO"
getchemicalattr(::Ammonia, ::Val{:SMILES}; kwargs...) = "N"
getchemicalattr(::MethylAmine, ::Val{:SMILES}; kwargs...) = "NC"
getchemicalattr(::DimethylAmine, ::Val{:SMILES}; kwargs...) = "NC(C)"
getchemicalattr(::TrimethylAmine, ::Val{:SMILES}; kwargs...) = "NC(C)(C)"
getchemicalattr(::HydrogenSulfide, ::Val{:SMILES}; kwargs...) = "S"
getchemicalattr(::Methanethiol, ::Val{:SMILES}; kwargs...) = "SC"
getchemicalattr(::HydrogenCyanide, ::Val{:SMILES}; kwargs...) = "C#N"
getchemicalattr(::PhosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)O"
getchemicalattr(::DiphosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)OP(=O)(O)O"
getchemicalattr(::TriphosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)OP(=O)(O)OP(=O)(O)O" 
getchemicalattr(::SulfuricAcid, ::Val{:SMILES}; kwargs...) = "OS(=O)(=O)O"
getchemicalattr(::NitricAcid, ::Val{:SMILES}; kwargs...) = "O[N+](=O)[O-]"

getchemicalattr(::Hydrogen, ::Val{:abbreviation}; kwargs...) = "H"
getchemicalattr(::Ethyl, ::Val{:abbreviation}; kwargs...) = "Et"
getchemicalattr(::Ethoxy, ::Val{:abbreviation}; kwargs...) = "OEt"
getchemicalattr(::Methyl, ::Val{:abbreviation}; kwargs...) = "Me"
getchemicalattr(::Methoxy, ::Val{:abbreviation}; kwargs...) = "OMe"
getchemicalattr(::Formyl, ::Val{:abbreviation}; kwargs...) = "Fo"
getchemicalattr(::Oformyl, ::Val{:abbreviation}; kwargs...) = "OFo"
getchemicalattr(::Nformyl, ::Val{:abbreviation}; kwargs...) = "NFo"
getchemicalattr(::Acetyl, ::Val{:abbreviation}; kwargs...) = "Ac"
getchemicalattr(::Oacetyl, ::Val{:abbreviation}; kwargs...) = "OAc"
getchemicalattr(::Nacetyl, ::Val{:abbreviation}; kwargs...) = "NAc"
getchemicalattr(::Bromo, ::Val{:abbreviation}; kwargs...) = "Br"
getchemicalattr(::Chloro, ::Val{:abbreviation}; kwargs...) = "Cl"
getchemicalattr(::Fluoro, ::Val{:abbreviation}; kwargs...) = "F"
getchemicalattr(::Iodo, ::Val{:abbreviation}; kwargs...) = "I"
getchemicalattr(::OxygenAtom, ::Val{:abbreviation}; kwargs...) = "O"
getchemicalattr(::Hydroxy, ::Val{:abbreviation}; kwargs...) = "OH"
getchemicalattr(::OLinkage, ::Val{:abbreviation}; kwargs...) = "O"
getchemicalattr(::CarboxylicAcidGroup, ::Val{:abbreviation}; kwargs...) = "COOH"
getchemicalattr(::CarboxylicLinkage, ::Val{:abbreviation}; kwargs...) = "CO"
getchemicalattr(::Oxo, ::Val{:abbreviation}; kwargs...) = "oxo"
getchemicalattr(::Ether, ::Val{:abbreviation}; kwargs...) = "oxy"
getchemicalattr(::Hydroperoxyl, ::Val{:abbreviation}; kwargs...) = "OOH"
getchemicalattr(::Epoxy, ::Val{:abbreviation}; kwargs...) = "Ep"
getchemicalattr(::Peroxy, ::Val{:abbreviation}; kwargs...) = "OO"
getchemicalattr(::Amino, ::Val{:abbreviation}; kwargs...) = "NH2"
getchemicalattr(::MethylAmino, ::Val{:abbreviation}; kwargs...) = "NMe"
getchemicalattr(::DimethylAmino, ::Val{:abbreviation}; kwargs...) = "NMe2"
getchemicalattr(::TrimethylAmino, ::Val{:abbreviation}; kwargs...) = "NMe3"
getchemicalattr(::NLinkage, ::Val{:abbreviation}; kwargs...) = "N"
getchemicalattr(::Sulfanyl, ::Val{:abbreviation}; kwargs...) = "SH"
getchemicalattr(::SLinkage, ::Val{:abbreviation}; kwargs...) = "S"
getchemicalattr(::MethanethiolGroup, ::Val{:abbreviation}; kwargs...) = "SMe"
getchemicalattr(::Cyano, ::Val{:abbreviation}; kwargs...) = "CN"
getchemicalattr(::Phosphoryl, ::Val{:abbreviation}; kwargs...) = "P"
getchemicalattr(::Substituent{Phosphoryl, Dehydroxy}, ::Val{:abbreviation}; kwargs...) = "P"
getchemicalattr(::Diphosphoryl, ::Val{:abbreviation}; kwargs...) = "DP"
getchemicalattr(::Triphosphoryl, ::Val{:abbreviation}; kwargs...) = "TP"
getchemicalattr(::Sulfo, ::Val{:abbreviation}; kwargs...) = "S"
# getchemicalattr(::Sulfate, ::Val{:abbreviation}; kwargs...) = "OSO3"
getchemicalattr(::Nitro, ::Val{:abbreviation}; kwargs...) = "NO2"
getchemicalattr(::Amide, ::Val{:abbreviation}; kwargs...) = "CONH"
getchemicalattr(::Ester, ::Val{:abbreviation}; kwargs...) = "COO"
getchemicalattr(::Thioester, ::Val{:abbreviation}; kwargs...) = "COS"
getchemicalattr(x::T, ::Val{:name}; group = true, kwargs...) where {T <: FunctionalGroup{<: BasicCompound}} = group ? string(T, " Group") : string(T)
function getchemicalattr(x::T, ::Val{:elements}; kwargs...) where {T <: FunctionalGroup{<: BasicCompound}} 
	es = chemicalelements(parentchemical(x))
	ls = leavinggroupelements(leavinggroup(x))
	fn_f = length(ls) > 1  ? findlast : findfirst
	for (e, n) in ls 
		if n > 0
			i = fn_f(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = fn_f(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	es
end
getchemicalattr(x::T, ::Val{:formula}; unique = false, kwargs...) where {T <: FunctionalGroup{<: BasicCompound}} = chemicalformula(chemicalelements(x); unique, kwargs...)

getchemicalattr(x::OxygenAtom, ::Val{:elements}; kwargs...) = ["O" => 1]
getchemicalattr(x::OxygenAtom, ::Val{:formula}; kwargs...) = "O"
getchemicalattr(x::OLinkage, ::Val{:elements}; kwargs...) = ["O" => 1]
getchemicalattr(x::OLinkage, ::Val{:formula}; kwargs...) = "O"
getchemicalattr(x::NLinkage, ::Val{:elements}; kwargs...) = ["N" => 1, "H" => 1]
getchemicalattr(x::NLinkage, ::Val{:formula}; kwargs...) = "NH"
getchemicalattr(x::SLinkage, ::Val{:elements}; kwargs...) = ["S" => 1]
getchemicalattr(x::SLinkage, ::Val{:formula}; kwargs...) = "S"
getchemicalattr(x::CarboxylicAcidGroup, ::Val{:elements}; kwargs...) = ["C" => -1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(x::CarboxylicAcidGroup, ::Val{:formula}; kwargs...) = "OOH"
getchemicalattr(x::CarboxylicLinkage, ::Val{:elements}; kwargs...) = ["C" => -1, "C" => 1, "O" => 1]
getchemicalattr(x::CarboxylicLinkage, ::Val{:formula}; kwargs...) = "O"

getchemicalattr(::Hydrogen, ::Val{:SMILES}; kwargs...) = "(H)"
getchemicalattr(::Ethyl, ::Val{:SMILES}; kwargs...) = "(CC)"
getchemicalattr(::Ethoxy, ::Val{:SMILES}; kwargs...) = "(OCC)"
getchemicalattr(::Methyl, ::Val{:SMILES}; kwargs...) = "(C)"
getchemicalattr(::Methoxy, ::Val{:SMILES}; kwargs...) = "(OC)"
getchemicalattr(::Formyl, ::Val{:SMILES}; kwargs...) = "(C(=O))"
getchemicalattr(::Oformyl, ::Val{:SMILES}; kwargs...) = "(OC(=O))"
getchemicalattr(::Nformyl, ::Val{:SMILES}; kwargs...) = "(NC(=O))"
getchemicalattr(::Acetyl, ::Val{:SMILES}; kwargs...) = "(C(=O)C)"
getchemicalattr(::Oacetyl, ::Val{:SMILES}; kwargs...) = "(OC(=O)C)"
getchemicalattr(::Nacetyl, ::Val{:SMILES}; kwargs...) = "(NC(=O)C)"
getchemicalattr(::Bromo, ::Val{:SMILES}; kwargs...) = "(Br)"
getchemicalattr(::Chloro, ::Val{:SMILES}; kwargs...) = "(Cl)"
getchemicalattr(::Fluoro, ::Val{:SMILES}; kwargs...) = "(F)"
getchemicalattr(::Iodo, ::Val{:SMILES}; kwargs...) = "(I)"
getchemicalattr(::OxygenAtom, ::Val{:SMILES}; kwargs...) = "(O)"
getchemicalattr(::Hydroxy, ::Val{:SMILES}; kwargs...) = "(O)"
getchemicalattr(::OLinkage, ::Val{:SMILES}; kwargs...) = "O"
getchemicalattr(::CarboxylicAcidGroup, ::Val{:SMILES}; kwargs...) = "(C(=O)O)"
getchemicalattr(::CarboxylicLinkage, ::Val{:SMILES}; kwargs...) = "C(=O)"
getchemicalattr(::Oxo, ::Val{:SMILES}; kwargs...) = "(=O)"
getchemicalattr(::Ether, ::Val{:SMILES}; kwargs...) = "O"
getchemicalattr(::Hydroperoxyl, ::Val{:SMILES}; kwargs...) = "(OO)"
getchemicalattr(::Epoxy, ::Val{:SMILES}; kwargs...) = "(O1)X1"
getchemicalattr(::Peroxy, ::Val{:SMILES}; kwargs...) = "(OO1)X1"
getchemicalattr(::Amino, ::Val{:SMILES}; kwargs...) = "(N)"
getchemicalattr(::MethylAmino, ::Val{:SMILES}; kwargs...) = "(NC)"
getchemicalattr(::DimethylAmino, ::Val{:SMILES}; kwargs...) = "(NC(C))"
getchemicalattr(::TrimethylAmino, ::Val{:SMILES}; kwargs...) = "([N+]C(C)(C))"
getchemicalattr(::NLinkage, ::Val{:SMILES}; kwargs...) = "N"
getchemicalattr(::Sulfanyl, ::Val{:SMILES}; kwargs...) = "(S)"
getchemicalattr(::SLinkage, ::Val{:SMILES}; kwargs...) = "S"
getchemicalattr(::MethanethiolGroup, ::Val{:SMILES}; kwargs...) = "(SC)"
getchemicalattr(::Cyano, ::Val{:SMILES}; kwargs...) = "(C#N)"
getchemicalattr(::Phosphoryl, ::Val{:SMILES}; kwargs...) = "(P(=O)(O)O)"
getchemicalattr(::Substituent{Phosphoryl, Dehydroxy}, ::Val{:SMILES}; kwargs...) = "(P(=O)(O))"
getchemicalattr(::Diphosphoryl, ::Val{:SMILES}; kwargs...) = "(P(=O)(O)OP(=O)(O)O)"
getchemicalattr(::Triphosphoryl, ::Val{:SMILES}; kwargs...) = "(P(=O)(O)OP(=O)(O)OP(=O)(O)O)"
getchemicalattr(::Sulfo, ::Val{:SMILES}; kwargs...) = "(S(=O)(=O)O)"
# getchemicalattr(::Sulfate, ::Val{:SMILES}; kwargs...) = "(OS(=O)(=O)O)"
getchemicalattr(::Nitro, ::Val{:SMILES}; kwargs...) = "([N+](=O)[O-])"

isdissociated(::Hydrogen) = false
isdissociated(::Ethyl) = false
isdissociated(::Ethoxy) = false
isdissociated(::Methyl) = false
isdissociated(::Methoxy) = false
isdissociated(::Formyl) = false
isdissociated(::Oformyl) = false
isdissociated(::Nformyl) = false
isdissociated(::Acetyl) = false
isdissociated(::Oacetyl) = false
isdissociated(::Nacetyl) = false
isdissociated(::Bromo) = false
isdissociated(::Chloro) = false
isdissociated(::Fluoro) = false
isdissociated(::Iodo) = false
isdissociated(::OxygenAtom) = false
isdissociated(::Hydroxy) = false
isdissociated(::OLinkage) = false
isdissociated(::CarboxylicAcidGroup) = true
isdissociated(::CarboxylicLinkage) = false
isdissociated(::Oxo) = false
isdissociated(::Ether) = false
isdissociated(::Hydroperoxyl) = false
isdissociated(::Epoxy) = false
isdissociated(::Peroxy) = false
isdissociated(::Amino) = true
isdissociated(::NLinkage) = true
isdissociated(::Sulfanyl) = true
isdissociated(::SLinkage) = false
isdissociated(::MethanethiolGroup) = false
isdissociated(::Cyano) = false
isdissociated(::Phosphoryl) = true
isdissociated(::Substituent{Phosphoryl, Dehydroxy}) = true
isdissociated(::Diphosphoryl) = true
isdissociated(::Triphosphoryl) = true
isdissociated(::Sulfo) = true
isdissociated(::Nitro) = false
isdissociated(::Phenolhydroxy) = true
isdissociated(::Amide) = false
isdissociated(::Ester) = false
isdissociated(::Thioester) = false

dehydrogenposition(::BasicCompound) = nothing
dehydroxyposition(::BasicCompound) = missing
dehydroxyposition(::Methanol) = nothing
dehydroxyposition(::Ethanol) = nothing
dehydroxyposition(::FormicAcid) = nothing
dehydroxyposition(::AceticAcid) = nothing
dehydroxyposition(::PhosphoricAcid) = nothing
dehydroxyposition(::DiphosphoricAcid) = nothing
dehydroxyposition(::TriphosphoricAcid) = nothing
dehydroxyposition(::SulfuricAcid) = nothing
dehydroxyposition(::NitricAcid) = nothing

dehydrogengroup(::Dihydrogen; position = nothing) = Hydrogen()
dehydrogengroup(::Ethane; position = nothing) = Ethyl()
dehydrogengroup(::Ethanol; position = nothing) = Ethoxy()
dehydrogengroup(::Methane; position = nothing) = Methyl()
dehydrogengroup(::Methanol; position = nothing) = Methoxy()
dehydroxygroup(::FormicAcid; position = nothing) = Formyl()
dehydrogengroup(::FormicAcid; position = nothing) = Oformyl()
dehydrogengroup(::Formamide; position = nothing) = Nformyl()
dehydroxygroup(::AceticAcid; position = nothing) = Acetyl()
dehydrogengroup(::AceticAcid; position = nothing) = Oacetyl()
dehydrogengroup(::Acetamide; position = nothing) = Nacetyl()
dehydrogengroup(::HydrogenBromide; position = nothing) = Bromo()
dehydrogengroup(::HydrogenChloride; position = nothing) = Chloro()
dehydrogengroup(::HydrogenFluoride; position = nothing) = Fluoro()
dehydrogengroup(::HydrogenIodide; position = nothing) = Iodo()
dehydrogengroup(::HydrogenOxide; position = nothing) = Hydroxy()
dehydrogengroup(::HydrogenPeroxide; position = nothing) = Hydroperoxyl()
dehydrogengroup(::Ammonia; position = nothing) = Amino()
dehydrogengroup(::MethylAmine; position = nothing) = MethylAmino()
dehydrogengroup(::DimethylAmine; position = nothing) = DimethylAmino()
dehydrogengroup(::TrimethylAmine; position = nothing) = TrimethylAmino()
dehydrogengroup(::HydrogenSulfide; position = nothing) = Sulfanyl()
dehydrogengroup(::Methanethiol; position = nothing) = MethanethiolGroup()
dehydrogengroup(::HydrogenCyanide; position = nothing) = Cyano()
dehydroxygroup(::PhosphoricAcid; position = nothing) = Phosphoryl()
dehydroxygroup(::DiphosphoricAcid; position = nothing) = Diphosphoryl()
dehydroxygroup(::TriphosphoricAcid; position = nothing) = Triphosphoryl()
dehydroxygroup(::SulfuricAcid; position = nothing) = Sulfo()
# dehydrogengroup(::SulfuricAcid; position = nothing) = Sulfate()
dehydroxygroup(::NitricAcid; position = nothing) = Nitro()

nlinkage(::CarboxylicAcidGroup) = 3
nlinkage(::OxygenAtom) = 0
nlinkage(::Ether) = 1
nbridge(::Ether) = 1
nlinkage(::Epoxy) = 1
nbridge(::Epoxy) = 1
nlinkage(::Peroxy) = 1
nbridge(::Peroxy) = 1
nlinkage(::XLinkedFunctionalGroup{CarboxylicLinkage}) = 3

dehydroxyfunctionalgroup(::FormicAcid; position = nothing) = CarboxylicAcidGroup()
dehydroxyfunctionalgroup(::AceticAcid; position = nothing) = CarboxylicAcidGroup()

end