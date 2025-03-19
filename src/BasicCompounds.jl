module BasicCompounds
using Reexport
using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
import MassSpecChemicals: getchemicalattr
import ..MassSpecBioChemicals: dehydroxyposition, dehydrogenposition, leavinggroup, leavinggroupelements
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
       Alkoxy, # oxy 2
       HydrogenPeroxide, 
       Hydroperoxyl, # OOH
       Epoxy, # Ep 2
       Peroxy, # OO 2
       Ammonia, 
       Amino, # NH2/N
       NLinkage, # N...
       HydrogenSulfide, 
       Sulfanyl, # SH
       HydrogenCyanide, 
       Cyano, # CN
       PhosphoricAcid, # P
       Phosphate, # P/MP
       DiphosphoricAcid, 
       Diphosphate, # DP
       TriphosphoricAcid, 
       Triphosphate, # TP
       SulfuricAcid, # S
       Sulfate, # S
       SulfurousAcid, 
       Sulfite,
       NitricAcid, 
       Nitro # NO2


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
struct CarboxylicAcidGroup <: FunctionalGroup{FormicAcid, Demethyl} end # COOH
struct CarboxylicLinkage <: FunctionalGroup{CarboxylicAcidGroup, Dehydroxy} end # CO...
struct Oxo <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxo 2
struct Alkoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxy 2
struct HydrogenPeroxide <: BasicCompound end
struct Hydroperoxyl <: FunctionalGroup{HydrogenPeroxide, Dehydrogen} end # OOH
struct Epoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # Ep 2
struct Peroxy <: FunctionalGroup{HydrogenPeroxide, Didehydrogen} end # OO 2
struct Ammonia <: BasicCompound end
struct Amino <: FunctionalGroup{Ammonia, Dehydrogen} end # NH2/N
struct NLinkage <: FunctionalGroup{Amino, Dehydrogen} end # N...
struct HydrogenSulfide <: BasicCompound end
struct Sulfanyl <: FunctionalGroup{HydrogenSulfide, Dehydrogen} end # SH
struct HydrogenCyanide <: BasicCompound end
struct Cyano <: FunctionalGroup{HydrogenCyanide, Dehydrogen} end # CN
struct PhosphoricAcid <: BasicCompound end # P
struct Phosphate <: FunctionalGroup{PhosphoricAcid, Dehydroxy} end # P/MP
struct DiphosphoricAcid <: BasicCompound end
struct Diphosphate <: FunctionalGroup{DiphosphoricAcid, Dehydroxy} end # DP
struct TriphosphoricAcid <: BasicCompound end
struct Triphosphate <: FunctionalGroup{TriphosphoricAcid, Dehydroxy} end # TP
struct SulfuricAcid <: BasicCompound end # S
struct Sulfate <: FunctionalGroup{SulfuricAcid, Dehydroxy} end # S
struct SulfurousAcid <: BasicCompound end
struct Sulfite <: FunctionalGroup{SulfurousAcid, Dehydroxy} end
struct NitricAcid <: BasicCompound end
struct Nitro <: FunctionalGroup{NitricAcid, Dehydroxy} end # NO2

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
getchemicalattr(::HydrogenSulfide, ::Val{:name}; kwargs...) = "Hydrogen sulfide"
getchemicalattr(::HydrogenCyanide, ::Val{:name}; kwargs...) = "Hydrogen cyanide"
getchemicalattr(::PhosphoricAcid, ::Val{:name}; kwargs...) = "Phosphoric acid"
getchemicalattr(::DiphosphoricAcid, ::Val{:name}; kwargs...) = "Diphosphoric acid"
getchemicalattr(::TriphosphoricAcid, ::Val{:name}; kwargs...) = "Triphosphoric acid" 
getchemicalattr(::SulfuricAcid, ::Val{:name}; kwargs...) = "Sulfuric acid"
getchemicalattr(::SulfurousAcid, ::Val{:name}; kwargs...) = "Sulfurous acid"
getchemicalattr(::NitricAcid, ::Val{:name}; kwargs...) = "Nitric acid"

getchemicalattr(::Dihydrogen, ::Val{:formula}; kwargs...) = "H2"
getchemicalattr(::Ethane, ::Val{:formula}; kwargs...) = "CH3CH3"
getchemicalattr(::Ethanol, ::Val{:formula}; kwargs...) = "CH3CH2OH"
getchemicalattr(::Methane, ::Val{:formula}; kwargs...) = "CH4" 
getchemicalattr(::Methanol, ::Val{:formula}; kwargs...) = "CH3OH"
getchemicalattr(::FormicAcid, ::Val{:formula}; kwargs...) = "HCOOH" 
getchemicalattr(::Formamide, ::Val{:formula}; kwargs...) = "HCONH2"
getchemicalattr(::AceticAcid, ::Val{:formula}; kwargs...) = "CH3COOH"
getchemicalattr(::Acetamide, ::Val{:formula}; kwargs...) = "CH3CONH2"
getchemicalattr(::HydrogenBromide, ::Val{:formula}; kwargs...) = "HBr"
getchemicalattr(::HydrogenChloride, ::Val{:formula}; kwargs...) = "HCl"
getchemicalattr(::HydrogenFluoride, ::Val{:formula}; kwargs...) = "HF"
getchemicalattr(::HydrogenIodide, ::Val{:formula}; kwargs...) = "HI"
getchemicalattr(::HydrogenOxide, ::Val{:formula}; kwargs...) = "H2O"
getchemicalattr(::HydrogenPeroxide, ::Val{:formula}; kwargs...) = "HOOH"
getchemicalattr(::Ammonia, ::Val{:formula}; kwargs...) = "NH3"
getchemicalattr(::HydrogenSulfide, ::Val{:formula}; kwargs...) = "H2S"
getchemicalattr(::HydrogenCyanide, ::Val{:formula}; kwargs...) = "HCN"
getchemicalattr(::PhosphoricAcid, ::Val{:formula}; kwargs...) = "H3PO4"
getchemicalattr(::DiphosphoricAcid, ::Val{:formula}; kwargs...) = "H4P2O7"
getchemicalattr(::TriphosphoricAcid, ::Val{:formula}; kwargs...) = "H5P3O10" 
getchemicalattr(::SulfuricAcid, ::Val{:formula}; kwargs...) = "H2SO4"
getchemicalattr(::SulfurousAcid, ::Val{:formula}; kwargs...) = "H2SO3"
getchemicalattr(::NitricAcid, ::Val{:formula}; kwargs...) = "HNO3"

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
getchemicalattr(::HydrogenSulfide, ::Val{:SMILES}; kwargs...) = "S"
getchemicalattr(::HydrogenCyanide, ::Val{:SMILES}; kwargs...) = "C#N"
getchemicalattr(::PhosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)O"
getchemicalattr(::DiphosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)OP(=O)(O)O"
getchemicalattr(::TriphosphoricAcid, ::Val{:SMILES}; kwargs...) = "OP(=O)(O)OP(=O)(O)OP(=O)(O)O" 
getchemicalattr(::SulfuricAcid, ::Val{:SMILES}; kwargs...) = "OS(=O)(=O)O"
getchemicalattr(::SulfurousAcid, ::Val{:SMILES}; kwargs...) = "OS(=O)O"
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
getchemicalattr(::Alkoxy, ::Val{:abbreviation}; kwargs...) = "oxy"
getchemicalattr(::Hydroperoxyl, ::Val{:abbreviation}; kwargs...) = "OOH"
getchemicalattr(::Epoxy, ::Val{:abbreviation}; kwargs...) = "Ep"
getchemicalattr(::Peroxy, ::Val{:abbreviation}; kwargs...) = "OO"
getchemicalattr(::Amino, ::Val{:abbreviation}; kwargs...) = "NH2"
getchemicalattr(::NLinkage, ::Val{:abbreviation}; kwargs...) = "N"
getchemicalattr(::Sulfanyl, ::Val{:abbreviation}; kwargs...) = "SH"
getchemicalattr(::Cyano, ::Val{:abbreviation}; kwargs...) = "CN"
getchemicalattr(::Phosphate, ::Val{:abbreviation}; kwargs...) = "P"
getchemicalattr(::Diphosphate, ::Val{:abbreviation}; kwargs...) = "DP"
getchemicalattr(::Triphosphate, ::Val{:abbreviation}; kwargs...) = "TP"
getchemicalattr(::Sulfate, ::Val{:abbreviation}; kwargs...) = "S"
getchemicalattr(::Sulfite, ::Val{:abbreviation}; kwargs...) = "HSO3"
getchemicalattr(::Nitro, ::Val{:abbreviation}; kwargs...) = "NO2"
getchemicalattr(x::T, ::Val{:name}; kwargs...) where {T <: FunctionalGroup{<: BasicCompound}} = string(T, " Group")
getchemicalattr(x::T, ::Val{:elements}; kwargs...) where {T <: FunctionalGroup{<: BasicCompound}} = vcat(chemicalelements(parentchemical(x)), leavinggroupelements(leavinggroup(x)))

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
getchemicalattr(::Alkoxy, ::Val{:SMILES}; kwargs...) = "O"
getchemicalattr(::Hydroperoxyl, ::Val{:SMILES}; kwargs...) = "(OO)"
getchemicalattr(::Epoxy, ::Val{:SMILES}; kwargs...) = "(O1)X1"
getchemicalattr(::Peroxy, ::Val{:SMILES}; kwargs...) = "OO"
getchemicalattr(::Amino, ::Val{:SMILES}; kwargs...) = "(N)"
getchemicalattr(::NLinkage, ::Val{:SMILES}; kwargs...) = "N"
getchemicalattr(::Sulfanyl, ::Val{:SMILES}; kwargs...) = "(S)"
getchemicalattr(::Cyano, ::Val{:SMILES}; kwargs...) = "(C#N)"
getchemicalattr(::Phosphate, ::Val{:SMILES}; kwargs...) = "(OP(=O)(O)O)"
getchemicalattr(::Diphosphate, ::Val{:SMILES}; kwargs...) = "(OP(=O)(O)OP(=O)(O)O)"
getchemicalattr(::Triphosphate, ::Val{:SMILES}; kwargs...) = "(OP(=O)(O)OP(=O)(O)OP(=O)(O)O)"
getchemicalattr(::Sulfate, ::Val{:SMILES}; kwargs...) = "(OS(=O)(=O)O)"
getchemicalattr(::Sulfite, ::Val{:SMILES}; kwargs...) = "(OS(=O)O)"
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
isdissociated(::Alkoxy) = false
isdissociated(::Hydroperoxyl) = false
isdissociated(::Epoxy) = false
isdissociated(::Peroxy) = false
isdissociated(::Amino) = true
isdissociated(::NLinkage) = true
isdissociated(::Sulfanyl) = true
isdissociated(::Cyano) = false
isdissociated(::Phosphate) = true
isdissociated(::Diphosphate) = true
isdissociated(::Triphosphate) = true
isdissociated(::Sulfate) = true
isdissociated(::Sulfite) = true
isdissociated(::Nitro) = false

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
dehydroxyposition(::SulfurousAcid) = nothing
dehydroxyposition(::NitricAcid) = nothing

end