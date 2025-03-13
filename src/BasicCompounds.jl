module BasicCompounds
using Reexport
@reexport using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
import MassSpecChemicals: chemicalname, chemicalformula, chemicalabbr, repr_smiles
import ..MassSpecBioChemicals: dehydroxyposition, dehydrogenposition
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

chemicalname(::Dihydrogen) = "H2"
chemicalname(::Ethane) = "Ethane"
chemicalname(::Ethanol) = "EtOH"
chemicalname(::Methane) = "Methane" 
chemicalname(::Methanol) = "MeOH"
chemicalname(::FormicAcid) = "Formic acid" 
chemicalname(::Formamide) = "Formamide"
chemicalname(::AceticAcid) = "Acetic acid"
chemicalname(::Acetamide) = "Acetamide"
chemicalname(::HydrogenBromide) = "Hydrogen bromide"
chemicalname(::HydrogenChloride) = "Hydrogen chloride"
chemicalname(::HydrogenFluoride) = "Hydrogen fluoride"
chemicalname(::HydrogenIodide) = "Hydrogen iodide"
chemicalname(::HydrogenOxide) = "Hydrogen oxide"
chemicalname(::HydrogenPeroxide) = "Hydrogen peroxide"
chemicalname(::Ammonia) = "Ammonia"
chemicalname(::HydrogenSulfide) = "Hydrogen sulfide"
chemicalname(::HydrogenCyanide) = "Hydrogen cyanide"
chemicalname(::PhosphoricAcid) = "Phosphoric acid"
chemicalname(::DiphosphoricAcid) = "Diphosphoric acid"
chemicalname(::TriphosphoricAcid) = "Triphosphoric acid" 
chemicalname(::SulfuricAcid) = "Sulfuric acid"
chemicalname(::SulfurousAcid) = "Sulfurous acid"
chemicalname(::NitricAcid) = "Nitric acid"

chemicalformula(::Dihydrogen) = "H2"
chemicalformula(::Ethane) = "CH3CH3"
chemicalformula(::Ethanol) = "CH3CH2OH"
chemicalformula(::Methane) = "CH4" 
chemicalformula(::Methanol) = "CH3OH"
chemicalformula(::FormicAcid) = "HCOOH" 
chemicalformula(::Formamide) = "HCONH2"
chemicalformula(::AceticAcid) = "CH3COOH"
chemicalformula(::Acetamide) = "CH3CONH2"
chemicalformula(::HydrogenBromide) = "HBr"
chemicalformula(::HydrogenChloride) = "HCl"
chemicalformula(::HydrogenFluoride) = "HF"
chemicalformula(::HydrogenIodide) = "HI"
chemicalformula(::HydrogenOxide) = "H2O"
chemicalformula(::HydrogenPeroxide) = "HOOH"
chemicalformula(::Ammonia) = "NH3"
chemicalformula(::HydrogenSulfide) = "H2S"
chemicalformula(::HydrogenCyanide) = "HCN"
chemicalformula(::PhosphoricAcid) = "H3PO4"
chemicalformula(::DiphosphoricAcid) = "H4P2O7"
chemicalformula(::TriphosphoricAcid) = "H5P3O10" 
chemicalformula(::SulfuricAcid) = "H2SO4"
chemicalformula(::SulfurousAcid) = "H2SO3"
chemicalformula(::NitricAcid) = "HNO3"

chemicalabbr(::Hydrogen) = "H"
chemicalabbr(::Ethyl) = "Et"
chemicalabbr(::Ethoxy) = "OEt"
chemicalabbr(::Methyl) = "Me"
chemicalabbr(::Methoxy) = "OMe"
chemicalabbr(::Formyl) = "Fo"
chemicalabbr(::Oformyl) = "OFo"
chemicalabbr(::Nformyl) = "NFo"
chemicalabbr(::Acetyl) = "Ac"
chemicalabbr(::Oacetyl) = "OAc"
chemicalabbr(::Nacetyl) = "NAc"
chemicalabbr(::Bromo) = "Br"
chemicalabbr(::Chloro) = "Cl"
chemicalabbr(::Fluoro) = "F"
chemicalabbr(::Iodo) = "I"
chemicalabbr(::OxygenAtom) = "O"
chemicalabbr(::Hydroxy) = "OH"
chemicalabbr(::OLinkage) = "O"
chemicalabbr(::CarboxylicAcidGroup) = "COOH"
chemicalabbr(::CarboxylicLinkage) = "CO"
chemicalabbr(::Oxo) = "oxo"
chemicalabbr(::Alkoxy) = "oxy"
chemicalabbr(::Hydroperoxyl) = "OOH"
chemicalabbr(::Epoxy) = "Ep"
chemicalabbr(::Peroxy) = "OO"
chemicalabbr(::Amino) = "NH2"
chemicalabbr(::NLinkage) = "N"
chemicalabbr(::Sulfanyl) = "SH"
chemicalabbr(::Cyano) = "CN"
chemicalabbr(::Phosphate) = "P"
chemicalabbr(::Diphosphate) = "DP"
chemicalabbr(::Triphosphate) = "TP"
chemicalabbr(::Sulfate) = "S"
chemicalabbr(::Sulfite) = "HSO3"
chemicalabbr(::Nitro) = "NO2"

repr_smiles(::Hydrogen) = "(H)"
repr_smiles(::Ethyl) = "(CC)"
repr_smiles(::Ethoxy) = "(OCC)"
repr_smiles(::Methyl) = "(C)"
repr_smiles(::Methoxy) = "(OC)"
repr_smiles(::Formyl) = "(C(=O))"
repr_smiles(::Oformyl) = "(OC(=O))"
repr_smiles(::Nformyl) = "(NC(=O))"
repr_smiles(::Acetyl) = "(C(=O)C)"
repr_smiles(::Oacetyl) = "(OC(=O)C)"
repr_smiles(::Nacetyl) = "(NC(=O)C)"
repr_smiles(::Bromo) = "(Br)"
repr_smiles(::Chloro) = "(Cl)"
repr_smiles(::Fluoro) = "(F)"
repr_smiles(::Iodo) = "(I)"
repr_smiles(::OxygenAtom) = "(O)"
repr_smiles(::Hydroxy) = "(O)"
repr_smiles(::OLinkage) = "O"
repr_smiles(::CarboxylicAcidGroup) = "(C(=O)O)"
repr_smiles(::CarboxylicLinkage) = "C(=O)"
repr_smiles(::Oxo) = "(=O)"
repr_smiles(::Alkoxy) = "O"
repr_smiles(::Hydroperoxyl) = "(OO)"
repr_smiles(::Epoxy) = "(O1)X1"
repr_smiles(::Peroxy) = "OO"
repr_smiles(::Amino) = "(N)"
repr_smiles(::NLinkage) = "N"
repr_smiles(::Sulfanyl) = "(S)"
repr_smiles(::Cyano) = "(C#N)"
repr_smiles(::Phosphate) = "(OP(=O)(O)O)"
repr_smiles(::Diphosphate) = "(OP(=O)(O)OP(=O)(O)O)"
repr_smiles(::Triphosphate) = "(OP(=O)(O)OP(=O)(O)OP(=O)(O)O)"
repr_smiles(::Sulfate) = "(OS(=O)(=O)O)"
repr_smiles(::Sulfite) = "(OS(=O)O)"
repr_smiles(::Nitro) = "([N+](=O)[O-])"

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