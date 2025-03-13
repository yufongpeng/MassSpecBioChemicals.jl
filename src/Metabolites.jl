module Metabolites
using Reexport
@reexport using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
import ..MassSpecBioChemicals: originalmolecule, leavinggroup, conjugation, dehydroxyposition, dehydrogenposition
import MassSpecChemicals: chemicalname, chemicalformula, chemicalabbr, repr_smiles
export Metabolite, 
       Ethanolamine, 
       Nmethylethanolamine, 
       NNdimethylethanolamine, 
       GABA, 
       Dopamine, 
       Taurine,
       Tauryl, # T/Tau
       Carnitine, 
       Choline, 
       CoA, 
       GlycolicAcid, 
       Glycoyl, # Gc
       LacticAcid, 
       Lactyl, # Lt
       PyruvicAcid, 
       Pyruvyl, # Py
       Glycerol,
       Glyceryl # Gr

abstract type Metabolite <: AbstractChemical end
struct Ethanolamine <: Metabolite end
struct Nmethylethanolamine <: Metabolite end
struct NNdimethylethanolamine <: Metabolite end
struct GABA <: Metabolite end
struct Dopamine <: Metabolite end
struct Taurine <: Metabolite end # T
struct Tauryl <: FunctionalGroup{Taurine, Dehydrogen} end # T/Tau
struct Carnitine <: Metabolite end
struct Choline <: Metabolite end
struct CoA <: Metabolite end
struct GlycolicAcid <: Metabolite end
struct Glycoyl <: FunctionalGroup{GlycolicAcid, Dehydroxy} end # Gc
struct LacticAcid <: Metabolite end
struct Lactyl <: FunctionalGroup{LacticAcid, Dehydroxy} end # Lt
struct PyruvicAcid <: Metabolite end
struct Pyruvyl <: FunctionalGroup{PyruvicAcid, Dehydroxy} end # Py
struct Glycerol <: Metabolite end
struct Glyceryl <: FunctionalGroup{Glycerol, Dehydrogen} end # Gr
conjugation(t::Taurine) = Tauryl()

chemicalname(::Ethanolamine) = "Ethanolamine"
chemicalname(::Nmethylethanolamine) = "N-Methylethanolamine"
chemicalname(::NNdimethylethanolamine) = "Dimethylethanolamine"
chemicalname(::GABA) = "GABA"
chemicalname(::Dopamine) = "Dopamine"
chemicalname(::Taurine) = "Taurine"
chemicalname(::Carnitine) = "Carnitine"
chemicalname(::Choline) = "Choline"
chemicalname(::CoA) = "CoA"
chemicalname(::GlycolicAcid) = "Glycolic acid"
chemicalname(::LacticAcid) = "Lactic acid"
chemicalname(::PyruvicAcid) = "Pyruvic acid"
chemicalname(::Glycerol) = "Glycerol"

chemicalformula(::Ethanolamine) = "HOCH2CH2NH2"
chemicalformula(::Nmethylethanolamine) = "CH3NHCH2CH2OH"
chemicalformula(::NNdimethylethanolamine) = "(CH3)2NCH2CH2OH"
chemicalformula(::GABA) = "H2N(CH2)3COOH"
chemicalformula(::Dopamine) = "C8H11NO2"
chemicalformula(::Taurine) = "H3NCH2CH2SO3"
chemicalformula(::Carnitine) = "(CH3)3NCH2CHOHCH2COO"
chemicalformula(::Choline) = "(CH3)3NCH2CH2OH"
chemicalformula(::CoA) = "C21H36N7O16P3S"
chemicalformula(::GlycolicAcid) = "HOCH2COOH"
chemicalformula(::LacticAcid) = "CH3CHOHCOOH"
chemicalformula(::PyruvicAcid) = "CH3COCOOH"
chemicalformula(::Glycerol) = "HOCH2CHOHCH2OH"

chemicalabbr(::Tauryl) = "Tau"
chemicalabbr(::Glycoyl) = "Gc"
chemicalabbr(::Lactyl) = "Lt"
chemicalabbr(::Pyruvyl) = "Py"
chemicalabbr(::Glyceryl) = "Gr"


repr_smiles(::Ethanolamine) = "NCCO"
repr_smiles(::Nmethylethanolamine) = "CNCCO"
repr_smiles(::NNdimethylethanolamine) = "CN(C)CCO"
repr_smiles(::GABA) = "NCCCC(=O)[O-]"
repr_smiles(::Dopamine) = "NCCc1cc(O)c(O)cc1"
repr_smiles(::Taurine) = "NCCS(=O)(=O)O"
repr_smiles(::Carnitine) = "C[N+](C)(C)CC(O)CC(=O)[O-]"
repr_smiles(::Choline) = "C[N+](C)(C)CCO"
repr_smiles(::CoA) = "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O"
repr_smiles(::GlycolicAcid) = "OCC(=O)[O-]"
repr_smiles(::LacticAcid) = "OC(C)C(=O)[O-]"
repr_smiles(::PyruvicAcid) = "CC(=O)C(=O)[O-]"
repr_smiles(::Glycerol) = "C(O)C(O)C(O)"

dehydrogenposition(::Ethanolamine) = nothing
dehydroxyposition(::Ethanolamine) = nothing
dehydrogenposition(::Nmethylethanolamine) = nothing
dehydroxyposition(::Nmethylethanolamine) = nothing
dehydrogenposition(::NNdimethylethanolamine) = missing
dehydroxyposition(::NNdimethylethanolamine) = nothing
dehydrogenposition(::GABA) = nothing
dehydroxyposition(::GABA) = nothing
dehydrogenposition(::Dopamine) = nothing
dehydroxyposition(::Dopamine) = missing
dehydrogenposition(::Taurine) = nothing
dehydroxyposition(::Taurine) = missing
dehydrogenposition(::Choline) = missing
dehydroxyposition(::Choline) = nothing
dehydrogenposition(::CoA) = nothing
dehydroxyposition(::CoA) = missing
dehydrogenposition(::GlycolicAcid) = nothing
dehydroxyposition(::GlycolicAcid) = nothing
dehydrogenposition(::LacticAcid) = nothing
dehydroxyposition(::LacticAcid) = nothing
dehydrogenposition(::PyruvicAcid) = missing
dehydroxyposition(::PyruvicAcid) = nothing
dehydrogenposition(::Glycerol) = 0x01
dehydroxyposition(::Glycerol) = 0x03
end