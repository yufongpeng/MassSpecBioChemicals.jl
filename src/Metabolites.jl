module Metabolites
using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
import ..MassSpecBioChemicals: parentchemical, leavinggroup, conjugation, dehydroxyposition, dehydrogenposition, leavinggroupelements
import MassSpecChemicals: getchemicalattr
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

getchemicalattr(::Ethanolamine, ::Val{:name}; kwargs...) = "Ethanolamine"
getchemicalattr(::Nmethylethanolamine, ::Val{:name}; kwargs...) = "N-Methylethanolamine"
getchemicalattr(::NNdimethylethanolamine, ::Val{:name}; kwargs...) = "Dimethylethanolamine"
getchemicalattr(::GABA, ::Val{:name}; kwargs...) = "GABA"
getchemicalattr(::Dopamine, ::Val{:name}; kwargs...) = "Dopamine"
getchemicalattr(::Taurine, ::Val{:name}; kwargs...) = "Taurine"
getchemicalattr(::Carnitine, ::Val{:name}; kwargs...) = "Carnitine"
getchemicalattr(::Choline, ::Val{:name}; kwargs...) = "Choline"
getchemicalattr(::CoA, ::Val{:name}; kwargs...) = "CoA"
getchemicalattr(::GlycolicAcid, ::Val{:name}; kwargs...) = "Glycolic acid"
getchemicalattr(::LacticAcid, ::Val{:name}; kwargs...) = "Lactic acid"
getchemicalattr(::PyruvicAcid, ::Val{:name}; kwargs...) = "Pyruvic acid"
getchemicalattr(::Glycerol, ::Val{:name}; kwargs...) = "Glycerol"

getchemicalattr(::Ethanolamine, ::Val{:formula}; kwargs...) = "HOCH2CH2NH2"
getchemicalattr(::Nmethylethanolamine, ::Val{:formula}; kwargs...) = "CH3NHCH2CH2OH"
getchemicalattr(::NNdimethylethanolamine, ::Val{:formula}; kwargs...) = "(CH3)2NCH2CH2OH"
getchemicalattr(::GABA, ::Val{:formula}; kwargs...) = "H2N(CH2)3COOH"
getchemicalattr(::Dopamine, ::Val{:formula}; kwargs...) = "C8H11NO2"
getchemicalattr(::Taurine, ::Val{:formula}; kwargs...) = "H3NCH2CH2SO3"
getchemicalattr(::Carnitine, ::Val{:formula}; kwargs...) = "(CH3)3NCH2CHOHCH2COO"
getchemicalattr(::Choline, ::Val{:formula}; kwargs...) = "(CH3)3NCH2CH2OH"
getchemicalattr(::CoA, ::Val{:formula}; kwargs...) = "C21H36N7O16P3S"
getchemicalattr(::GlycolicAcid, ::Val{:formula}; kwargs...) = "HOCH2COOH"
getchemicalattr(::LacticAcid, ::Val{:formula}; kwargs...) = "CH3CHOHCOOH"
getchemicalattr(::PyruvicAcid, ::Val{:formula}; kwargs...) = "CH3COCOOH"
getchemicalattr(::Glycerol, ::Val{:formula}; kwargs...) = "HOCH2CHOHCH2OH"

getchemicalattr(::Ethanolamine, ::Val{:SMILES}; kwargs...) = "NCCO"
getchemicalattr(::Nmethylethanolamine, ::Val{:SMILES}; kwargs...) = "CNCCO"
getchemicalattr(::NNdimethylethanolamine, ::Val{:SMILES}; kwargs...) = "CN(C)CCO"
getchemicalattr(::GABA, ::Val{:SMILES}; kwargs...) = "NCCCC(=O)[O-]"
getchemicalattr(::Dopamine, ::Val{:SMILES}; kwargs...) = "NCCc1cc(O)c(O)cc1"
getchemicalattr(::Taurine, ::Val{:SMILES}; kwargs...) = "NCCS(=O)(=O)O"
getchemicalattr(::Carnitine, ::Val{:SMILES}; kwargs...) = "C[N+](C)(C)CC(O)CC(=O)[O-]"
getchemicalattr(::Choline, ::Val{:SMILES}; kwargs...) = "C[N+](C)(C)CCO"
getchemicalattr(::CoA, ::Val{:SMILES}; kwargs...) = "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O"
getchemicalattr(::GlycolicAcid, ::Val{:SMILES}; kwargs...) = "OCC(=O)[O-]"
getchemicalattr(::LacticAcid, ::Val{:SMILES}; kwargs...) = "OC(C)C(=O)[O-]"
getchemicalattr(::PyruvicAcid, ::Val{:SMILES}; kwargs...) = "CC(=O)C(=O)[O-]"
getchemicalattr(::Glycerol, ::Val{:SMILES}; kwargs...) = "C(O)C(O)C(O)"

getchemicalattr(::Tauryl, ::Val{:abbreviation}; kwargs...) = "Tau"
getchemicalattr(::Glycoyl, ::Val{:abbreviation}; kwargs...) = "Gc"
getchemicalattr(::Lactyl, ::Val{:abbreviation}; kwargs...) = "Lt"
getchemicalattr(::Pyruvyl, ::Val{:abbreviation}; kwargs...) = "Py"
getchemicalattr(::Glyceryl, ::Val{:abbreviation}; kwargs...) = "Gr"
getchemicalattr(x::T, ::Val{:name}; kwargs...) where {T <: FunctionalGroup{<: Metabolite}} = string(T, " Group")
getchemicalattr(x::T, ::Val{:elements}; kwargs...) where {T <: FunctionalGroup{<: Metabolite}} = vcat(chemicalelements(parentchemical(x)), leavinggroupelements(leavinggroup(x)))

getchemicalattr(::Tauryl, ::Val{:SMILES}; kwargs...) = "(NCCS(=O)(=O)O)"
getchemicalattr(::Glycoyl, ::Val{:SMILES}; kwargs...) = "(C(=O)CO)"
getchemicalattr(::Lactyl, ::Val{:SMILES}; kwargs...) = "(C(=O)C(O)C)"
getchemicalattr(::Pyruvyl, ::Val{:SMILES}; kwargs...) = "(C(=O)C(=O)C)"
getchemicalattr(::Glyceryl, ::Val{:SMILES}; kwargs...) = "(OCC(O)C(O))"

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