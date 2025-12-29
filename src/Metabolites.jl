module Metabolites
using ..MassSpecBioChemicals
using MassSpecChemicals: AbstractChemical
using ..MassSpecBioChemicals: GlyceraldehydeSystem
import ..MassSpecBioChemicals: parentchemical, leavinggroup, conjugation, dehydroxyposition, dehydrogenposition, leavinggroupelements, dehydrogengroup, dehydroxygroup, dehydrogenfunctionalgroup, dehydroxyfunctionalgroup, chiralchemical, deisomerize, decompose
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
       Glycolyl, # Gc
       LacticAcid, 
       Lactyl, # Lt
       PyruvicAcid, 
       Pyruvyl, # Py
       Glycerol,
       Glyceryl, # Gr
       Serotonin, 
       Melatonin

abstract type Metabolite <: AbstractChemical end
struct Ethanolamine <: Metabolite end
struct Nmethylethanolamine <: Metabolite end
struct NNdimethylethanolamine <: Metabolite end
struct GABA <: Metabolite end
struct Dopamine <: Metabolite end
struct Taurine <: Metabolite end # T
struct Tauryl <: FunctionalGroup{Taurine, Dehydrogen} end # T/Tau
struct Carnitine{L} <: Metabolite end
struct Choline <: Metabolite end
struct CoA <: Metabolite end
struct GlycolicAcid <: Metabolite end
struct Glycolyl <: FunctionalGroup{GlycolicAcid, Dehydroxy} end # Gc
struct LacticAcid{L} <: Metabolite end
struct Lactyl{L} <: FunctionalGroup{LacticAcid{L}, Dehydroxy} end # Lt
struct PyruvicAcid <: Metabolite end
struct Pyruvyl <: FunctionalGroup{PyruvicAcid, Dehydroxy} end # Py
struct Glycerol{R} <: Metabolite end
struct Glyceryl{R} <: FunctionalGroup{Glycerol{R}, Dehydrogen} end # Gr
struct Serotonin{R} <: Metabolite end # 5HT 
struct Melatonin{R} <: Metabolite end # MT 
# TCA cycle
# Seotonin, melanin
# Norepinephrine, epinephrine, 

chiralchemical(::Type{Carnitine}, ::Type{I}) where {I <: GlyceraldehydeSystem} = Carnitine{I}()
chiralchemical(::Type{Carnitine}, ::Type{RSChirality}) = Carnitine{DLForm}()
chiralchemical(::Type{Carnitine}, ::Type{SChirality}) = Carnitine{DForm}()
chiralchemical(::Type{Carnitine}, ::Type{RChirality}) = Carnitine{LForm}()
chiralchemical(::Type{LacticAcid}, ::Type{I}) where {I <: GlyceraldehydeSystem} = LacticAcid{I}()
chiralchemical(::Type{LacticAcid}, ::Type{RSChirality}) = LacticAcid{DLForm}()
chiralchemical(::Type{LacticAcid}, ::Type{RChirality}) = LacticAcid{DForm}()
chiralchemical(::Type{LacticAcid}, ::Type{SChirality}) = LacticAcid{LForm}()
chiralchemical(::Type{Lactyl}, ::Type{I}) where {I <: GlyceraldehydeSystem} = Lactyl{I}()
chiralchemical(::Type{Lactyl}, ::Type{RSChirality}) = Lactyl{DLForm}()
chiralchemical(::Type{Lactyl}, ::Type{RChirality}) = Lactyl{DForm}()
chiralchemical(::Type{Lactyl}, ::Type{SChirality}) = Lactyl{LForm}()
chiralchemical(::Type{Glycerol}, ::Type{I}) where {I <: RSChirality} = Glycerol{I}()
chiralchemical(::Type{Glyceryl}, ::Type{I}) where {I <: RSChirality} = Glyceryl{I}()
chiralchemical(::Type{Serotonin}, ::Type{I}) where {I <: RSChirality} = Serotonin{I}()
chiralchemical(::Type{Melatonin}, ::Type{I}) where {I <: RSChirality} = Melatonin{I}()

dehydrogengroup(::Taurine; position = nothing) = Tauryl()
dehydroxygroup(::GlycolicAcid; position = nothing) = Glycolyl()
dehydroxygroup(::LacticAcid{T}; position = nothing) where T = Lactyl{T}()
dehydroxygroup(::PyruvicAcid; position = nothing) = Pyruvyl()
dehydrogengroup(::Glycerol{T}; position = nothing) where T = Glyceryl{T}()
conjugation(::Taurine) = Tauryl()

getchemicalattr(::Ethanolamine, ::Val{:name}; kwargs...) = "Ethanolamine"
getchemicalattr(::Nmethylethanolamine, ::Val{:name}; kwargs...) = "N-Methylethanolamine"
getchemicalattr(::NNdimethylethanolamine, ::Val{:name}; kwargs...) = "Dimethylethanolamine"
getchemicalattr(::GABA, ::Val{:name}; kwargs...) = "GABA"
getchemicalattr(::Dopamine, ::Val{:name}; kwargs...) = "Dopamine"
getchemicalattr(::Taurine, ::Val{:name}; kwargs...) = "Taurine"
getchemicalattr(::Carnitine, ::Val{:name}; kwargs...) = "Carnitine"
getchemicalattr(::Carnitine{LForm}, ::Val{:name}; kwargs...) = "L-Carnitine"
getchemicalattr(::Carnitine{DForm}, ::Val{:name}; kwargs...) = "D-Carnitine"
getchemicalattr(::Choline, ::Val{:name}; kwargs...) = "Choline"
getchemicalattr(::CoA, ::Val{:name}; kwargs...) = "CoA"
getchemicalattr(::GlycolicAcid, ::Val{:name}; kwargs...) = "Glycolic acid"
getchemicalattr(::LacticAcid, ::Val{:name}; kwargs...) = "Lactic acid"
getchemicalattr(::LacticAcid{LForm}, ::Val{:name}; kwargs...) = "L-Lactic acid"
getchemicalattr(::LacticAcid{DForm}, ::Val{:name}; kwargs...) = "D-Lactic acid"
getchemicalattr(::PyruvicAcid, ::Val{:name}; kwargs...) = "Pyruvic acid"
getchemicalattr(::Glycerol, ::Val{:name}; kwargs...) = "Glycerol"
getchemicalattr(::Glycerol{RChirality}, ::Val{:name}; kwargs...) = "(R)-Glycerol"
getchemicalattr(::Glycerol{SChirality}, ::Val{:name}; kwargs...) = "(S)-Glycerol"
getchemicalattr(::Serotonin, ::Val{:name}; kwargs...) = "Serotonin"
getchemicalattr(::Melatonin, ::Val{:name}; kwargs...) = "Melatonin"


getchemicalattr(::Ethanolamine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "O" => 1, "H" => 1]
getchemicalattr(::Nmethylethanolamine, ::Val{:elements}; kwargs...) = ["H" => 1, "C" => 1, "H" => 3, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "O" => 1, "H" => 1]
getchemicalattr(::NNdimethylethanolamine, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "O" => 1, "H" => 1]
getchemicalattr(::GABA, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Dopamine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "C" => 6, "H" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Taurine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "S" => 1, "O" => 2, "O" => 1, "H" => 1]
getchemicalattr(::Carnitine, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 3, "H" => 5, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1]
getchemicalattr(::Choline, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 2, "H" => 3, "O" => 1, "H" => 1]
getchemicalattr(::CoA, ::Val{:elements}; kwargs...) = ["H" => 1, "S" => 1, "C" => 2, "H" => 4, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 2, "H" => 4, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 5, "H" => 9, "O" => 1, "H" => 1, "O" => 1, "P" => 2, "O" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 5, "H" => 8, "O" => 3, "O" => 1, "P" => 1, "O" => 1, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 5, "H" => 2, "N" => 4, "N" => 1, "H" => 2]
getchemicalattr(::GlycolicAcid, ::Val{:elements}; kwargs...) = ["H" => 1, "O" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::LacticAcid, ::Val{:elements}; kwargs...) = ["H" => 1, "O" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::PyruvicAcid, ::Val{:elements}; kwargs...) = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Glycerol, ::Val{:elements}; kwargs...) = ["H" => 1, "O" => 1, "C" => 3, "H" => 5, "O" => 1, "H" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Serotonin, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 8, "H" => 5, "N" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Melatonin, ::Val{:elements}; kwargs...) = ["H" => 1, "N" => 1, "C" => 1, "O" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 8, "H" => 5, "N" => 1, "O" => 1, "C" => 1, "H" => 3]
getchemicalattr(m::Metabolite, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)

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
function getchemicalattr(::Glycerol, ::Val{:SMILES}; chain = (:O1, :C3))
    if chain == (:O1, :C3)
        "OCC(O)C(O)"
    elseif chain == (:C3, :O1)
        "C(O)C(O)CO"
    elseif chain == (:C3, :C1)
        "C(O)C(O)C(O)"
    elseif chain == (:C1, :C3)
        "C(O)C(O)C(O)"
    elseif chain == (:O3, :C1)
        "OCC(O)C(O)"
    elseif chain == (:C1, :O3)
        "C(O)C(O)CO"
    elseif chain == (:O3, :O1)
        "OCC(O)CO"
    elseif chain == (:O1, :O3)
        "OCC(O)CO"
    end
end
getchemicalattr(::Serotonin, ::Val{:SMILES}; kwargs...) = "NCCC1=CNC2=C1C=C(C=C2)O"
getchemicalattr(::Melatonin, ::Val{:SMILES}; kwargs...) = "N(C(=O)C)CCC1=CNC2=C1C=C(C=C2)OC"


getchemicalattr(::Tauryl, ::Val{:abbreviation}; kwargs...) = "Tau"
getchemicalattr(::Glycolyl, ::Val{:abbreviation}; kwargs...) = "Gc"
getchemicalattr(::Lactyl, ::Val{:abbreviation}; kwargs...) = "Lt"
getchemicalattr(::Pyruvyl, ::Val{:abbreviation}; kwargs...) = "Py"
getchemicalattr(::Glyceryl, ::Val{:abbreviation}; kwargs...) = "Gr"
getchemicalattr(::Tauryl, ::Val{:name}; kwargs...) = "Tauryl"
getchemicalattr(::Glycolyl, ::Val{:name}; kwargs...) = "Glycolyl"
getchemicalattr(::Lactyl, ::Val{:name}; kwargs...) = "Lactyl"
getchemicalattr(::Lactyl{DForm}, ::Val{:name}; kwargs...) = "D-Lactyl"
getchemicalattr(::Lactyl{LForm}, ::Val{:name}; kwargs...) = "L-Lactyl"
getchemicalattr(::Pyruvyl, ::Val{:name}; kwargs...) = "Pyruvyl"
getchemicalattr(::Glyceryl, ::Val{:name}; kwargs...) = "Glyceryl"
getchemicalattr(::Glyceryl{RChirality}, ::Val{:name}; kwargs...) = "(R)-Glyceryl"
getchemicalattr(::Glyceryl{SChirality}, ::Val{:name}; kwargs...) = "(S)-Glyceryl"
getchemicalattr(x::T, ::Val{:name}; kwargs...) where {T <: FunctionalGroup{<: Metabolite}} = string(T, " Group")
function getchemicalattr(x::T, ::Val{:elements}; kwargs...) where {T <: FunctionalGroup{<: Metabolite}}
    es = chemicalelements(parentchemical(x))
	ls = leavinggroupelements(leavinggroup(x))
	for (e, n) in ls 
		if n > 0
			i = findfirst(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = findfirst(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	es
end 
getchemicalattr(x::T, ::Val{:formula}; unique = false, kwargs...) where {T <: FunctionalGroup{<: Metabolite}} = chemicalformula(chemicalelements(x); unique, kwargs...)

getchemicalattr(::Tauryl, ::Val{:SMILES}; kwargs...) = "(NCCS(=O)(=O)O)"
getchemicalattr(::Glycolyl, ::Val{:SMILES}; kwargs...) = "(C(=O)CO)"
getchemicalattr(::Lactyl, ::Val{:SMILES}; kwargs...) = "(C(=O)C(O)C)"
getchemicalattr(::Pyruvyl, ::Val{:SMILES}; kwargs...) = "(C(=O)C(=O)C)"
getchemicalattr(::Glyceryl, ::Val{:SMILES}; kwargs...) = "(OCC(O)C(O))"

dehydrogenposition(::Ethanolamine) = nothing
dehydroxyposition(::Ethanolamine) = nothing
dehydrogenfunctionalgroup(::Ethanolamine; position = nothing) = Amino()
dehydroxyfunctionalgroup(::Ethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::Nmethylethanolamine) = nothing
dehydroxyposition(::Nmethylethanolamine) = nothing
dehydrogenfunctionalgroup(::Nmethylethanolamine; position = nothing) = MethylAmino()
dehydroxyfunctionalgroup(::Nmethylethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::NNdimethylethanolamine) = missing
dehydroxyposition(::NNdimethylethanolamine) = nothing
dehydroxyfunctionalgroup(::NNdimethylethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::GABA) = nothing
dehydroxyposition(::GABA) = nothing
dehydrogenfunctionalgroup(::GABA; position = nothing) = Amino()
dehydroxyfunctionalgroup(::GABA; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::Dopamine) = nothing
dehydroxyposition(::Dopamine) = missing
dehydrogenfunctionalgroup(::Dopamine; position = nothing) = Amino()
dehydrogenposition(::Taurine) = nothing
dehydrogenfunctionalgroup(::Taurine; position = nothing) = Sulfo()
dehydroxyposition(::Taurine) = missing
dehydrogenposition(::Choline) = missing
dehydroxyposition(::Choline) = nothing
dehydroxyfunctionalgroup(::Choline; position = nothing) = Hydroxy()
dehydrogenposition(::CoA) = nothing
dehydroxyposition(::CoA) = missing
dehydrogenfunctionalgroup(::CoA; position = nothing) = Sulfo()
dehydrogenposition(::GlycolicAcid) = nothing
dehydroxyposition(::GlycolicAcid) = nothing
dehydrogenfunctionalgroup(::GlycolicAcid; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::GlycolicAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::LacticAcid) = nothing
dehydroxyposition(::LacticAcid) = nothing
dehydrogenfunctionalgroup(::LacticAcid; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::LacticAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::PyruvicAcid) = missing
dehydroxyposition(::PyruvicAcid) = nothing
dehydroxyfunctionalgroup(::PyruvicAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::Glycerol) = 0x01
dehydroxyposition(::Glycerol) = 0x03
dehydrogenfunctionalgroup(::LacticAcid; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::LacticAcid; position = nothing) = Hydroxy()

deisomerize(::Carnitine) = Carnitine{DLForm}()
deisomerize(::LacticAcid) = LacticAcid{DLForm}()
deisomerize(::Glycerol) = Glycerol{RSChirality}()

decompose(::Ethanolamine) = Dict(Hydroxy() => 0x01, Amino() => 0x01)
decompose(::Nmethylethanolamine) = Dict(Hydroxy() => 0x01, MethylAmino() => 0x01) # N-alky
decompose(::NNdimethylethanolamine) = Dict(Hydroxy() => 0x01, DimethylAmino() => 0x01)
decompose(::GABA) = Dict(CarboxylicAcidGroup() => 0x01, Amino() => 0x01)
decompose(::Dopamine) = Dict(Phenolhydroxy() => 0x02, Amino() => 0x01)
decompose(::Taurine) = Dict(Sulfo() => 0x01, Amino() => 0x01)
decompose(::Carnitine) = Dict(TrimethylAmino() => 0x01, Hydroxy() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Choline) = Dict(TrimethylAmino() => 0x01, Hydroxy() => 0x01)
decompose(::CoA) = Dict(CoA() => 0x01)
decompose(::GlycolicAcid) = Dict(CarboxylicAcidGroup() => 0x01, Hydroxy() => 0x01)
decompose(::LacticAcid) = Dict(CarboxylicAcidGroup() => 0x01, Hydroxy() => 0x01)
decompose(::PyruvicAcid) = Dict(CarboxylicAcidGroup() => 0x01, Oxo() => 0x01)
decompose(::Glycerol) = Dict(Hydroxy() => 0x03)
decompose(::Serotonin) = Dict(Phenolhydroxy() => 0x01, Amino() => 0x01, MethylAmino() => 0x01)
decompose(::Melatonin) = Dict(Methoxy() => 0x01, Nacetyl() => 0x01, MethylAmino() => 0x01)

end