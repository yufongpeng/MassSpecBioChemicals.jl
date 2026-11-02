module Metabolites
using ..MassSpecBioChemicals
using MassSpecChemicals: 
    AbstractChemical, AttributedChemical, ChemicalAttribute, make_attributedchemical
import MassSpecChemicals: 
    chemicalname, chemicalformula, chemicalelements, chemicalabbr, chemicalsmiles, 
    defaultchemicalkwargs
using ..MassSpecBioChemicals: GlyceraldehydeSystem
import ..MassSpecBioChemicals: 
        parentchemical, 
        leavinggroup, 
        conjugation, 
        dehydroxyposition, 
        dehydrogenposition, 
        leavinggroupelements, 
        dehydrogengroup, 
        dehydroxygroup, 
        dehydrogenfunctionalgroup, 
        dehydroxyfunctionalgroup, 
        chiralchemical, 
        chiralitize, 
        deisomerize, 
        decompose
        
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

abstract type Metabolite <: AttributedChemical end

struct Ethanolamine <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Ethanolamine(; kwargs...) = make_attributedchemical(Ethanolamine; kwargs...)

struct Nmethylethanolamine <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nmethylethanolamine(; kwargs...) = make_attributedchemical(Nmethylethanolamine; kwargs...)

struct NNdimethylethanolamine <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
NNdimethylethanolamine(; kwargs...) = make_attributedchemical(NNdimethylethanolamine; kwargs...)

struct GABA <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
GABA(; kwargs...) = make_attributedchemical(GABA; kwargs...)

struct Dopamine <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Dopamine(; kwargs...) = make_attributedchemical(Dopamine; kwargs...)

struct Taurine <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Taurine(; kwargs...) = make_attributedchemical(Taurine; kwargs...)
struct Tauryl <: FunctionalGroup{Taurine, Dehydrogen} end # T/Tau

struct Carnitine{L} <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Carnitine(; chirality = DLForm, kwargs...) = Core.apply_type(Carnitine, chirality)(; kwargs...)
Carnitine{L}(; kwargs...) where {L} = make_attributedchemical(Carnitine{L}; kwargs...)

struct Choline <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Choline(; kwargs...) = make_attributedchemical(Choline; kwargs...)

struct CoA <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
CoA(; kwargs...) = make_attributedchemical(CoA; kwargs...)

struct GlycolicAcid <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
GlycolicAcid(; kwargs...) = make_attributedchemical(GlycolicAcid; kwargs...)
struct Glycolyl <: FunctionalGroup{GlycolicAcid, Dehydroxy} end # Gc

struct LacticAcid{L} <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
LacticAcid(; chirality = DLForm, kwargs...) = Core.apply_type(LacticAcid, chirality)(; kwargs...)
LacticAcid{L}(; kwargs...) where {L} = make_attributedchemical(LacticAcid{L}; kwargs...)
struct Lactyl{L} <: FunctionalGroup{LacticAcid{L}, Dehydroxy} end # Lt

struct PyruvicAcid <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
PyruvicAcid(; kwargs...) = make_attributedchemical(PyruvicAcid; kwargs...)
struct Pyruvyl <: FunctionalGroup{PyruvicAcid, Dehydroxy} end # Py

struct Glycerol{R} <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glycerol(; chirality = RSChirality, kwargs...) = Core.apply_type(Glycerol, chirality)(; kwargs...)
Glycerol{R}(; kwargs...) where {R} = make_attributedchemical(Glycerol{R}; kwargs...)
struct Glyceryl{R} <: FunctionalGroup{Glycerol{R}, Dehydrogen} end # Gr

struct Serotonin{R} <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Serotonin(; chirality = RSChirality, kwargs...) = Core.apply_type(Serotonin, chirality)(; kwargs...)
Serotonin{R}(; kwargs...) where {R} = make_attributedchemical(Serotonin{R}; kwargs...)

struct Melatonin{R} <: Metabolite
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Melatonin(; chirality = RSChirality, kwargs...) = Core.apply_type(Melatonin, chirality)(; kwargs...)
Melatonin{R}(; kwargs...) where {R} = make_attributedchemical(Melatonin{R}; kwargs...)
# TCA cycle
# Seotonin, melanin
# Norepinephrine, epinephrine, 

include("msc.jl")
include("default.jl")
include("interface.jl")

end
