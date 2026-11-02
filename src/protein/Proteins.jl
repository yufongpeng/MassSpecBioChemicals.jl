module Proteins
using Reexport
import Base: length
using ..MassSpecBioChemicals
using ..MassSpecBioChemicals: lk, GlyceraldehydeSystem, RSSystem
using MassSpecChemicals: 
    AbstractChemical, AttributedChemical, ChemicalAttribute, AbstractChemicalParser,
    make_attributedchemical, buildattribute
import MassSpecChemicals: 
    chemicalname, chemicalformula, chemicalelements, chemicalabbr, chemicalsmiles, 
    defaultchemicalkwargs, parse_chemical
import ..MassSpecBioChemicals: 
        parentchemical, 
        leavinggroup, 
        conjugation, 
        chainedchemical, 
        getchaincomponent, 
        getchainlinkage, 
        getchainconfig, 
        composition, 
        ischainedchemical, 
        requirelinkage, 
        requireconfig, 
        dehydroxyposition, 
        dehydrogenposition, 
        dehydrogengroup, 
        dehydroxygroup, 
        dehydrogenfunctionalgroup, 
        dehydroxyfunctionalgroup, 
        chiralchemical, 
        chiralitize, 
        losswaterelements, 
        deisomerize, 
        decompose
        
export AbstractPeptide,
        αAminoAcid, 
       Alanine,  
       Alanyl,  
       Arginine,  
       Arginyl,  
       Aspargine,  
       Asparginyl,  
       Aspartate,  
       Aspartyl,  
       Cysteine,  
       Cysteinyl,  
       Glutamine,  
       Glutaminyl,  
       Glutamate,  
       Glutamyl,  
       Glycine,  
       Glycyl,  
       Histidine,  
       Histidinyl,  
       Isoleucine,  
       Isoleucyl,  
       Leucine,  
       Leucyl,  
       Lysine,  
       Lysyl,  
       Methionine,  
       Methionyl,  
       Phenylalanine,  
       Phenylalanyl,  
       Proline,  
       Prolyl,  
       Serine,  
       Seryl,  
       Threonine,  
       Threonyl,  
       Tryptophan,  
       Tryptophanyl,  
       Tyrosine,  
       Tyrosyl,  
       Valine,  
       Valyl,  
       Selenocysteine,  
       Selenocysteinyl,  
       Pyrrolysine,  
       Pyrrolysyl,  
       Ornithine,  
       Ornithyl,  
       fiveHydroxytryptophan,  
       DOPA,
       Peptide,

       letter3_abbr,
       letter1_abbr,
       alpha_rs

abstract type AbstractPeptide <: AttributedChemical end
abstract type αAminoAcid{DL} <: AbstractPeptide end

struct Alanine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Alanine(; chirality = DLForm, kwargs...) = Core.apply_type(Alanine, chirality)(; kwargs...)
Alanine{DL}(; kwargs...) where {DL} = make_attributedchemical(Alanine{DL}; kwargs...)
struct Alanyl{DL} <: FunctionalGroup{Alanine{DL}, Dehydroxy} end # Ala
struct Arginine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Arginine(; chirality = DLForm, kwargs...) = Core.apply_type(Arginine, chirality)(; kwargs...)
Arginine{DL}(; kwargs...) where {DL} = make_attributedchemical(Arginine{DL}; kwargs...)
struct Arginyl{DL} <: FunctionalGroup{Arginine{DL}, Dehydroxy} end # Arg
struct Aspargine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Aspargine(; chirality = DLForm, kwargs...) = Core.apply_type(Aspargine, chirality)(; kwargs...)
Aspargine{DL}(; kwargs...) where {DL} = make_attributedchemical(Aspargine{DL}; kwargs...)
struct Asparginyl{DL} <: FunctionalGroup{Aspargine{DL}, Dehydroxy} end # Asn
struct Aspartate{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Aspartate(; chirality = DLForm, kwargs...) = Core.apply_type(Aspartate, chirality)(; kwargs...)
Aspartate{DL}(; kwargs...) where {DL} = make_attributedchemical(Aspartate{DL}; kwargs...)
struct Aspartyl{DL} <: FunctionalGroup{Aspartate{DL}, Dehydroxy} end # Asp
struct Cysteine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Cysteine(; chirality = DLForm, kwargs...) = Core.apply_type(Cysteine, chirality)(; kwargs...)
Cysteine{DL}(; kwargs...) where {DL} = make_attributedchemical(Cysteine{DL}; kwargs...)
struct Cysteinyl{DL} <: FunctionalGroup{Cysteine{DL}, Dehydroxy} end # Cys
struct Glutamine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glutamine(; chirality = DLForm, kwargs...) = Core.apply_type(Glutamine, chirality)(; kwargs...)
Glutamine{DL}(; kwargs...) where {DL} = make_attributedchemical(Glutamine{DL}; kwargs...)
struct Glutaminyl{DL} <: FunctionalGroup{Glutamine{DL}, Dehydroxy} end # Gln
struct Glutamate{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glutamate(; chirality = DLForm, kwargs...) = Core.apply_type(Glutamate, chirality)(; kwargs...)
Glutamate{DL}(; kwargs...) where {DL} = make_attributedchemical(Glutamate{DL}; kwargs...)
struct Glutamyl{DL} <: FunctionalGroup{Glutamate{DL}, Dehydroxy} end # Glu
struct Glycine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glycine(; chirality = NoDLForm, kwargs...) = Core.apply_type(Glycine, chirality)(; kwargs...)
Glycine{DL}(; kwargs...) where {DL} = make_attributedchemical(Glycine{DL}; kwargs...)
struct Glycyl{DL} <: FunctionalGroup{Glycine{DL}, Dehydroxy} end # Gly
struct Histidine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Histidine(; chirality = DLForm, kwargs...) = Core.apply_type(Histidine, chirality)(; kwargs...)
Histidine{DL}(; kwargs...) where {DL} = make_attributedchemical(Histidine{DL}; kwargs...)
struct Histidinyl{DL} <: FunctionalGroup{Histidine{DL}, Dehydroxy} end # His
struct Isoleucine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Isoleucine(; chirality = DLForm, kwargs...) = Core.apply_type(Isoleucine, chirality)(; kwargs...)
Isoleucine{DL}(; kwargs...) where {DL} = make_attributedchemical(Isoleucine{DL}; kwargs...)
struct Isoleucyl{DL} <: FunctionalGroup{Isoleucine{DL}, Dehydroxy} end # Ile
struct Leucine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Leucine(; chirality = DLForm, kwargs...) = Core.apply_type(Leucine, chirality)(; kwargs...)
Leucine{DL}(; kwargs...) where {DL} = make_attributedchemical(Leucine{DL}; kwargs...)
struct Leucyl{DL} <: FunctionalGroup{Leucine{DL}, Dehydroxy} end # Leu
struct Lysine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Lysine(; chirality = DLForm, kwargs...) = Core.apply_type(Lysine, chirality)(; kwargs...)
Lysine{DL}(; kwargs...) where {DL} = make_attributedchemical(Lysine{DL}; kwargs...)
struct Lysyl{DL} <: FunctionalGroup{Lysine{DL}, Dehydroxy} end # Lys
struct Methionine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Methionine(; chirality = DLForm, kwargs...) = Core.apply_type(Methionine, chirality)(; kwargs...)
Methionine{DL}(; kwargs...) where {DL} = make_attributedchemical(Methionine{DL}; kwargs...)
struct Methionyl{DL} <: FunctionalGroup{Methionine{DL}, Dehydroxy} end # Met
struct Phenylalanine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Phenylalanine(; chirality = DLForm, kwargs...) = Core.apply_type(Phenylalanine, chirality)(; kwargs...)
Phenylalanine{DL}(; kwargs...) where {DL} = make_attributedchemical(Phenylalanine{DL}; kwargs...)
struct Phenylalanyl{DL} <: FunctionalGroup{Phenylalanine{DL}, Dehydroxy} end # Phe
struct Proline{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Proline(; chirality = DLForm, kwargs...) = Core.apply_type(Proline, chirality)(; kwargs...)
Proline{DL}(; kwargs...) where {DL} = make_attributedchemical(Proline{DL}; kwargs...)
struct Prolyl{DL} <: FunctionalGroup{Proline{DL}, Dehydroxy} end # Pro
struct Serine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Serine(; chirality = DLForm, kwargs...) = Core.apply_type(Serine, chirality)(; kwargs...)
Serine{DL}(; kwargs...) where {DL} = make_attributedchemical(Serine{DL}; kwargs...)
struct Seryl{DL} <: FunctionalGroup{Serine{DL}, Dehydroxy} end # Ser
struct Threonine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Threonine(; chirality = DLForm, kwargs...) = Core.apply_type(Threonine, chirality)(; kwargs...)
Threonine{DL}(; kwargs...) where {DL} = make_attributedchemical(Threonine{DL}; kwargs...)
struct Threonyl{DL} <: FunctionalGroup{Threonine{DL}, Dehydroxy} end # Thr
struct Tryptophan{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Tryptophan(; chirality = DLForm, kwargs...) = Core.apply_type(Tryptophan, chirality)(; kwargs...)
Tryptophan{DL}(; kwargs...) where {DL} = make_attributedchemical(Tryptophan{DL}; kwargs...)
struct Tryptophanyl{DL} <: FunctionalGroup{Tryptophan{DL}, Dehydroxy} end # Trp
struct Tyrosine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Tyrosine(; chirality = DLForm, kwargs...) = Core.apply_type(Tyrosine, chirality)(; kwargs...)
Tyrosine{DL}(; kwargs...) where {DL} = make_attributedchemical(Tyrosine{DL}; kwargs...)
struct Tyrosyl{DL} <: FunctionalGroup{Tyrosine{DL}, Dehydroxy} end # Tyr
struct Valine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Valine(; chirality = DLForm, kwargs...) = Core.apply_type(Valine, chirality)(; kwargs...)
Valine{DL}(; kwargs...) where {DL} = make_attributedchemical(Valine{DL}; kwargs...)
struct Valyl{DL} <: FunctionalGroup{Valine{DL}, Dehydroxy} end # Val
struct Selenocysteine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Selenocysteine(; chirality = DLForm, kwargs...) = Core.apply_type(Selenocysteine, chirality)(; kwargs...)
Selenocysteine{DL}(; kwargs...) where {DL} = make_attributedchemical(Selenocysteine{DL}; kwargs...)
struct Selenocysteinyl{DL} <: FunctionalGroup{Selenocysteine{DL}, Dehydroxy} end # Sec
struct Pyrrolysine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Pyrrolysine(; chirality = DLForm, kwargs...) = Core.apply_type(Pyrrolysine, chirality)(; kwargs...)
Pyrrolysine{DL}(; kwargs...) where {DL} = make_attributedchemical(Pyrrolysine{DL}; kwargs...)
struct Pyrrolysyl{DL} <: FunctionalGroup{Pyrrolysine{DL}, Dehydroxy} end # Pyl
struct Ornithine{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Ornithine(; chirality = DLForm, kwargs...) = Core.apply_type(Ornithine, chirality)(; kwargs...)
Ornithine{DL}(; kwargs...) where {DL} = make_attributedchemical(Ornithine{DL}; kwargs...)
struct Ornithyl{DL} <: FunctionalGroup{Ornithine{DL}, Dehydroxy} end # Orn
struct fiveHydroxytryptophan{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
fiveHydroxytryptophan(; chirality = DLForm, kwargs...) = Core.apply_type(fiveHydroxytryptophan, chirality)(; kwargs...)
fiveHydroxytryptophan{DL}(; kwargs...) where {DL} = make_attributedchemical(fiveHydroxytryptophan{DL}; kwargs...)

struct DOPA{DL} <: αAminoAcid{DL}
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
DOPA(; chirality = DLForm, kwargs...) = Core.apply_type(DOPA, chirality)(; kwargs...)
DOPA{DL}(; kwargs...) where {DL} = make_attributedchemical(DOPA{DL}; kwargs...)

struct Peptide{T} <: AbstractPeptide
    chain::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end

function Peptide(chain; new_attr = true, kwargs...)
    kw = (; pairs(defaultchemicalkwargs(Peptide, chain))..., kwargs...)
    attr = Base.structdiff(kw, NamedTuple{(:name, :elements, :abbreviation)})
    Peptide(chain, String(kw.name), collect(Pair{String, Int}, kw.elements),
            String(kw.abbreviation), buildattribute(new_attr, pairs(attr)))
end

include("base.jl")
include("msc.jl")
include("default.jl")
include("interface.jl")
include("io.jl")

end
