"""
    NeuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

Neuraminic acid.
"""
struct NeuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
NeuraminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    NeuraminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
NeuraminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(NeuraminicAcid{D, P, T}, substituent; kwargs...)
"""
    NacetylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

N-acetylneuraminic acid (nana).
"""
struct NacetylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
NacetylneuraminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    NacetylneuraminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
NacetylneuraminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(NacetylneuraminicAcid{D, P, T}, substituent; kwargs...)
"""
    NglycolylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

N-glycolylneuraminic acid.
"""
struct NglycolylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
NglycolylneuraminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    NglycolylneuraminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
NglycolylneuraminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(NglycolylneuraminicAcid{D, P, T}, substituent; kwargs...)
"""
    Kdn{D, P, T} <: SialicAcid{D, P, T}

2-Keto-3-deoxy-nononic acid. 
"""
struct Kdn{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Kdn(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Kdn{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Kdn{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Kdn{D, P, T}, substituent; kwargs...)
# Kdn, Neu = Kdn5N

# 9dKdn5,7N: Pse, Leg, Aci, 4eLeg, 8eLeg, 8eAci

"""
   DiaminodideoxynonulosonicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

Diamino-dideoxy-nonulosonicAcid.
"""
struct DiaminodideoxynonulosonicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
DiaminodideoxynonulosonicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    DiaminodideoxynonulosonicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
DiaminodideoxynonulosonicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(DiaminodideoxynonulosonicAcid{D, P, T}, substituent; kwargs...)
"""
   LegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

Legionaminic acid.
"""
struct LegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
LegionaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    LegionaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
LegionaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(LegionaminicAcid{D, P, T}, substituent; kwargs...)
# Leg = 9dKdn5,7N
"""
   FourepilegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

4-epi-legionaminic acid.
"""
struct FourepilegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
FourepilegionaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    FourepilegionaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
FourepilegionaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(FourepilegionaminicAcid{D, P, T}, substituent; kwargs...)
"""
   EightepilegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

8-epi-legionaminic acid.
"""
struct EightepilegionaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
EightepilegionaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    EightepilegionaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
EightepilegionaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(EightepilegionaminicAcid{D, P, T}, substituent; kwargs...)
"""
   AcinetaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

Acinetaminic acid.
"""
struct AcinetaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
AcinetaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    AcinetaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
AcinetaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(AcinetaminicAcid{D, P, T}, substituent; kwargs...)
# Aci = 9d,7,8eKdn5,7N
"""
   EightepiacinetaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

8-epei-acinetaminic acid.
"""
struct EightepiacinetaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
EightepiacinetaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    EightepiacinetaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
EightepiacinetaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(EightepiacinetaminicAcid{D, P, T}, substituent; kwargs...)
"""
   PseudaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}

Pseudaminic acid.
"""
struct PseudaminicAcid{D, P, T} <: AbstractdiaminodideoxynonulosonicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
PseudaminicAcid(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    PseudaminicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
PseudaminicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(PseudaminicAcid{D, P, T}, substituent; kwargs...)
# Pse = 9d,5,7,8eKdn5,7N
