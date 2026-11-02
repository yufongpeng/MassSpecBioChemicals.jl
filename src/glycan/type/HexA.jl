"""
    HexuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Hexuronate.
"""
struct HexuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
HexuronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    HexuronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
HexuronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(HexuronicAcid{D, P, T}, substituent; kwargs...)
"""
    GlucuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Glucuronic acid.
"""
struct GlucuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
GlucuronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    GlucuronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
GlucuronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(GlucuronicAcid{D, P, T}, substituent; kwargs...)
"""
    MannuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Mannuronic acid.
"""
struct MannuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
MannuronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    MannuronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
MannuronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(MannuronicAcid{D, P, T}, substituent; kwargs...)
"""
    GalacturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Galactcuronic acid.
"""
struct GalacturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
GalacturonicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    GalacturonicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
GalacturonicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(GalacturonicAcid{D, P, T}, substituent; kwargs...)
"""
    GuluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

GuluronicAcid.
"""
struct GuluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
GuluronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    GuluronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
GuluronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(GuluronicAcid{D, P, T}, substituent; kwargs...)
"""
    AlturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

AlturonicAcid.
"""
struct AlturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
AlturonicAcid(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    AlturonicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
AlturonicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(AlturonicAcid{D, P, T}, substituent; kwargs...)
"""
    AlluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

AlluronicAcid.
"""
struct AlluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
AlluronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    AlluronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
AlluronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(AlluronicAcid{D, P, T}, substituent; kwargs...)
"""
    TaluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

TaluronicAcid.
"""
struct TaluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
TaluronicAcid(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    TaluronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
TaluronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(TaluronicAcid{D, P, T}, substituent; kwargs...)
"""
    IduronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

IduronicAcid.
"""
struct IduronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
IduronicAcid(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    IduronicAcid{DL, LFP, typeof(substituent)}(substituent; kwargs...)
IduronicAcid{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(IduronicAcid{D, P, T}, substituent; kwargs...)
