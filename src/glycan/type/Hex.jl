"""
    Hexose{D, P, T} <: AbstractHexose{D, P, T}

Hexose.
"""
struct Hexose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Hexose(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Hexose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Hexose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Hexose{D, P, T}, substituent; kwargs...)
"""
    Glucose{D, P, T} <: AbstractHexose{D, P, T}

Glucose.
"""
struct Glucose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glucose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Glucose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Glucose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Glucose{D, P, T}, substituent; kwargs...)
"""
    Mannose{D, P, T} <: AbstractHexose{D, P, T}

Mannose.
"""
struct Mannose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Mannose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Mannose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Mannose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Mannose{D, P, T}, substituent; kwargs...)
"""
    Galactose{D, P, T} <: AbstractHexose{D, P, T}

Galactose.
"""
struct Galactose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Galactose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Galactose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Galactose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Galactose{D, P, T}, substituent; kwargs...)
"""
    Gulose{D, P, T} <: AbstractHexose{D, P, T}

Gulose.
"""
struct Gulose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Gulose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Gulose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Gulose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Gulose{D, P, T}, substituent; kwargs...)
"""
    Altose{D, P, T} <: AbstractHexose{D, P, T}

Altose.
"""
struct Altose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Altose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Altose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Altose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Altose{D, P, T}, substituent; kwargs...)
"""
    Allose{D, P, T} <: AbstractHexose{D, P, T}

Allose.
"""
struct Allose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Allose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Allose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Allose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Allose{D, P, T}, substituent; kwargs...)
"""
    Talose{D, P, T} <: AbstractHexose{D, P, T}

Talose.
"""
struct Talose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Talose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Talose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Talose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Talose{D, P, T}, substituent; kwargs...)
"""
    Idose{D, P, T} <: AbstractHexose{D, P, T}

Idose.
"""
struct Idose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Idose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Idose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Idose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Idose{D, P, T}, substituent; kwargs...)

# 3-MethylMethyl : Api, Fru, Tag, Sor, Psi
"""
    Apiose{D, P, T} <: AbstractHexose{D, P, T}

Apiose.
"""
struct Apiose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Apiose(substituent = nothing; LFP = FuranoseForm, DL = LForm, kwargs...) =
    Apiose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Apiose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Apiose{D, P, T}, substituent; kwargs...)
"""
    Fructose{D, P, T} <: AbstractHexose{D, P, T}

Fructose.
"""
struct Fructose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Fructose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Fructose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Fructose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Fructose{D, P, T}, substituent; kwargs...)
"""
    Tagatose{D, P, T} <: AbstractHexose{D, P, T}

Tagatose.
"""
struct Tagatose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Tagatose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Tagatose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Tagatose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Tagatose{D, P, T}, substituent; kwargs...)
"""
    Sorbose{D, P, T} <: AbstractHexose{D, P, T}

Sorbose.
"""
struct Sorbose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Sorbose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Sorbose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Sorbose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Sorbose{D, P, T}, substituent; kwargs...)
"""
    Psicose{D, P, T} <: AbstractHexose{D, P, T}

Psicose.
"""
struct Psicose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Psicose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Psicose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Psicose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Psicose{D, P, T}, substituent; kwargs...)
