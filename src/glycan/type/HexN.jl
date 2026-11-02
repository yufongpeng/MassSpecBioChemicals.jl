"""
    Hexosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Hexosamine.
"""
struct Hexosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Hexosamine(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Hexosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Hexosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Hexosamine{D, P, T}, substituent; kwargs...)
"""
    Glucosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Glucosamine.
"""
struct Glucosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Glucosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Glucosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Glucosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Glucosamine{D, P, T}, substituent; kwargs...)
"""
    Mannosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Mannosamine.
"""
struct Mannosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Mannosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Mannosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Mannosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Mannosamine{D, P, T}, substituent; kwargs...)
"""
    Galactosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Galactosamine.
"""
struct Galactosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Galactosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Galactosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Galactosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Galactosamine{D, P, T}, substituent; kwargs...)
"""
    Gulosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Gulosamine.
"""
struct Gulosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Gulosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Gulosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Gulosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Gulosamine{D, P, T}, substituent; kwargs...)
"""
    Altosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Altosamine.
"""
struct Altosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Altosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Altosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Altosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Altosamine{D, P, T}, substituent; kwargs...)
"""
    Allosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Allosamine.
"""
struct Allosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Allosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Allosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Allosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Allosamine{D, P, T}, substituent; kwargs...)
"""
    Talosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Talosamine.
"""
struct Talosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Talosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Talosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Talosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Talosamine{D, P, T}, substituent; kwargs...)
"""
    Idosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Idosamine.
"""
struct Idosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Idosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Idosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Idosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Idosamine{D, P, T}, substituent; kwargs...)
