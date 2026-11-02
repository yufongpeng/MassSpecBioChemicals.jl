"""
    Pentose{D, P, T} <: AbstractPentose{D, P, T}

Pentose.
"""
struct Pentose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Pentose(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Pentose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Pentose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Pentose{D, P, T}, substituent; kwargs...)
"""
    Arabinose{D, P, T} <: AbstractPentose{D, P, T}

Arabinose.
"""
struct Arabinose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Arabinose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Arabinose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Arabinose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Arabinose{D, P, T}, substituent; kwargs...)
"""
    Lyxose{D, P, T} <: AbstractPentose{D, P, T}

Lyxose.
"""
struct Lyxose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Lyxose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Lyxose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Lyxose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Lyxose{D, P, T}, substituent; kwargs...)
"""
    Xylose{D, P, T} <: AbstractPentose{D, P, T}

Xylose.
"""
struct Xylose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Xylose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Xylose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Xylose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Xylose{D, P, T}, substituent; kwargs...)
"""
    Ribose{D, P, T} <: AbstractPentose{D, P, T}

Ribose.
"""
struct Ribose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Ribose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Ribose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Ribose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Ribose{D, P, T}, substituent; kwargs...)

"""
    Deoxypentose{D, P, T} <: AbstractDeoxypentose{D, P, T}

Deoxypentose.
"""
struct Deoxypentose{D, P, T} <: AbstractDeoxypentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Deoxypentose(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Deoxypentose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Deoxypentose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Deoxypentose{D, P, T}, substituent; kwargs...)
"""
    Deoxyribose{D, P, T} <: AbstractDeoxypentose{D, P, T}

Deoxyribose.
"""
struct Deoxyribose{D, P, T} <: AbstractDeoxypentose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Deoxyribose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Deoxyribose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Deoxyribose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Deoxyribose{D, P, T}, substituent; kwargs...)
