"""
    Deoxyhexose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Deoxyhexose.
"""
struct Deoxyhexose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Deoxyhexose(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Deoxyhexose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Deoxyhexose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Deoxyhexose{D, P, T}, substituent; kwargs...)
"""
    Quinovose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Quinovose.
"""
struct Quinovose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Quinovose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Quinovose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Quinovose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Quinovose{D, P, T}, substituent; kwargs...)
# Qui = 6dGlc
# SQui = 6SQui
"""
    Rhamnose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Rhamnose.
"""
struct Rhamnose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Rhamnose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Rhamnose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Rhamnose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Rhamnose{D, P, T}, substituent; kwargs...)
# Rha = 6dMan(L)
"""
    Fucose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Fucose.
"""
struct Fucose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Fucose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Fucose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Fucose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Fucose{D, P, T}, substituent; kwargs...)
# Fuc = 6dGal(L)
"""
    Sixdeoxygulose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxygulose.
"""
struct Sixdeoxygulose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Sixdeoxygulose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Sixdeoxygulose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Sixdeoxygulose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Sixdeoxygulose{D, P, T}, substituent; kwargs...)
# 6dGul (no NAc)
"""
    Sixdeoxyaltose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxyaltose.
"""
struct Sixdeoxyaltose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Sixdeoxyaltose(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Sixdeoxyaltose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Sixdeoxyaltose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Sixdeoxyaltose{D, P, T}, substituent; kwargs...)
"""
    Sixdeoxytalose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxytalose.
"""
struct Sixdeoxytalose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Sixdeoxytalose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Sixdeoxytalose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Sixdeoxytalose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Sixdeoxytalose{D, P, T}, substituent; kwargs...)
