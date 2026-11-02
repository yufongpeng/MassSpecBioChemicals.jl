"""
    Dideoxyhexose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Dideoxyhexose.
"""
struct Dideoxyhexose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Dideoxyhexose(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Dideoxyhexose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Dideoxyhexose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Dideoxyhexose{D, P, T}, substituent; kwargs...)
"""
    Olivose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Olivose.
"""
struct Olivose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Olivose(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Olivose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Olivose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Olivose{D, P, T}, substituent; kwargs...)
# Oli = 2,6dGlc/Man
"""
    Tyvelose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Tyvelose.
"""
struct Tyvelose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Tyvelose(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Tyvelose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Tyvelose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Tyvelose{D, P, T}, substituent; kwargs...)
# Tyv = 3,6dMan/Alt(D)
"""
    Abequose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Abequose.
"""
struct Abequose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Abequose(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Abequose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Abequose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Abequose{D, P, T}, substituent; kwargs...)
# Abe = 3,6dGul/Gal
"""
    Colitose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Colitose.
"""
struct Colitose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Colitose(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Colitose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Colitose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Colitose{D, P, T}, substituent; kwargs...)
# Col = L-Abe
"""
    Paratose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Paratose.
"""
struct Paratose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Paratose(substituent = nothing; LFP = PyranoseForm, DL = NoDLForm, kwargs...) =
    Paratose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Paratose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Paratose{D, P, T}, substituent; kwargs...)
# Par = 3,6dGlc/All 
"""
    Digitoxose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Digitoxose.
"""
struct Digitoxose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Digitoxose(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Digitoxose{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Digitoxose{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Digitoxose{D, P, T}, substituent; kwargs...)
# Dig = 2,6dAll/Alt(D) or L
