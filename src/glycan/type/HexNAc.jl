"""
    Nacetylhexosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylhexosamine.
"""
struct Nacetylhexosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylhexosamine(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Nacetylhexosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylhexosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylhexosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylglucosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylglucosamine.
"""
struct Nacetylglucosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylglucosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylglucosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylglucosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylglucosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylmannosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylmannosamine.
"""
struct Nacetylmannosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylmannosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylmannosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylmannosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylmannosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylgalactosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylgalactosamine.
"""
struct Nacetylgalactosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylgalactosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylgalactosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylgalactosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylgalactosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylgulosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylgulosamine.
"""
struct Nacetylgulosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylgulosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylgulosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylgulosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylgulosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylaltosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylaltosamine.
"""
struct Nacetylaltosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylaltosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Nacetylaltosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylaltosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylaltosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylallosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylallosamine.
"""
struct Nacetylallosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylallosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylallosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylallosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylallosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetyltalosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetyltalosamine.
"""
struct Nacetyltalosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetyltalosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetyltalosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetyltalosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetyltalosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylidosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylidosamine.
"""
struct Nacetylidosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylidosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Nacetylidosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylidosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylidosamine{D, P, T}, substituent; kwargs...)
