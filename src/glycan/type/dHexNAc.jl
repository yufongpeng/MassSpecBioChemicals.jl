"""
    Nacetyldeoxyhexosamine{D, P, T} <: AbstractNacetyldeoxyhexosamine{D, P, T}

N-acetyldeoxyhexosamine.
"""
struct Nacetyldeoxyhexosamine{D, P, T} <: AbstractNacetyldeoxyhexosamine{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetyldeoxyhexosamine(substituent = nothing; LFP = PyranoseForm, DL = DLForm, kwargs...) =
    Nacetyldeoxyhexosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetyldeoxyhexosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetyldeoxyhexosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylquinovosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylquinovose.
"""
struct Nacetylquinovosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylquinovosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylquinovosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylquinovosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylquinovosamine{D, P, T}, substituent; kwargs...)
# Qui = 6dGlc
# SQui = 6SQui
"""
    Nacetylrhamnosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylrhamnosamine.
"""
struct Nacetylrhamnosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylrhamnosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Nacetylrhamnosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylrhamnosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylrhamnosamine{D, P, T}, substituent; kwargs...)
# Rha = 6dMan(L)
"""
    Nacetylfucosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylfucosamine.
"""
struct Nacetylfucosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylfucosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Nacetylfucosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylfucosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylfucosamine{D, P, T}, substituent; kwargs...)
# Fuc = 6dGal(L)
"""
    Nacetylsixdeoxyaltosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetyl-6-deoxyaltosamine.
"""
struct Nacetylsixdeoxyaltosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylsixdeoxyaltosamine(substituent = nothing; LFP = PyranoseForm, DL = LForm, kwargs...) =
    Nacetylsixdeoxyaltosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylsixdeoxyaltosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylsixdeoxyaltosamine{D, P, T}, substituent; kwargs...)
"""
    Nacetylsixdeoxytalosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetyl-6-deoxytalosamine.
"""
struct Nacetylsixdeoxytalosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
Nacetylsixdeoxytalosamine(substituent = nothing; LFP = PyranoseForm, DL = DForm, kwargs...) =
    Nacetylsixdeoxytalosamine{DL, LFP, typeof(substituent)}(substituent; kwargs...)
Nacetylsixdeoxytalosamine{D, P, T}(substituent = nothing; kwargs...) where {D, P, T} =
    make_monosaccharide(Nacetylsixdeoxytalosamine{D, P, T}, substituent; kwargs...)
