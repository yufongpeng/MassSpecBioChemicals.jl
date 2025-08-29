"""
    Nacetyldeoxyhexosamine{D, P, T} <: AbstractNacetyldeoxyhexosamine{D, P, T}

N-acetyldeoxyhexosamine.
"""
struct Nacetyldeoxyhexosamine{D, P, T} <: AbstractNacetyldeoxyhexosamine{D, P, T}
    substituent::T
end
Nacetyldeoxyhexosamine(s::T; LFP = PyranoseForm, DL = DLForm) where T = Nacetyldeoxyhexosamine{DL, LFP, T}(s)
Nacetyldeoxyhexosamine(; LFP = PyranoseForm, DL = DLForm) = Nacetyldeoxyhexosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylquinovosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylquinovose.
"""
struct Nacetylquinovosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Nacetylquinovosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylquinovosamine{DL, LFP, T}(s)
Nacetylquinovosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylquinovosamine{DL, LFP, Nothing}(nothing)
# Qui = 6dGlc
# SQui = 6SQui
"""
    Nacetylrhamnosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylrhamnosamine.
"""
struct Nacetylrhamnosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Nacetylrhamnosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Nacetylrhamnosamine{DL, LFP, T}(s)
Nacetylrhamnosamine(; LFP = PyranoseForm, DL = LForm) = Nacetylrhamnosamine{DL, LFP, Nothing}(nothing)
# Rha = 6dMan(L)
"""
    Nacetylfucosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetylfucosamine.
"""
struct Nacetylfucosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Nacetylfucosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Nacetylfucosamine{DL, LFP, T}(s)
Nacetylfucosamine(; LFP = PyranoseForm, DL = LForm) = Nacetylfucosamine{DL, LFP, Nothing}(nothing)
# Fuc = 6dGal(L)
"""
    Nacetylsixdeoxyaltosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetyl-6-deoxyaltosamine.
"""
struct Nacetylsixdeoxyaltosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Nacetylsixdeoxyaltosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Nacetylsixdeoxyaltosamine{DL, LFP, T}(s)
Nacetylsixdeoxyaltosamine(; LFP = PyranoseForm, DL = LForm) = Nacetylsixdeoxyaltosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylsixdeoxytalosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}

N-acetyl-6-deoxytalosamine.
"""
struct Nacetylsixdeoxytalosamine{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Nacetylsixdeoxytalosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylsixdeoxytalosamine{DL, LFP, T}(s)
Nacetylsixdeoxytalosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylsixdeoxytalosamine{DL, LFP, Nothing}(nothing)
