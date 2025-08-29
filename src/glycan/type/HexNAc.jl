"""
    Nacetylhexosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylhexosamine.
"""
struct Nacetylhexosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
end
Nacetylhexosamine(s::T; LFP = PyranoseForm, DL = DLForm) where T = Nacetylhexosamine{DL, LFP, T}(s)
Nacetylhexosamine(; LFP = PyranoseForm, DL = DLForm) = Nacetylhexosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylglucosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylglucosamine.
"""
struct Nacetylglucosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}
    substituent::T
end
Nacetylglucosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylglucosamine{DL, LFP, T}(s)
Nacetylglucosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylglucosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylmannosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylmannosamine.
"""
struct Nacetylmannosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylmannosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylmannosamine{DL, LFP, T}(s)
Nacetylmannosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylmannosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylgalactosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylgalactosamine.
"""
struct Nacetylgalactosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylgalactosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylgalactosamine{DL, LFP, T}(s)
Nacetylgalactosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylgalactosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylgulosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylgulosamine.
"""
struct Nacetylgulosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylgulosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylgulosamine{DL, LFP, T}(s)
Nacetylgulosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylgulosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylaltosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylaltosamine.
"""
struct Nacetylaltosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylaltosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Nacetylaltosamine{DL, LFP, T}(s)
Nacetylaltosamine(; LFP = PyranoseForm, DL = LForm) = Nacetylaltosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylallosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylallosamine.
"""
struct Nacetylallosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylallosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetylallosamine{DL, LFP, T}(s)
Nacetylallosamine(; LFP = PyranoseForm, DL = DForm) = Nacetylallosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetyltalosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetyltalosamine.
"""
struct Nacetyltalosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetyltalosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Nacetyltalosamine{DL, LFP, T}(s)
Nacetyltalosamine(; LFP = PyranoseForm, DL = DForm) = Nacetyltalosamine{DL, LFP, Nothing}(nothing)
"""
    Nacetylidosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}

N-acetylidosamine.
"""
struct Nacetylidosamine{D, P, T} <: AbstractNacetylhexosamine{D, P, T}    
    substituent::T
end
Nacetylidosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Nacetylidosamine{DL, LFP, T}(s)
Nacetylidosamine(; LFP = PyranoseForm, DL = LForm) = Nacetylidosamine{DL, LFP, Nothing}(nothing)
