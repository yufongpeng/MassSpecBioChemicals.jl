"""
    Dideoxyhexose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Dideoxyhexose.
"""
struct Dideoxyhexose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
end
Dideoxyhexose(s::T; LFP = PyranoseForm, DL = DLForm) where T = Dideoxyhexose{DL, LFP, T}(s)
Dideoxyhexose(; LFP = PyranoseForm, DL = DLForm) = Dideoxyhexose{DL, LFP, Nothing}(nothing)
"""
    Olivose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Olivose.
"""
struct Olivose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
end
Olivose(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Olivose{DL, LFP, T}(s)
Olivose(; LFP = PyranoseForm, DL = NoDLForm) = Olivose{DL, LFP, Nothing}(nothing)
# Oli = 2,6dGlc/Man
"""
    Tyvelose{D, P, T} <: AbstractDideoxyhexose{D, P, T}

Tyvelose.
"""
struct Tyvelose{D, P, T} <: AbstractDideoxyhexose{D, P, T}
    substituent::T
end
Tyvelose(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Tyvelose{DL, LFP, T}(s)
Tyvelose(; LFP = PyranoseForm, DL = NoDLForm) = Tyvelose{DL, LFP, Nothing}(nothing)
# Tyv = 3,6dMan/Alt(D)
"""
    Abequose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Abequose.
"""
struct Abequose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Abequose(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Abequose{DL, LFP, T}(s)
Abequose(; LFP = PyranoseForm, DL = NoDLForm) = Abequose{DL, LFP, Nothing}(nothing)
# Abe = 3,6dGul/Gal
"""
    Colitose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Colitose.
"""
struct Colitose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Colitose(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Colitose{DL, LFP, T}(s)
Colitose(; LFP = PyranoseForm, DL = NoDLForm) = Colitose{DL, LFP, Nothing}(nothing)
# Col = L-Abe
"""
    Paratose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Paratose.
"""
struct Paratose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Paratose(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Paratose{DL, LFP, T}(s)
Paratose(; LFP = PyranoseForm, DL = NoDLForm) = Paratose{DL, LFP, Nothing}(nothing)
# Par = 3,6dGlc/All 
"""
    Digitoxose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Digitoxose.
"""
struct Digitoxose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Digitoxose(s::T; LFP = PyranoseForm, DL = DForm) where T = Digitoxose{DL, LFP, T}(s)
Digitoxose(; LFP = PyranoseForm, DL = DForm) = Digitoxose{DL, LFP, Nothing}(nothing)
# Dig = 2,6dAll/Alt(D) or L
