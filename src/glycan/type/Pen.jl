"""
    Pentose{D, P, T} <: AbstractPentose{D, P, T}

Pentose.
"""
struct Pentose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
end
Pentose(s::T; LFP = PyranoseForm, DL = DLForm) where T = Pentose{DL, LFP, T}(s)
Pentose(; LFP = PyranoseForm, DL = DLForm) = Pentose{DL, LFP, Nothing}(nothing)
"""
    Arabinose{D, P, T} <: AbstractPentose{D, P, T}

Arabinose.
"""
struct Arabinose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
end
Arabinose(s::T; LFP = PyranoseForm, DL = LForm) where T = Arabinose{DL, LFP, T}(s)
Arabinose(; LFP = PyranoseForm, DL = LForm) = Arabinose{DL, LFP, Nothing}(nothing)
"""
    Lyxose{D, P, T} <: AbstractPentose{D, P, T}

Lyxose.
"""
struct Lyxose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
end
Lyxose(s::T; LFP = PyranoseForm, DL = DForm) where T = Lyxose{DL, LFP, T}(s)
Lyxose(; LFP = PyranoseForm, DL = DForm) = Lyxose{DL, LFP, Nothing}(nothing)
"""
    Xylose{D, P, T} <: AbstractPentose{D, P, T}

Xylose.
"""
struct Xylose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
end
Xylose(s::T; LFP = PyranoseForm, DL = DForm) where T = Xylose{DL, LFP, T}(s)
Xylose(; LFP = PyranoseForm, DL = DForm) = Xylose{DL, LFP, Nothing}(nothing)
"""
    Ribose{D, P, T} <: AbstractPentose{D, P, T}

Ribose.
"""
struct Ribose{D, P, T} <: AbstractPentose{D, P, T}
    substituent::T
end
Ribose(s::T; LFP = PyranoseForm, DL = DForm) where T = Ribose{DL, LFP, T}(s)
Ribose(; LFP = PyranoseForm, DL = DForm) = Ribose{DL, LFP, Nothing}(nothing)

"""
    Deoxypentose{D, P, T} <: AbstractDeoxypentose{D, P, T}

Deoxypentose.
"""
struct Deoxypentose{D, P, T} <: AbstractDeoxypentose{D, P, T}
    substituent::T
end
Deoxypentose(s::T; LFP = PyranoseForm, DL = DLForm) where T = Deoxypentose{DL, LFP, T}(s)
Deoxypentose(; LFP = PyranoseForm, DL = DLForm) = Deoxypentose{DL, LFP, Nothing}(nothing)
"""
    Deoxyribose{D, P, T} <: AbstractDeoxypentose{D, P, T}

Deoxyribose.
"""
struct Deoxyribose{D, P, T} <: AbstractDeoxypentose{D, P, T}
    substituent::T
end
Deoxyribose(s::T; LFP = PyranoseForm, DL = DForm) where T = Deoxyribose{DL, LFP, T}(s)
Deoxyribose(; LFP = PyranoseForm, DL = DForm) = Deoxyribose{DL, LFP, Nothing}(nothing)
