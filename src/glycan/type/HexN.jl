"""
    Hexosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Hexosamine.
"""
struct Hexosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
end
Hexosamine(s::T; LFP = PyranoseForm, DL = DLForm) where T = Hexosamine{DL, LFP, T}(s)
Hexosamine(; LFP = PyranoseForm, DL = DLForm) = Hexosamine{DL, LFP, Nothing}(nothing)
"""
    Glucosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Glucosamine.
"""
struct Glucosamine{D, P, T} <: AbstractHexosamine{D, P, T}
    substituent::T
end
Glucosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Glucosamine{DL, LFP, T}(s)
Glucosamine(; LFP = PyranoseForm, DL = DForm) = Glucosamine{DL, LFP, Nothing}(nothing)
"""
    Mannosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Mannosamine.
"""
struct Mannosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Mannosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Mannosamine{DL, LFP, T}(s)
Mannosamine(; LFP = PyranoseForm, DL = DForm) = Mannosamine{DL, LFP, Nothing}(nothing)
"""
    Galactosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Galactosamine.
"""
struct Galactosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Galactosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Galactosamine{DL, LFP, T}(s)
Galactosamine(; LFP = PyranoseForm, DL = DForm) = Galactosamine{DL, LFP, Nothing}(nothing)
"""
    Gulosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Gulosamine.
"""
struct Gulosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Gulosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Gulosamine{DL, LFP, T}(s)
Gulosamine(; LFP = PyranoseForm, DL = DForm) = Gulosamine{DL, LFP, Nothing}(nothing)
"""
    Altosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Altosamine.
"""
struct Altosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Altosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Altosamine{DL, LFP, T}(s)
Altosamine(; LFP = PyranoseForm, DL = LForm) = Altosamine{DL, LFP, Nothing}(nothing)
"""
    Allosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Allosamine.
"""
struct Allosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Allosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Allosamine{DL, LFP, T}(s)
Allosamine(; LFP = PyranoseForm, DL = DForm) = Allosamine{DL, LFP, Nothing}(nothing)
"""
    Talosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Talosamine.
"""
struct Talosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Talosamine(s::T; LFP = PyranoseForm, DL = DForm) where T = Talosamine{DL, LFP, T}(s)
Talosamine(; LFP = PyranoseForm, DL = DForm) = Talosamine{DL, LFP, Nothing}(nothing)
"""
    Idosamine{D, P, T} <: AbstractHexosamine{D, P, T}

Idosamine.
"""
struct Idosamine{D, P, T} <: AbstractHexosamine{D, P, T}    
    substituent::T
end
Idosamine(s::T; LFP = PyranoseForm, DL = LForm) where T = Idosamine{DL, LFP, T}(s)
Idosamine(; LFP = PyranoseForm, DL = LForm) = Idosamine{DL, LFP, Nothing}(nothing)
