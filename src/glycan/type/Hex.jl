"""
    Hexose{D, P, T} <: AbstractHexose{D, P, T}

Hexose.
"""
struct Hexose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
end
Hexose(s::T; LFP = PyranoseForm, DL = DLForm) where T = Hexose{DL, LFP, T}(s)
Hexose(; LFP = PyranoseForm, DL = DLForm) = Hexose{DL, LFP, Nothing}(nothing)
"""
    Glucose{D, P, T} <: AbstractHexose{D, P, T}

Glucose.
"""
struct Glucose{D, P, T} <: AbstractHexose{D, P, T}
    substituent::T
end
Glucose(s::T; LFP = PyranoseForm, DL = DForm) where T = Glucose{DL, LFP, T}(s)
Glucose(; LFP = PyranoseForm, DL = DForm) = Glucose{DL, LFP, Nothing}(nothing)
"""
    Mannose{D, P, T} <: AbstractHexose{D, P, T}

Mannose.
"""
struct Mannose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Mannose(s::T; LFP = PyranoseForm, DL = DForm) where T = Mannose{DL, LFP, T}(s)
Mannose(; LFP = PyranoseForm, DL = DForm) = Mannose{DL, LFP, Nothing}(nothing)
"""
    Galactose{D, P, T} <: AbstractHexose{D, P, T}

Galactose.
"""
struct Galactose{D, P, T} <: AbstractHexose{D, P, T} 
    substituent::T
end
Galactose(s::T; LFP = PyranoseForm, DL = DForm) where T = Galactose{DL, LFP, T}(s)
Galactose(; LFP = PyranoseForm, DL = DForm) = Galactose{DL, LFP, Nothing}(nothing)
"""
    Gulose{D, P, T} <: AbstractHexose{D, P, T}

Gulose.
"""
struct Gulose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Gulose(s::T; LFP = PyranoseForm, DL = DForm) where T = Gulose{DL, LFP, T}(s)
Gulose(; LFP = PyranoseForm, DL = DForm) = Gulose{DL, LFP, Nothing}(nothing)
"""
    Altose{D, P, T} <: AbstractHexose{D, P, T}

Altose.
"""
struct Altose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Altose(s::T; LFP = PyranoseForm, DL = LForm) where T = Altose{DL, LFP, T}(s)
Altose(; LFP = PyranoseForm, DL = LForm) = Altose{DL, LFP, Nothing}(nothing)
"""
    Allose{D, P, T} <: AbstractHexose{D, P, T}

Allose.
"""
struct Allose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Allose(s::T; LFP = PyranoseForm, DL = DForm) where T = Allose{DL, LFP, T}(s)
Allose(; LFP = PyranoseForm, DL = DForm) = Allose{DL, LFP, Nothing}(nothing)
"""
    Talose{D, P, T} <: AbstractHexose{D, P, T}

Talose.
"""
struct Talose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Talose(s::T; LFP = PyranoseForm, DL = DForm) where T = Talose{DL, LFP, T}(s)
Talose(; LFP = PyranoseForm, DL = DForm) = Talose{DL, LFP, Nothing}(nothing)
"""
    Idose{D, P, T} <: AbstractHexose{D, P, T}

Idose.
"""
struct Idose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Idose(s::T; LFP = PyranoseForm, DL = LForm) where T = Idose{DL, LFP, T}(s)
Idose(; LFP = PyranoseForm, DL = LForm) = Idose{DL, LFP, Nothing}(nothing)

# 3-MethylMethyl : Api, Fru, Tag, Sor, Psi
"""
    Apiose{D, P, T} <: AbstractHexose{D, P, T}

Apiose.
"""
struct Apiose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Apiose(s::T; LFP = FuranoseForm, DL = LForm) where T = Apiose{DL, LFP, T}(s)
Apiose(; LFP = FuranoseForm, DL = LForm) = Apiose{DL, LFP, Nothing}(nothing)
"""
    Fructose{D, P, T} <: AbstractHexose{D, P, T}

Fructose.
"""
struct Fructose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Fructose(s::T; LFP = PyranoseForm, DL = LForm) where T = Fructose{DL, LFP, T}(s)
Fructose(; LFP = PyranoseForm, DL = LForm) = Fructose{DL, LFP, Nothing}(nothing)
"""
    Tagatose{D, P, T} <: AbstractHexose{D, P, T}

Tagatose.
"""
struct Tagatose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Tagatose(s::T; LFP = PyranoseForm, DL = DForm) where T = Tagatose{DL, LFP, T}(s)
Tagatose(; LFP = PyranoseForm, DL = DForm) = Tagatose{DL, LFP, Nothing}(nothing)
"""
    Sorbose{D, P, T} <: AbstractHexose{D, P, T}

Sorbose.
"""
struct Sorbose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Sorbose(s::T; LFP = PyranoseForm, DL = LForm) where T = Sorbose{DL, LFP, T}(s)
Sorbose(; LFP = PyranoseForm, DL = LForm) = Sorbose{DL, LFP, Nothing}(nothing)
"""
    Psicose{D, P, T} <: AbstractHexose{D, P, T}

Psicose.
"""
struct Psicose{D, P, T} <: AbstractHexose{D, P, T}    
    substituent::T
end
Psicose(s::T; LFP = PyranoseForm, DL = DForm) where T = Psicose{DL, LFP, T}(s)
Psicose(; LFP = PyranoseForm, DL = DForm) = Psicose{DL, LFP, Nothing}(nothing)
