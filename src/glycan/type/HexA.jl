"""
    HexuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Hexuronate.
"""
struct HexuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
end
HexuronicAcid(s::T; LFP = PyranoseForm, DL = DLForm) where T = HexuronicAcid{DL, LFP, T}(s)
HexuronicAcid(; LFP = PyranoseForm, DL = DLForm) = HexuronicAcid{DL, LFP, Nothing}(nothing)
"""
    GlucuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Glucuronic acid.
"""
struct GlucuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
end
GlucuronicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = GlucuronicAcid{DL, LFP, T}(s)
GlucuronicAcid(; LFP = PyranoseForm, DL = DForm) = GlucuronicAcid{DL, LFP, Nothing}(nothing)
"""
    MannuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Mannuronic acid.
"""
struct MannuronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
end
MannuronicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = MannuronicAcid{DL, LFP, T}(s)
MannuronicAcid(; LFP = PyranoseForm, DL = DForm) = MannuronicAcid{DL, LFP, Nothing}(nothing)
"""
    GalacturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

Galactcuronic acid.
"""
struct GalacturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}
    substituent::T
end
GalacturonicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = GalacturonicAcid{DL, LFP, T}(s)
GalacturonicAcid(; LFP = PyranoseForm, DL = DForm) = GalacturonicAcid{DL, LFP, Nothing}(nothing)
"""
    GuluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

GuluronicAcid.
"""
struct GuluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}    
    substituent::T
end
GuluronicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = GuluronicAcid{DL, LFP, T}(s)
GuluronicAcid(; LFP = PyranoseForm, DL = DForm) = GuluronicAcid{DL, LFP, Nothing}(nothing)
"""
    AlturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

AlturonicAcid.
"""
struct AlturonicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}    
    substituent::T
end
AlturonicAcid(s::T; LFP = PyranoseForm, DL = LForm) where T = AlturonicAcid{DL, LFP, T}(s)
AlturonicAcid(; LFP = PyranoseForm, DL = LForm) = AlturonicAcid{DL, LFP, Nothing}(nothing)
"""
    AlluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

AlluronicAcid.
"""
struct AlluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}    
    substituent::T
end
AlluronicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = AlluronicAcid{DL, LFP, T}(s)
AlluronicAcid(; LFP = PyranoseForm, DL = DForm) = AlluronicAcid{DL, LFP, Nothing}(nothing)
"""
    TaluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

TaluronicAcid.
"""
struct TaluronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}    
    substituent::T
end
TaluronicAcid(s::T; LFP = PyranoseForm, DL = DForm) where T = TaluronicAcid{DL, LFP, T}(s)
TaluronicAcid(; LFP = PyranoseForm, DL = DForm) = TaluronicAcid{DL, LFP, Nothing}(nothing)
"""
    IduronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}

IduronicAcid.
"""
struct IduronicAcid{D, P, T} <: AbstractHexuronicAcid{D, P, T}    
    substituent::T
end
IduronicAcid(s::T; LFP = PyranoseForm, DL = LForm) where T = IduronicAcid{DL, LFP, T}(s)
IduronicAcid(; LFP = PyranoseForm, DL = LForm) = IduronicAcid{DL, LFP, Nothing}(nothing)
