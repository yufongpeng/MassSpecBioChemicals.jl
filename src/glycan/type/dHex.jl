"""
    Deoxyhexose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Deoxyhexose.
"""
struct Deoxyhexose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Deoxyhexose(s::T; LFP = PyranoseForm, DL = DLForm) where T = Deoxyhexose{DL, LFP, T}(s)
Deoxyhexose(; LFP = PyranoseForm, DL = DLForm) = Deoxyhexose{DL, LFP, Nothing}(nothing)
"""
    Quinovose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Quinovose.
"""
struct Quinovose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Quinovose(s::T; LFP = PyranoseForm, DL = DForm) where T = Quinovose{DL, LFP, T}(s)
Quinovose(; LFP = PyranoseForm, DL = DForm) = Quinovose{DL, LFP, Nothing}(nothing)
# Qui = 6dGlc
# SQui = 6SQui
"""
    Rhamnose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Rhamnose.
"""
struct Rhamnose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Rhamnose(s::T; LFP = PyranoseForm, DL = LForm) where T = Rhamnose{DL, LFP, T}(s)
Rhamnose(; LFP = PyranoseForm, DL = LForm) = Rhamnose{DL, LFP, Nothing}(nothing)
# Rha = 6dMan(L)
"""
    Fucose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

Fucose.
"""
struct Fucose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Fucose(s::T; LFP = PyranoseForm, DL = LForm) where T = Fucose{DL, LFP, T}(s)
Fucose(; LFP = PyranoseForm, DL = LForm) = Fucose{DL, LFP, Nothing}(nothing)
# Fuc = 6dGal(L)
"""
    Sixdeoxygulose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxygulose.
"""
struct Sixdeoxygulose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Sixdeoxygulose(s::T; LFP = PyranoseForm, DL = DForm) where T = Sixdeoxygulose{DL, LFP, T}(s)
Sixdeoxygulose(; LFP = PyranoseForm, DL = DForm) = Sixdeoxygulose{DL, LFP, Nothing}(nothing)
# 6dGul (no NAc)
"""
    Sixdeoxyaltose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxyaltose.
"""
struct Sixdeoxyaltose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Sixdeoxyaltose(s::T; LFP = PyranoseForm, DL = LForm) where T = Sixdeoxyaltose{DL, LFP, T}(s)
Sixdeoxyaltose(; LFP = PyranoseForm, DL = LForm) = Sixdeoxyaltose{DL, LFP, Nothing}(nothing)
"""
    Sixdeoxytalose{D, P, T} <: AbstractDeoxyhexose{D, P, T}

6-deoxytalose.
"""
struct Sixdeoxytalose{D, P, T} <: AbstractDeoxyhexose{D, P, T}
    substituent::T
end
Sixdeoxytalose(s::T; LFP = PyranoseForm, DL = DForm) where T = Sixdeoxytalose{DL, LFP, T}(s)
Sixdeoxytalose(; LFP = PyranoseForm, DL = DForm) = Sixdeoxytalose{DL, LFP, Nothing}(nothing)
