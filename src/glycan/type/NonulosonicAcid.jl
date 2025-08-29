"""
    NeuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

Neuraminic acid.
"""
struct NeuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
end
NeuraminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = NeuraminicAcid{DL, LFP, T}(s)
NeuraminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = NeuraminicAcid{DL, LFP, Nothing}(nothing)
"""
    NacetylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

N-acetylneuraminic acid (nana).
"""
struct NacetylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
end
NacetylneuraminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = NacetylneuraminicAcid{DL, LFP, T}(s)
NacetylneuraminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = NacetylneuraminicAcid{DL, LFP, Nothing}(nothing)
"""
    NglycolylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}

N-glycolylneuraminic acid.
"""
struct NglycolylneuraminicAcid{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
end
NglycolylneuraminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = NglycolylneuraminicAcid{DL, LFP, T}(s)
NglycolylneuraminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = NglycolylneuraminicAcid{DL, LFP, Nothing}(nothing)
"""
    Kdn{D, P, T} <: SialicAcid{D, P, T}

2-Keto-3-deoxy-nononic acid. 
"""
struct Kdn{D, P, T} <: SialicAcid{D, P, T}
    substituent::T
end
Kdn(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = Kdn{DL, LFP, T}(s)
Kdn(; LFP = PyranoseForm, DL = NoDLForm) = Kdn{DL, LFP, Nothing}(nothing)
# Kdn, Neu = Kdn5N

# 9dKdn5,7N: Pse, Leg, Aci, 4eLeg, 8eLeg, 8eAci

"""
   LegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

Legionaminic acid.
"""
struct LegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
LegionaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = LegionaminicAcid{DL, LFP, T}(s)
LegionaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = LegionaminicAcid{DL, LFP, Nothing}(nothing)
# Leg = 9dKdn5,7N
"""
   FourepilegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

4-epi-legionaminic acid.
"""
struct FourepilegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
FourepilegionaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = FourepilegionaminicAcid{DL, LFP, T}(s)
FourepilegionaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = FourepilegionaminicAcid{DL, LFP, Nothing}(nothing)
"""
   EightepilegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

8-epi-legionaminic acid.
"""
struct EightepilegionaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
EightepilegionaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = EightepilegionaminicAcid{DL, LFP, T}(s)
EightepilegionaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = EightepilegionaminicAcid{DL, LFP, Nothing}(nothing)
"""
   AcinetaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

Acinetaminic acid.
"""
struct AcinetaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
AcinetaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = AcinetaminicAcid{DL, LFP, T}(s)
AcinetaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = AcinetaminicAcid{DL, LFP, Nothing}(nothing)
# Aci = 9d,7,8eKdn5,7N
"""
   EightepiacinetaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

8-epei-acinetaminic acid.
"""
struct EightepiacinetaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
EightepiacinetaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = EightepiacinetaminicAcid{DL, LFP, T}(s)
EightepiacinetaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = EightepiacinetaminicAcid{DL, LFP, Nothing}(nothing)
"""
   PseudaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}

Pseudaminic acid.
"""
struct PseudaminicAcid{D, P, T} <: DideoxynonulosonicAcid{D, P, T}
    substituent::T
end
PseudaminicAcid(s::T; LFP = PyranoseForm, DL = NoDLForm) where T = PseudaminicAcid{DL, LFP, T}(s)
PseudaminicAcid(; LFP = PyranoseForm, DL = NoDLForm) = PseudaminicAcid{DL, LFP, Nothing}(nothing)
# Pse = 9d,5,7,8eKdn5,7N
