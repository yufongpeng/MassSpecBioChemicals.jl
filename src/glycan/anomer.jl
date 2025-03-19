abstract type AbstractAnomerposition <: AbstractLinkageposition end
struct Anomerposition <: AbstractAnomerposition
    position::UInt8
end
ap(x) = Anomerposition(UInt8(x))
struct Alphaposition <: AbstractAnomerposition
    position::UInt8
end
α(x) = Alphaposition(UInt8(x))
struct Betaposition <: AbstractAnomerposition
    position::UInt8
end
β(x) = Betaposition(UInt8(x))