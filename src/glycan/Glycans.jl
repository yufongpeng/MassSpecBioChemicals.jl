module Glycans

using Reexport, MLStyle
using ..MassSpecBioChemicals
@reexport using ..MassSpecBioChemicals.BasicCompounds, ..MassSpecBioChemicals.Metabolites
using MassSpecChemicals: AbstractChemical
import MassSpecChemicals: parse_chemical, getchemicalattr
using ..MassSpecBioChemicals.Proteins: parse_aa_fg
using ..MassSpecBioChemicals: lk
using IterTools
import Base: show
import ..MassSpecBioChemicals: repr_linkage, makelinkage, transformlinkage, getchaincomponent, getchainlinkage, ischainedchemical, dehydroxyposition, dehydrogenposition
export Saccharide, Monosaccharide, 
        AbstractHexose, Hexose, Glucose, Galactose, Mannose, 
        AbstractDeoxyhexose, Deoxyhexose, Fucose, 
        AbstractHexosamine, Hexosamine, Glucosamine, Galactosamine, Mannosamine, 
        AbstractNacetylhexosamine, Nacetylhexosamine, Nacetylgluosamine, Nacetylgalactosamine, Nacetylmannosamine, 
        AbstractHexuronicAcid, HexuronicAcid, GlucuronicAcid, GalacturonicAcid, MannuronicAcid,
        SialicAcid, NeuraminicAcid, NacetylneuraminicAcid, NglycolylneuraminicAcid, Kdn,
        AbstractPentose, Pentose, Ribose, Arabinose, Xylose, AbstractDeoxypentose, Deoxypentose, Deoxyribose, 
        Inositol, Sulfoquinovose, 
        AbstractGlycan, Glycan, GlyComp,
        AbstractAnomerposition, Anomerposition, Alphaposition, Betaposition,
        Ganglioseries,
        GM4,
        SM4,
        Lac,
        SM3,
        SM2,
        SM1a,
        SM1b,
        SB1a,
        SB1,
        Ganglio0series,
        GA2,
        GA1,
        GM1b,
        GM1α,
        GD1c,
        GD1α,
        GD1e,
        GanglioAseries,
        GM3,
        GM2,
        GM1a,
        GD1a,
        GD1aα,
        GD1aa,
        GT1a,
        GT1aα,
        GT1aa,
        GanglioBseries,
        GD3,
        GD2,
        GD1b,
        GT1b,
        GT1bα,
        GT1ba,
        GQ1b,
        GQ1bα,
        GQ1ba,
        GanglioCseries,
        GT3,
        GT2,
        GT1c,
        GQ1c,
        GQ1cα,
        GQ1ca,
        GP1c,
        GP1cα,
        GP1ca,
        SM1,
        GM1,
        GD1,
        GT1,
        GQ1,
        GP1,
        Globoseries,
        Gb3,
        Gb4,
        Isogloboseries,
        iGb3,
        iGb4,
        Lactoseries,
        Lc3,
        LM1
        
abstract type Saccharide <: AbstractChemical end
abstract type Monosaccharide{T} <: Saccharide end
abstract type AbstractHexose{T} <: Monosaccharide{T} end
abstract type AbstractDeoxyhexose{T} <: Monosaccharide{T} end
abstract type AbstractHexosamine{T} <: Monosaccharide{T} end
abstract type AbstractNacetylhexosamine{T} <: Monosaccharide{T} end
abstract type AbstractHexuronicAcid{T} <: Monosaccharide{T} end
abstract type SialicAcid{T} <: Monosaccharide{T} end
abstract type AbstractPentose{T} <: Monosaccharide{T} end
abstract type AbstractDeoxypentose{T} <: Monosaccharide{T} end

include("anomer.jl")
#=
T
Nothing
Vector{<: Pair{<: FunctionalGroup, UInt8}}
Vector{<: Pair{UInt8, <: FunctionalGroup}}
=#
"""
    Hexose{T} <: AbstractHexose{T}

Hexose with or without substituents.
"""
struct Hexose{T} <: AbstractHexose{T}
    substituent::T
end
Hexose() = Hexose(nothing)
"""
    Glucose{T} <: AbstractHexose{T}

Glucose with or without substituents.
"""
struct Glucose{T} <: AbstractHexose{T}
    substituent::T
end
Glucose() = Glucose(nothing)
"""
    Galactose{T} <: AbstractHexose{T}

Galactose with or without substituents.
"""
struct Galactose{T} <: AbstractHexose{T}    
    substituent::T
end
Galactose() = Galactose(nothing)
"""
    Mannose{T} <: AbstractHexose{T}

Mannose with or without substituents.
"""
struct Mannose{T} <: AbstractHexose{T}    
    substituent::T
end
Mannose() = Mannose(nothing)

"""
    Hexosamine{T} <: AbstractHexosamine{T}

Hexosamine with or without substituents.
"""
struct Hexosamine{T} <: AbstractHexosamine{T}
    substituent::T
end
Hexosamine() = Hexosamine(nothing)
"""
    Glucosamine{T} <: AbstractHexosamine{T}

Glucosamine with or without substituents.
"""
struct Glucosamine{T} <: AbstractHexosamine{T}
    substituent::T
end
Glucosamine() = Glucosamine(nothing)
"""
    Galactosamine{T} <: AbstractHexosamine{T}

Galactosamine with or without substituents.
"""
struct Galactosamine{T} <: AbstractHexosamine{T}    
    substituent::T
end
Galactosamine() = Galactosamine(nothing)
"""
    Mannosamine{T} <: AbstractHexosamine{T}

Mannosamine with or without substituents.
"""
struct Mannosamine{T} <: AbstractHexosamine{T}    
    substituent::T
end
Mannosamine() = Mannosamine(nothing)

"""
    Nacetylhexosamine{T} <: AbstractNacetylhexosamine{T}

N-acetylhexosamine with or without substituents.
"""
struct Nacetylhexosamine{T} <: AbstractNacetylhexosamine{T}
    substituent::T
end
Nacetylhexosamine() = Nacetylhexosamine(nothing)
"""
    Nacetylgluosamine{T} <: AbstractNacetylhexosamine{T}

N-acetylglucosamine with or without substituents.
"""
struct Nacetylgluosamine{T} <: AbstractNacetylhexosamine{T}
    substituent::T
end
Nacetylgluosamine() = Nacetylgluosamine(nothing)
"""
    Nacetylgalactosamine{T} <: AbstractNacetylhexosamine{T}

N-acetylgalactosamine with or without substituents.
"""
struct Nacetylgalactosamine{T} <: AbstractNacetylhexosamine{T}    
    substituent::T
end
Nacetylgalactosamine() = Nacetylgalactosamine(nothing)
"""
    Nacetylmannosamine{T} <: AbstractNacetylhexosamine{T}

N-acetylmannosamine with or without substituents.
"""
struct Nacetylmannosamine{T} <: AbstractNacetylhexosamine{T}    
    substituent::T
end
Nacetylmannosamine() = Nacetylmannosamine(nothing)

"""
    HexuronicAcid{T} <: AbstractHexuronicAcid{T}

Hexuronate with or without substituents.
"""
struct HexuronicAcid{T} <: AbstractHexuronicAcid{T}
    substituent::T
end
HexuronicAcid() = HexuronicAcid(nothing)
"""
    GlucuronicAcid{T} <: AbstractHexuronicAcid{T}

Glucuronic acid with or without substituents.
"""
struct GlucuronicAcid{T} <: AbstractHexuronicAcid{T}
    substituent::T
end
GlucuronicAcid() = GlucuronicAcid(nothing)
"""
    GalacturonicAcid{T} <: AbstractHexuronicAcid{T}

Galactcuronic acid with or without substituents.
"""
struct GalacturonicAcid{T} <: AbstractHexuronicAcid{T}
    substituent::T
end
GalacturonicAcid() = GalacturonicAcid(nothing)
"""
    MannuronicAcid{T} <: AbstractHexuronicAcid{T}

Mannuronic acid with or without substituents.
"""
struct MannuronicAcid{T} <: AbstractHexuronicAcid{T}
    substituent::T
end
MannuronicAcid() = MannuronicAcid(nothing)

"""
    Deoxyhexose{T} <: AbstractDeoxyhexose{T}

Deoxyhexose with or without substituents.
"""
struct Deoxyhexose{T} <: AbstractDeoxyhexose{T}
    substituent::T
end
Deoxyhexose() = Deoxyhexose(nothing)
"""
    Fucose{T} <: AbstractDeoxyhexose{T}

Fucose with or without substituents.
"""
struct Fucose{T} <: AbstractDeoxyhexose{T}
    substituent::T
end
Fucose() = Fucose(nothing)

"""
    NeuraminicAcid{T} <: SialicAcid{T}

Neuraminic acid with or without substituents.
"""
struct NeuraminicAcid{T} <: SialicAcid{T}
    substituent::T
end
NeuraminicAcid() = NeuraminicAcid(nothing)
"""
    NacetylneuraminicAcid{T} <: SialicAcid{T}

N-acetylneuraminic acid (nana) with or without substituents.
"""
struct NacetylneuraminicAcid{T} <: SialicAcid{T}
    substituent::T
end
NacetylneuraminicAcid() = NacetylneuraminicAcid(nothing)
"""
    NglycolylneuraminicAcid{T} <: SialicAcid{T}

N-glycolylneuraminic acid with or without substituents.
"""
struct NglycolylneuraminicAcid{T} <: SialicAcid{T}
    substituent::T
end
NglycolylneuraminicAcid() = NglycolylneuraminicAcid(nothing)
"""
    Kdn <: SialicAcid

Keto-deoxy-glycero-galactonononic acid with or without substituents.
"""
struct Kdn{T} <: SialicAcid{T}
    substituent::T
end
Kdn() = Kdn(nothing)

"""
    Pentose{T} <: AbstractPentose{T}

Pentose with or without substituents.
"""
struct Pentose{T} <: AbstractPentose{T}
    substituent::T
end
Pentose() = Pentose(nothing)
"""
    Arabinose{T} <: AbstractPentose{T}

Arabinose with or without substituents.
"""
struct Arabinose{T} <: AbstractPentose{T}
    substituent::T
end
Arabinose() = Arabinose(nothing)
"""
    Xylose{T} <: AbstractPentose{T}

Xylose with or without substituents.
"""
struct Xylose{T} <: AbstractPentose{T}
    substituent::T
end
Xylose() = Xylose(nothing)
"""
    Ribose{T} <: AbstractPentose{T}

Ribose with or without substituents.
"""
struct Ribose{T} <: AbstractPentose{T}
    substituent::T
end
Ribose() = Ribose(nothing)

"""
    Deoxypentose{T} <: AbstractDeoxypentose{T}

Deoxypentose with or without substituents.
"""
struct Deoxypentose{T} <: AbstractDeoxypentose{T}
    substituent::T
end
Deoxypentose() = Deoxypentose(nothing)
"""
    Deoxyribose{T} <: AbstractDeoxypentose{T}

Deoxyribose with or without substituents.
"""
struct Deoxyribose{T} <: AbstractDeoxypentose{T}
    substituent::T
end
Deoxyribose() = Deoxyribose(nothing)
"""
    Inositol{T} <: Monosaccharide{T}

Inositol with or without substituents.
"""
struct Inositol{T} <: Monosaccharide{T}
    substituent::T
end
Inositol() = Inositol(nothing)
"""
Sulfoquinovose{T} <: Monosaccharide{T}

Sulfoquinovose with or without substituents.
"""
struct Sulfoquinovose{T} <: Monosaccharide{T}
    substituent::T
end
Sulfoquinovose() = Sulfoquinovose(nothing)

abstract type AbstractGlycan <: Saccharide end
struct GM4 <: AbstractGlycan end
struct SM4 <: AbstractGlycan end
struct Lac <: AbstractGlycan end #?
abstract type Ganglioseries <: AbstractGlycan end
struct SM3 <: Ganglioseries end
struct SM2 <: Ganglioseries end
struct SM1a <: Ganglioseries end
struct SM1b <: Ganglioseries end
struct SB1a <: Ganglioseries end
const SB1 =  SB1a
abstract type Ganglio0series <: Ganglioseries end
struct GA2 <: Ganglio0series end
struct GA1 <: Ganglio0series end
struct GM1b <: Ganglio0series end
struct GM1α <: Ganglio0series end
struct GD1c <: Ganglio0series end
struct GD1α <: Ganglio0series end
const GD1e = GD1α
abstract type GanglioAseries <: Ganglioseries end
struct GM3 <: GanglioAseries end
struct GM2 <: GanglioAseries end
struct GM1a <: GanglioAseries end
struct GD1a <: GanglioAseries end
struct GD1aα <: GanglioAseries end
const GD1aa = GD1aα
struct GT1a <: GanglioAseries end
struct GT1aα <: GanglioAseries end
const GT1aa = GT1aα
abstract type GanglioBseries <: Ganglioseries end
struct GD3 <: GanglioBseries end
struct GD2 <: GanglioBseries end
struct GD1b <: GanglioBseries end
struct GT1b <: GanglioBseries end
struct GT1bα <: GanglioBseries end
const GT1ba = GT1bα
struct GQ1b <: GanglioBseries end
struct GQ1bα <: GanglioBseries end
const GQ1ba = GQ1bα
abstract type GanglioCseries <: Ganglioseries end
struct GT3 <: GanglioCseries end
struct GT2 <: GanglioCseries end
struct GT1c <: GanglioCseries end
struct GQ1c <: GanglioCseries end
struct GQ1cα <: GanglioCseries end
const GQ1ca = GQ1cα
struct GP1c <: GanglioCseries end
struct GP1cα <: GanglioCseries end
const GP1ca = GP1cα

struct SM1{T} <: Ganglioseries
    isomer::T
end # SM1a, SM1b, xxxxx

struct GM1{T} <: Ganglioseries
    isomer::T
end # GM1a, x, x, GM1b, x, x, x, GM1α

struct GD1{T} <: Ganglioseries
    isomer::T
end # GD1a, GD1b, x, GD1c, GD1aα, x, x, GD1α 
struct GT1{T} <: Ganglioseries
    isomer::T
end # GT1a, GT1b, GT1c, x, GT1aα, GT1bα, x, x
struct GQ1{T} <: Ganglioseries
    isomer::T
end # x, GQ1b, GQ1c, x, x, GQ1bα, GQ1cα, x
struct GP1{T} <: Ganglioseries
    isomer::T
end # x, x GP1c, x, x, x, GP1cα, x

function SM1()
    SM1(SM1a(), SM1b())
end
function GM1()
    GM1(GM1a(), GM1b(), GM1α())
end
function GD1()
    GD1(GD1a(), GD1aα(), GD1b(), GD1c(), GD1α())
end
function GT1()
    GT1(GT1a(), GT1aα(), GT1b(), GT1bα(), GT1c())
end
function GQ1()
    GQ1(GQ1b(), GQ1bα(), GQ1c(), GQ1cα())
end
function GP1()
    GP1(GP1c(), GP1cα())
end
function SM1(isomer...)
    SM1(isomer)
end
function GM1(isomer...)
    GM1(isomer)
end
function GD1(isomer...)
    GD1(isomer)
end
function GT1(isomer...)
    GT1(isomer)
end
function GQ1(isomer...)
    GQ1(isomer)
end
function GP1(isomer...)
    GP1(isomer)
end

abstract type Globoseries <: AbstractGlycan end
struct Gb3 <: Globoseries end
struct Gb4 <: Globoseries end
abstract type Isogloboseries <: AbstractGlycan end
struct iGb3 <: Isogloboseries end
struct iGb4 <: Isogloboseries end
abstract type Lactoseries <: AbstractGlycan end
struct Lc3 <: Lactoseries end
struct LM1 <: Lactoseries end

"""
    Glycan{T} <: Saccharide

Oligosaccarides or Polysaccharides with or without defined order.
"""
struct Glycan{T} <: AbstractGlycan
    chain::T
    linkage::Vector{Pair{AbstractAnomerposition, Linkageposition}}
end
Glycan(sugar::Vararg{<: Saccharide}) = Glycan(sugar, Pair{AbstractAnomerposition, UInt8}[])
#= Tuple: order, length of linkage == chain, defined reducing end anomer, 
S::UInt8

S::Pair{UInt8, UInt8} x, y => p, a = divrem(x, 3) ap-y, ay
a = 0, unknown anomer; a = 1, α; a = 2, β
p = 0, default reducing end 1/2
=#

struct GlyComp{T} <: AbstractGlycan
    comp::Vector{Pair{T, UInt8}}
end
# Vector{Pair{Monosaccharide, UInt8}}: type => number, no order

include("interface.jl")
include("io.jl")
end