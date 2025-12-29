module Glycans

using Reexport, MLStyle
using ..MassSpecBioChemicals
@reexport using ..MassSpecBioChemicals.BasicCompounds, ..MassSpecBioChemicals.Metabolites
using MassSpecChemicals: AbstractChemical
import MassSpecChemicals: parse_chemical, getchemicalattr
using ..MassSpecBioChemicals.Proteins: parse_aa_fg
using ..MassSpecBioChemicals: lk, AbstractFunctionalGroup, ntotallinkage, leavinggroupelements
using IterTools
import Base: show, length
import ..MassSpecBioChemicals: repr_linkage, makelinkage, transformlinkage, chainedchemical, getchaincomponent, getchainlinkage, getchainconfig, composition, ischainedchemical, requirelinkage, requireconfig, dehydroxyposition, dehydrogenposition, GlyceraldehydeSystem, losswaterelements, deisomerize, decompose
export Saccharide, Monosaccharide, 
        AbstractHexose, Hexose, Glucose, Mannose, Galactose, Gulose, Altose, Allose, Talose, Idose, Apiose, Fructose, Tagatose, Sorbose, Psicose, 
        AbstractHexosamine, Hexosamine, Glucosamine, Mannosamine, Galactosamine, Gulosamine, Altosamine, Allosamine, Talosamine, Idosamine, 
        AbstractNacetylhexosamine, Nacetylhexosamine, Nacetylgluosamine, Nacetylmannosamine, Nacetylgalactosamine, Nacetygalactosamine, Nacetygulosamine, Nacetyaltosamine, Nacetyallosamine, Nacetytalosamine, Nacetyidosamine, 
        AbstractDeoxyhexose, Deoxyhexose, Quinovose, Rhamnose, Fucose, Sixdeoxygulose, Sixdeoxyaltose, Sixdeoxytalose, 
        AbstractNacetyldeoxyhexosamine, Nacetyldeoxyhexosamine, Nacetylquinovosamine, Nacetylrhamnosamine, Nacetylfucosamine, Nacetylsixdeoxyaltosamine, Nacetylsixdeoxytalosamine, 
        AbstractDideoxyhexose, Dideoxyhexose, 
        AbstractHexuronicAcid, HexuronicAcid, GlucuronicAcid, MannuronicAcid, GalacturonicAcid, GuluronicAcid, AlturonicAcid, AlluronicAcid, TaluronicAcid, IduronicAcid, 
        NonulosonicAcid, SialicAcid, NeuraminicAcid, NacetylneuraminicAcid, NglycolylneuraminicAcid, Kdn,
        AbstractdiaminodideoxynonulosonicAcid, DiaminodideoxynonulosonicAcid, LegionaminicAcid, FourepilegionaminicAcid, EightepilegionaminicAcid, AcinetaminicAcid, EightepiacinetaminicAcid, PseudaminicAcid, 

        AbstractPentose, Pentose, Arabinose, Lyxose, Xylose, Ribose, AbstractDeoxypentose, Deoxypentose, Deoxyribose, 
        Inositol, 

        AbstractGlycan, SeriesGlycan, Glycan, GlyComp,
        AbstractAnomerposition, Anomerposition, Alphaposition, Betaposition,
        LinearForm, 
        PyranoseForm, 
        FuranoseForm, 
        CycliceForm, 
        Deoxy, 
        Epi, 
        
        SimpleGlcseries,
        GM4,
        SM4,
        Lac,
        Ganglioseries,
        Sulfoganglioseries,
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
        Gb5,
        CytolipinK,
        SSEA3Antigen,
        SSEA4Antigen,
        SSEA4α6Antigen,
        SSEA4a6Antigen,
        TypeIVHAntigen,
        TypeIVAAntigen,
        TypeIVBAntigen,
        ForssmanAntigen,
        ParaForssmanAntigen,
        BranchedForssmanAntigen,
        GloboLex9,

        Isogloboseries,
        iGb3,
        iGb4,
        iGb5,
        ForssmanlikeiGb4,
        CytolipinR,

        Lactoseries,
        Lc3,
        Lc4,
        Lea,
        Leb,
        TypeIHAntigen,
        TypeIAAntigen,
        Aleb,
        TypeIBAntigen,
        Bleb,
        LexA,
        LeyA,
        LexA9,
        LeyA9,
        LM1,
        sLea,
        dsLea,

        Neolactoseries,
        nLc4,
        nLc5,
        nLc6,
        IAntigen,
        TypeIIH7Antigen,
        TypeIIA8Antigen,
        TypeIIB8Antigen,
        nLM1,
        TypeIIHAntigen,
        TypeIIAAntigen,
        TypeIIBAntigen,
        iAntigen,
        Lex,
        Ley,
        Lex5,
        Ley6,
        SSEA1Antigen,
        Lex7,
        LexX,
        LexX8,
        DimericLex8,
        Ley8,
        LeyX,
        LeyX9,
        TrimericLey9,
        Lex9,
        LexX10,
        DimericLex10,
        LexXX,
        LexXX11,
        TrimericLex11,
        LeaX,
        LebX,
        sLex,
        sLex6,
        sLeaX,
        PAntigen

        
abstract type Saccharide <: AbstractChemical end
"""
"""
abstract type Monosaccharide{D, P, T} <: Saccharide end
abstract type AbstractHexose{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractHexosamine{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractNacetylhexosamine{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractHexuronicAcid{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractDeoxyhexose{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractNacetyldeoxyhexosamine{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractDideoxyhexose{D, P, T} <: Monosaccharide{D, P, T} end
abstract type NonulosonicAcid{D, P, T} <: Monosaccharide{D, P, T} end
abstract type SialicAcid{D, P, T} <: NonulosonicAcid{D, P, T} end
# abstract type DideoxynonulosonicAcid{D, P, T} <: NonulosonicAcid{D, P, T} end
abstract type AbstractdiaminodideoxynonulosonicAcid{D, P, T} <: NonulosonicAcid{D, P, T} end
abstract type AbstractPentose{D, P, T} <: Monosaccharide{D, P, T} end
abstract type AbstractDeoxypentose{D, P, T} <: Monosaccharide{D, P, T} end

include("anomer.jl")

abstract type MonosaccharideForm end 
struct LinearForm <: MonosaccharideForm end 
struct PyranoseForm <: MonosaccharideForm end 
struct FuranoseForm <: MonosaccharideForm end 
struct CycliceForm <: MonosaccharideForm end 

struct Deoxy <: AbstractFunctionalGroup end
struct Epi <: AbstractFunctionalGroup end

include(joinpath("type", "Hex.jl"))
include(joinpath("type", "HexNAc.jl"))
include(joinpath("type", "HexN.jl"))
include(joinpath("type", "HexA.jl"))
include(joinpath("type", "dHex.jl"))
include(joinpath("type", "dHexNAc.jl"))
include(joinpath("type", "ddHex.jl"))
include(joinpath("type", "NonulosonicAcid.jl"))
include(joinpath("type", "Pen.jl"))
#=
T
Nothing
Vector{<: Pair{<: FunctionalGroup, UInt8}}
Vector{<: Pair{UInt8, <: FunctionalGroup}}
=#

# Unknown 
# LDManHep, DDManHep = 6eLDManHep
# Kdo = 3d5eDDManHep1COOH or 3dDDGalHep1COOH, octanulosonic acid
# Dha: 3dGalA1COOH

# Bac: 2,4N,6dGlc
# Mur: 2N,3carboxyethylGlc

"""
    Inositol{D, P, T} <: Monosaccharide{D, P, T}

Inositol.
"""
struct Inositol{D, P, T} <: Monosaccharide{D, P, T}
    substituent::T
end
Inositol() = Inositol{Nothing, CycliceForm, Nothing}(nothing)
Inositol(x::T) where T = Inositol{Nothing, CycliceForm, T}(x)

abstract type AbstractGlycan <: Saccharide end
struct SeriesGlycan{G, T} <: AbstractGlycan 
    type::G
    substituent::T
end
# FG => # 
# MS => #
# p => MS
# GD1a(3[3]Neu5,?Ac/6[4]NeuGc)
# GT1?b/ba(Ac2)
# GD1b((NeuAc2)2)

abstract type GlycanType end
abstract type SimpleGlcseries <: GlycanType end
struct GM4 <: SimpleGlcseries end
struct SM4 <: SimpleGlcseries end
struct Lac <: SimpleGlcseries end #

abstract type Ganglioseries <: GlycanType end
abstract type Sulfoganglioseries <: Ganglioseries end
struct SM3 <: Sulfoganglioseries end
struct SM2 <: Sulfoganglioseries end
struct SM1a <: Sulfoganglioseries end
struct SM1b <: Sulfoganglioseries end
struct SB1a <: Sulfoganglioseries end
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

# GalNAcα1-4
# Galα1-3
# O-acetyl 
# NeuGc/Kdn
# Fucα1-2

abstract type Globoseries <: GlycanType end
struct Gb3 <: Globoseries end
struct Gb4 <: Globoseries end
struct Gb5 <: Globoseries end
const CytolipinK = Gb4
const SSEA3Antigen = Gb5
struct SSEA4Antigen <: Globoseries end
struct SSEA4α6Antigen <: Globoseries end
const SSEA4a6Antigen = SSEA4α6Antigen
struct TypeIVHAntigen <: Globoseries end
struct TypeIVAAntigen <: Globoseries end
struct TypeIVBAntigen <: Globoseries end
struct ForssmanAntigen <: Globoseries end
struct ParaForssmanAntigen <: Globoseries end
struct BranchedForssmanAntigen <: Globoseries end
struct GloboLex9 <: Globoseries end

abstract type Isogloboseries <: GlycanType end
struct iGb3 <: Isogloboseries end
struct iGb4 <: Isogloboseries end
struct iGb5 <: Isogloboseries end
struct ForssmanlikeiGb4 <: Isogloboseries end
const CytolipinR = iGb4

abstract type Lactoseries <: GlycanType end
struct Lc3 <: Lactoseries end
struct Lc4 <: Lactoseries end
struct Lea <: Lactoseries end
struct Leb <: Lactoseries end
struct TypeIHAntigen <: Lactoseries end
struct TypeIAAntigen <: Lactoseries end
struct Aleb <: Lactoseries end
struct TypeIBAntigen <: Lactoseries end
struct Bleb <: Lactoseries end
struct LexA <: Lactoseries end
struct LeyA <: Lactoseries end
const LexA9 = LexA
const LeyA9 = LeyA
struct LM1 <: Lactoseries end
struct sLea <: Lactoseries end
struct dsLea <: Lactoseries end

abstract type Neolactoseries <: GlycanType end
struct nLc4 <: Neolactoseries end
struct nLc5 <: Neolactoseries end
struct nLc6 <: Neolactoseries end
struct IAntigen <: Neolactoseries end
struct TypeIIH7Antigen <: Neolactoseries end
struct TypeIIA8Antigen <: Neolactoseries end
struct TypeIIB8Antigen <: Neolactoseries end
struct nLM1 <: Neolactoseries end
struct TypeIIHAntigen <: Neolactoseries end
struct TypeIIAAntigen <: Neolactoseries end
struct TypeIIBAntigen <: Neolactoseries end
const iAntigen = nLc6
struct Lex <: Neolactoseries end
struct Ley <: Neolactoseries end
const Lex5 = Lex
const Ley6 = Ley
const SSEA1Antigen = Lex5
struct Lex7 <: Neolactoseries end
struct LexX <: Neolactoseries end
const LexX8 = LexX
const DimericLex8 = LexX8
struct Ley8 <: Neolactoseries end
struct LeyX <: Neolactoseries end
const LeyX9 = LeyX
const TrimericLey9 = LeyX9
struct Lex9 <: Neolactoseries end
struct LexX10 <: Neolactoseries end
const DimericLex10 = LexX10
struct LexXX <: Neolactoseries end
const LexXX11 = LexXX
const TrimericLex11 = LexXX11
struct LeaX <: Neolactoseries end
struct LebX <: Neolactoseries end
struct sLex <: Neolactoseries end
const sLex6 = sLex
struct sLeaX <: Neolactoseries end
struct PAntigen <: Neolactoseries end


"""
    Glycan{T} <: AbstractGlycan

Oligosaccarides or Polysaccharides with or without defined order.
"""
struct Glycan{T} <: AbstractGlycan
    chain::T
    linkage::Vector{Pair{AbstractAnomerposition, Linkageposition}}
end
Glycan(sugar::Vararg{<: Saccharide}) = Glycan(sugar, Pair{AbstractAnomerposition, Linkageposition}[])
Glycan(sugar::T, linkage::Vector) where T = Glycan(sugar, convert(Vector{Pair{AbstractAnomerposition, Linkageposition}}, linkage))
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
GlyComp(sugar::SeriesGlycan, s = nothing) = GlyComp(generic_glycan(sugar), s)
function GlyComp(sugar::GlyComp, s = nothing)
    d = Dict{Monosaccharide, UInt8}()
    push_monosaccharide!(d, sugar)
    if !isnothing(s)
        for (m, n) in s 
            push_monosaccharide!(d, m, n)
        end
    end
    filter!(x -> last(x) > 0, d)
    GlyComp(collect(pairs(d)))
end
function GlyComp(sugar::Glycan, s = nothing)
    d = Dict{Monosaccharide, UInt8}()
    push_monosaccharide!(d, sugar)
    if !isnothing(s)
        for (m, n) in s 
            push_monosaccharide!(d, m, n)
        end
    end
    filter!(x -> last(x) > 0, d)
    GlyComp(collect(pairs(d)))
end
function push_monosaccharide!(d, s::Monosaccharide, x = 0x01) 
    n = get!(d, s, 0x00)
    d[s] = n + x
end
function push_monosaccharide!(d, s::Glycan) 
    for i in getchaincomponent(s)
        push_monosaccharide!(d, i)
    end
end
function push_monosaccharide!(d, s::GlyComp) 
    for (x, i) in getchaincomponent(s)
        push_monosaccharide!(d, x, i)
    end
end

# D for Abe, Bac, Dha, Kdo, Mur, Par, Tyv; L for Col; DD for Kdn, Neu, Leg, 4eLeg; LL for Pse, Aci; LD for 8eLeg; DL for 8eAci
include("glycan.jl")
include("interface.jl")
include("io.jl")
end