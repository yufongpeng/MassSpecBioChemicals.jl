function chiralitize(::Type{T}, args, kwargs) where {T <: Monosaccharide}
    if T === Base.typename(T).wrapper
        # keyword constructors resolve `LFP`, `DL` and `substituent` with type-specific defaults
        return T, args, kwargs
    end
    T, args, Base.structdiff(NamedTuple(kwargs), NamedTuple{(:LFP, :DL)})
end
chiralitize(::Type{T}, args, kwargs) where {T <: Union{Glycan, GlyComp, SeriesGlycan}} = T, args, kwargs

function losswaterelements(sugar::Monosaccharide, p, q)
    e = copy(chemicalelements(sugar))
    if p != 1 
        i = findfirst(x -> first(x) == "H", e)
        e[i] = "H" => (last(e[i]) - 1)
    end
    if q != 0
        i = findfirst(x -> first(x) == "O", e)
        e[i] = "O" => (last(e[i]) - 1)
        i = findfirst(x -> first(x) == "H", e)
        e[i] = "H" => (last(e[i]) - 1)
    end
    filter!(x -> last(x) != 0, e)
end

# MassSpecBioChemicals
ischainedchemical(::Glycan) = true
ischainedchemical(::Type{<: Glycan}) = true
requirelinkage(::Type{<: Glycan}) = true
requireconfig(::Type{<: Glycan}) = false
ischainedchemical(::GlyComp) = true # ?
ischainedchemical(::Type{<: GlyComp}) = true # ?
requirelinkage(::Type{<: GlyComp}) = false
requireconfig(::Type{<: GlyComp}) = false
getchaincomponent(sugar::Glycan) = sugar.chain
getchainlinkage(sugar::Glycan) = sugar.linkage
getchainconfig(sugar::Glycan) = missing 
getchaincomponent(sugar::GlyComp) = sugar.comp
getchainlinkage(sugar::GlyComp) = missing 
getchainconfig(sugar::GlyComp) = missing 
chainedchemical(::Type{<: Glycan}, chemicals; linkage = Vector{Pair{AbstractAnomerposition, Linkageposition}}[], kwargs...) = Glycan(chemicals, linkage)
chainedchemical(::Type{<: GlyComp}, chemicals; kwargs...) = GlyComp(chemicals)

dehydroxyposition(sugar::Saccharide) = 0x01
dehydrogenposition(sugar::Saccharide) = 0x00

repr_linkage(l::Anomerposition) = l.position > 0 ? string(Int(l.position)) : ""
repr_linkage(l::Alphaposition) = l.position > 0 ? string("α", Int(l.position)) : "α"
repr_linkage(l::Betaposition) = l.position > 0 ? string("β", Int(l.position)) : "β"

makelinkage(::Type{Glycan}, a, b) = lk(dehydroxyposition(a)) => lk(dehydrogenposition(b))
function transformlinkage(::Type{Glycan}, m::ChainedChemical)
    Tuple(first(first(ls)) => first(last(ls)) for ls in getchainlinkage(m))
end
function transformlinkage(::Type{ChainedChemical}, m::Glycan)
    Tuple((first(ls), Dehydroxy()) => (last(ls), Dehydrogen()) for ls in getchainlinkage(m))
end
# dissociate_group
getchaincomponent(sugar::SeriesGlycan) = getchaincomponent(generic_glycan(sugar))
composition(sugar::AbstractGlycan) = Dict(getchaincomponent(GlyComp(sugar))...)
composition(sugar::GlyComp) = Dict(getchaincomponent(sugar)...)
composition(mono::Monosaccharide) = Dict(mono => 0x01)
# concatchemical, makechemical for Glycomp
function deisomerize(sugar::AbstractGlycan) 
    comp = composition(sugar)
    d = Dict{Monosaccharide, UInt8}()
    for (k, v) in comp 
        kn = deisomerize(k)
        get!(d, kn, 0x00)
        d[kn] += v 
    end
    GlyComp(collect(pairs(d)))
end
deisomerize(sugar::AbstractHexose{D, P, T}) where {D, P, T} = Hexose{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractPentose{D, P, T}) where {D, P, T} = Pentose{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractDeoxypentose{D, P, T}) where {D, P, T} = Deoxypentose{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractHexuronicAcid{D, P, T}) where {D, P, T} = HexuronicAcid{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractHexosamine{D, P, T}) where {D, P, T} = Hexosamine{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractNacetylhexosamine{D, P, T}) where {D, P, T} = Nacetylhexosamine{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractDeoxyhexose{D, P, T}) where {D, P, T} = Deoxyhexose{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractDideoxyhexose{D, P, T}) where {D, P, T} = Dideoxyhexose{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::AbstractNacetyldeoxyhexosamine{D, P, T}) where {D, P, T} = Nacetyldeoxyhexosamine{DLForm, P, T}(sugar.substituent)
deisomerize(sugar::SialicAcid{D, P, T}) where {D, P, T} = sugar 
deisomerize(sugar::AbstractdiaminodideoxynonulosonicAcid{D, P, T}) where {D, P, T} = DiaminodideoxynonulosonicAcid{D, P, T}(sugar.substituent)
function decompose(sugar::AbstractHexose)
    d = Dict(Hydroxy() => 0x05, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractPentose)
    d = Dict(Hydroxy() => 0x04, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractDeoxypentose)
    d = Dict(Hydroxy() => 0x03, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractHexuronicAcid)
    d = Dict(Hydroxy() => 0x04, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractHexosamine)
    d = Dict(Hydroxy() => 0x04, Amine() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractNacetylhexosamine)
    d = Dict(Hydroxy() => 0x04, Nacetyl() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractDeoxyhexose)
    d = Dict(Hydroxy() => 0x04, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractDideoxyhexose)
    d = Dict(Hydroxy() => 0x04, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractNacetyldeoxyhexosamine)
    d = Dict(Hydroxy() => 0x03, Nacetyl() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::NeuraminicAcid)
    d = Dict(Hydroxy() => 0x05, Amine() => 0x01, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::NacetylneuraminicAcid)
    d = Dict(Hydroxy() => 0x05, Nacetyl() => 0x01, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::NglycolylneuraminicAcid)
    d = Dict(Hydroxy() => 0x05, XLinkedFunctionalGroup(NLinkage(), Glycolyl()) => 0x01, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::Kdn)
    d = Dict(Hydroxy() => 0x06, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end
function decompose(sugar::AbstractdiaminodideoxynonulosonicAcid)
    d = Dict(Hydroxy() => 0x03, Amine() => 0x02, CarboxylicAcidGroup() => 0x01, OLinkage() => 0x01)
    _decompose!(d, sugar)
end

function _decompose!(d, sugar)
    (isnothing(sugar.substituent) || isempty(sugar.substituent)) && return d 
    if first(sugar.substituent) isa Pair{UInt8}
        for (p, k) in sugar.substituent
            fg, lk = snfg_linkage(sugar, p)
            sub = decompose(snfg_add_linkage(k, lk))
            d[fg] -= 0x01
            for (kn, vn) in sub
                get!(d, kn, 0x00)
                d[kn] += vn
            end
        end
    else
        for (k, v) in sugar.substituent
            sub = decompose(snfg_add_linkage(k))
            d[Hydroxy()] -= v
            for (kn, vn) in sub
                get!(d, kn, 0x00)
                d[kn] += vn * v
            end
        end
    end
    d 
end
function decompose(sugar::AbstractGlycan) 
    comp = composition(sugar)
    d = Dict{AbstractChemical, UInt8}(Hydroxy() => 0x02, OLinkage() => 0x00)
    for (k, v) in comp 
        sub = decompose(k)
        for (kn, vn) in sub
            get!(d, kn, 0x00)
            d[kn] += v * vn
        end
        d[Hydroxy()] -= v * 0x02
        d[OLinkage()] += v
    end
    d[OLinkage()] -= 0x01
    d
end

snfg_add_linkage(x, l = OLinkage()) = XLinkedFunctionalGroup(x, l)
snfg_add_linkage(x::Amino, l = OLinkage()) = x
snfg_add_linkage(x::Tauryl, l = OLinkage()) = x
snfg_add_linkage(x::Nformyl, l = OLinkage()) = x
snfg_add_linkage(x::Nacetyl, l = OLinkage()) = x
snfg_add_linkage(x::Oformyl, l = OLinkage()) = x
snfg_add_linkage(x::Oacetyl, l = OLinkage()) = x
snfg_add_linkage(x::Methyl, l = OLinkage()) = l == OLinkage() ? Methoxy() : XLinkedFunctionalGroup(x, l)
snfg_add_linkage(x::Ethyl, l = OLinkage()) = l == OLinkage() ? Ethoxy() : XLinkedFunctionalGroup(x, l)
snfg_add_linkage(x::Formyl, l = OLinkage()) = l == OLinkage() ? Oformyl() : l == NLinkage() ? Nformyl() : XLinkedFunctionalGroup(x, l)
snfg_add_linkage(x::Acetyl, l = OLinkage()) = l == OLinkage() ? Oacetyl() : l == NLinkage() ? Nacetyl() : XLinkedFunctionalGroup(x, l)
snfg_linkage(sugar, p) = Hydroxy(), OLinkage()
snfg_linkage(sugar::AbstractHexosamine, p) = p == 0x02 ? (Amine(), NLinkage()) : (Hydroxy(), OLinkage())
snfg_linkage(sugar::NeuraminicAcid, p) = p == 0x05 ? (Amine(), NLinkage()) : (Hydroxy(), OLinkage())
snfg_linkage(sugar::AbstractdiaminodideoxynonulosonicAcid, p) = (p == 0x05 || p == 0x07) ? (Amine(), NLinkage()) : (Hydroxy(), OLinkage())

implicit_config(x) = false
implicit_config(::Abequose) = true
implicit_config(::Paratose) = true
implicit_config(::Tyvelose) = true
implicit_config(::Colitose) = true
implicit_config(::SialicAcid) = true
implicit_config(::AbstractdiaminodideoxynonulosonicAcid) = true
implicit_config(::SimpleGlcseries) = true
implicit_config(::Ganglioseries) = true
implicit_config(::Globoseries) = true
implicit_config(::Isogloboseries) = true
implicit_config(::Lactoseries) = true
implicit_config(::Neolactoseries) = true