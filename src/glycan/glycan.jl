# Plain data of `Glycan` used by defaults, so that no instance is needed to know the composition
struct GlycanSpec{C <: Tuple, L <: Vector}
    chain::C
    linkage::L
end
getchaincomponent(spec::GlycanSpec) = spec.chain
getchainlinkage(spec::GlycanSpec) = spec.linkage
_build_glycan(spec::GlycanSpec) = Glycan(map(_build_glycan, spec.chain), spec.linkage)
_build_glycan(x) = x
_glycanchain(sugar::Glycan) = collect(Saccharide, getchaincomponent(sugar))
_glycanchain(sugar::GlycanSpec) = collect(Any, getchaincomponent(sugar))
_rebuildglycan(::Glycan, chain, linkage) = Glycan((chain..., ), linkage)
_rebuildglycan(::GlycanSpec, chain, linkage) = GlycanSpec((chain..., ), linkage)
const GlycanOrSpec = Union{Glycan, GlycanSpec}

function _push_composition!(d, spec::GlycanSpec)
    for x in spec.chain
        _push_composition!(d, x)
    end
    d
end
function _push_composition!(d, m::Monosaccharide)
    d[m] = get(d, m, 0x00) + 0x01
    d
end
# composition of generic glycan without constructing `Glycan` or `GlyComp`
function _generic_composition(type::GlycanType)
    d = Dict{Monosaccharide, UInt8}()
    _push_composition!(d, _glycanspec(type))
    collect(pairs(d))
end

generic_glycan(type::GlycanType) = _build_glycan(_glycanspec(type))
generic_glycan(sugar::Glycan) = sugar
generic_glycan(sugar::SeriesGlycan) = _generic_glycan(sugar.type, sugar.substituent)
_generic_glycan(type, ::Nothing) = generic_glycan(type) 
function _generic_glycan(type, substituent::Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}})
    GlyComp(vcat(convert(Vector{Pair{Any, UInt8}}, getchaincomponent(GlyComp(generic_glycan(type)))), substituent))
end
function _generic_glycan(type, substituent::Vector{<: Pair{<: Monosaccharide, UInt8}})
    g = generic_glycan(type)
    if g isa Glycan && any(x -> first(x) isa SialicAcid, substituent)
        i = findall(x -> first(x) isa SialicAcid, substituent)
        sa = first(first(substituent[i]))
        if allequal(first, substituent[i]) && nsialicacid(g) == sum(last, substituent[i])
            replace_allsa!(g, sa)
        else 
            replace_sa_add_mono!(GlyComp(g), substituent)
        end
    else
        GlyComp(g, substituent)
    end
end
nsialicacid(sugar::Glycan) = sum(nsialicacid, getchaincomponent(sugar))
nsialicacid(sugar::Monosaccharide) = 0
nsialicacid(sugar::SialicAcid) = 1

function replace_allsa!(sugar::Glycan, sa)
    chain = collect(Saccharide, getchaincomponent(sugar))
    for (i, mono) in enumerate(chain)
        chain[i] = replace_allsa!(mono, sa)
    end
    Glycan((chain..., ), getchainlinkage(sugar))
end

replace_allsa!(sugar::SialicAcid, sa) = sa
replace_allsa!(sugar::Monosaccharide, sa) = sugar
function replace_sa_add_mono!(sugar::GlyComp, sas)
    d = Dict(getchaincomponent(sugar)...)
    for (sa, n) in sas 
        get!(d, sa, 0x00)
        if sa isa SialicAcid
            d[NacetylneuraminicAcid()] -= n
        end
        d[sa] += n
    end
    filter!(x -> last(x) > 0, d)
    GlyComp(collect(pairs(d)))
end

function _generic_glycan(type, substituent::Vector{<: Pair{<: Tuple{UInt8, <: Vector}}})
    g = generic_glycan(type) 
    s = [[deepcopy(collect(x)), y] for (x, y) in substituent]
    sort!(s; by = first ∘ first)
    replace_monosaccharide!(g, s)
end
function replace_monosaccharide!(sugar::GlycanOrSpec, s)
    chain = reverse!(_glycanchain(sugar))
    link = getchainlinkage(sugar)
    if length(link) == length(chain)
        link = link[begin:end - 1]
    end
    link = reverse(link)
    id = Int[]
    for (i, mono) in enumerate(chain)
        if !isempty(id)
            if mono isa Monosaccharide
                push!(id, i - 1)
                splits = [empty(s) for _ in id]
                slink = link[id]
                for x in s 
                    t = popfirst!(last(first(x)))
                    j = findfirst(x -> last(x).position == t, slink)
                    isnothing(j) && throw(ArgumentError("Cannot find specified glycan subchain `$(first(first(x)))[$(join([t, last(first(x))...], ","))]$(chemicalabbr(last(x)))` from $(repr_linkage.(slink))."))
                    push!(splits[j], deepcopy(x))
                end
                for (i, ss) in zip(id, splits)
                    chain[i + 1] = replace_monosaccharide!(chain[i + 1], ss)
                end
                s = last(splits)
                id = Int[]
            else
                push!(id, i - 1)
            end
            continue
        elseif mono isa GlycanOrSpec 
            push!(id, i - 1)
            continue
        else
            chain[i] = replace_monosaccharide!(chain[i], s)
        end
    end
    _rebuildglycan(sugar, reverse!(chain), getchainlinkage(sugar))
end

function replace_monosaccharide!(sugar::Monosaccharide, s)
    if isempty(s)
        return sugar 
    else
        for (p, v) in s 
            p[begin] -= 1
            # duplication
            if p[begin] == 0
                sugar = v
            end
        end
        filter!(x -> first(first(x)) > 0, s)
    end
    sugar
end

_glycanspec(::GM4) = GlycanSpec((NacetylneuraminicAcid(), Galactose()), [α(2) => lk(3)])
_glycanspec(::SM4) = GlycanSpec((Galactose([0x03 => Sulfo()]), ), [])
_glycanspec(::Lac) = GlycanSpec((Galactose(), Glucose()), [β(1) => lk(4)])
_glycanspec(::SM3) = GlycanSpec((Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(4)])
_glycanspec(::SM2) = GlycanSpec((Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::SM1a) = GlycanSpec((Galactose(), Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::SM1b) = GlycanSpec((Galactose([0x03 => Sulfo()]), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::SB1a) = GlycanSpec((Galactose([0x03 => Sulfo()]), Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GA2) = GlycanSpec((Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GA1) = GlycanSpec((Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GM1b) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GM1α) = GlycanSpec((Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GD1c) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GD1α) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((NacetylneuraminicAcid(),), [α(2) => lk(6)]), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GM3) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GM2) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GM1a) = GlycanSpec((Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GD1a) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GD1aα) = GlycanSpec((Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT1a) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT1aα) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GD3) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GD2) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GD1b) = GlycanSpec((Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT1b) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT1bα) = GlycanSpec((Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GQ1b) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GQ1bα) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT3) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT2) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GT1c) = GlycanSpec((Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GQ1c) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GQ1cα) = GlycanSpec((Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GP1c) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
_glycanspec(::GP1cα) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])

_generic_composition(::SM1) = ([Nacetylgalactosamine() => 0x01, Galactose([0x03 => Sulfo()]) => 0x01, Galactose() => 0x01, Glucose() => 0x01])
_generic_composition(::GM1) = ([NacetylneuraminicAcid() => 0x01, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
_generic_composition(::GD1) = ([NacetylneuraminicAcid() => 0x02, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
_generic_composition(::GT1) = ([NacetylneuraminicAcid() => 0x03, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
_generic_composition(::GQ1) = ([NacetylneuraminicAcid() => 0x04, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
_generic_composition(::GP1) = ([NacetylneuraminicAcid() => 0x05, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])

generic_glycan(type::Union{SM1, GM1, GD1, GT1, GQ1, GP1}) = GlyComp(_generic_composition(type))

_glycanspec(::Gb3) = GlycanSpec((Galactose(), Galactose(), Glucose()), [α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::Gb4) = GlycanSpec((Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::Gb5) = GlycanSpec((Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::ForssmanAntigen) = GlycanSpec((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::ParaForssmanAntigen) = GlycanSpec((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::BranchedForssmanAntigen) = GlycanSpec((Nacetylgalactosamine(), Nacetylgalactosamine(), GlycanSpec((Galactose(), Nacetylgalactosamine()), [β(1) => lk(3), β(1) => lk(4)]), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::SSEA4Antigen) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::SSEA4α6Antigen) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(2) => lk(6), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::TypeIVHAntigen) = GlycanSpec((Fucose(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::TypeIVAAntigen) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::TypeIVBAntigen) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
_glycanspec(::GloboLex9) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), GlycanSpec((Galactose(), ), [β(1) => lk(3)]), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(6), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])

_glycanspec(::iGb3) = GlycanSpec((Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(4)])
_glycanspec(::iGb4) = GlycanSpec((Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])
_glycanspec(::iGb5) = GlycanSpec((Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])
_glycanspec(::ForssmanlikeiGb4) = GlycanSpec((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])

_glycanspec(::Lc3) = GlycanSpec((Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Lc4) = GlycanSpec((Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Lea) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIHAntigen) = GlycanSpec((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Leb) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIAAntigen) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Aleb) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIBAntigen) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Bleb) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LexA) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LeyA) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LM1) = GlycanSpec((NacetylneuraminicAcid(), Nacetylgalactosamine(), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::sLea) = GlycanSpec((NacetylneuraminicAcid(), Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::dsLea) = GlycanSpec((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])

_glycanspec(::nLc4) = GlycanSpec((Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::nLc5) = GlycanSpec((Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIHAntigen) = GlycanSpec((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIAAntigen) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIBAntigen) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Lex5) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Ley6) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Lex7) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LexX8) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Ley8) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LeyX9) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::Lex9) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LexX10) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LexXX11) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::nLc6) = GlycanSpec((Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIH7Antigen) = GlycanSpec((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIA8Antigen) = GlycanSpec((Nacetylgalactosamine(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::TypeIIB8Antigen) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::IAntigen) = GlycanSpec((Galactose(), Nacetylglucosamine(), GlycanSpec((Galactose(), Nacetylglucosamine()), [β(1) => lk(4), β(1) => lk(6)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(6), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LeaX) = GlycanSpec((Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::LebX) = GlycanSpec((Fucose(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::sLex6) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::sLeaX) = GlycanSpec((NacetylneuraminicAcid(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), GlycanSpec((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::PAntigen) = GlycanSpec((Nacetylgalactosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
_glycanspec(::nLM1) = GlycanSpec((NacetylneuraminicAcid(), Nacetylgalactosamine(), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
