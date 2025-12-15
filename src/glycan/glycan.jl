generic_glycan(sugar::Glycan) = sugar
generic_glycan(sugar::SeriesGlycan{T, Nothing}) where T = generic_glycan(sugar.type) 
function generic_glycan(sugar::SeriesGlycan{G, T}) where {G, T <: Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}}}
    GlyComp(vcat(convert(Vector{Pair{Any, UInt8}}, getchaincomponent(GlyComp(SeriesGlycan(sugar.type, nothing)))), sugar.substituent))
end
function generic_glycan(sugar::SeriesGlycan{G, T}) where {G, T <: Vector{<: Pair{<: Monosaccharide, UInt8}}}
    g = generic_glycan(SeriesGlycan(sugar.type, nothing))
    if g isa Glycan && any(x -> first(x) isa SialicAcid, sugar.substituent)
        i = findall(x -> first(x) isa SialicAcid, sugar.substituent)
        sa = first(first(sugar.substituent[i]))
        if allequal(first, sugar.substituent[i]) && nsialicacid(g) == sum(last, sugar.substituent[i])
            replace_allsa!(g, sa)
        else 
            replace_sa_add_mono!(GlyComp(g), sugar.substituent)
        end
    else
        GlyComp(g, sugar.substituent)
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

function generic_glycan(sugar::SeriesGlycan{G, T}) where {G, T <: Vector{<: Pair{<: Tuple{UInt8, <: Vector}}}}
    g = generic_glycan(sugar.type) 
    s = [[deepcopy(collect(x)), y] for (x, y) in sugar.substituent]
    sort!(s; by = first ∘ first)
    replace_monosaccharide!(g, s)
end
function replace_monosaccharide!(sugar::Glycan, s)
    chain = reverse!(collect(Saccharide, getchaincomponent(sugar)))
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
        elseif mono isa Glycan 
            push!(id, i - 1)
            continue
        else
            chain[i] = replace_monosaccharide!(chain[i], s)
        end
    end
    Glycan((reverse!(chain)..., ), getchainlinkage(sugar))
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

generic_glycan(::GM4) = Glycan((NacetylneuraminicAcid(), Galactose()), [α(2) => lk(3)])
generic_glycan(::SM4) = Glycan((Galactose([0x03 => Sulfo()]), ), [])
generic_glycan(::Lac) = Glycan((Galactose(), Glucose()), [β(1) => lk(4)])
generic_glycan(::SM3) = Glycan((Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(4)])
generic_glycan(::SM2) = Glycan((Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::SM1a) = Glycan((Galactose(), Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::SM1b) = Glycan((Galactose([0x03 => Sulfo()]), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::SB1a) = Glycan((Galactose([0x03 => Sulfo()]), Nacetylgalactosamine(), Galactose([0x03 => Sulfo()]), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GA2) = Glycan((Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GA1) = Glycan((Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GM1b) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GM1α) = Glycan((Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GD1c) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GD1α) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((NacetylneuraminicAcid(),), [α(2) => lk(6)]), Nacetylgalactosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GM3) = Glycan((NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GM2) = Glycan((Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GM1a) = Glycan((Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GD1a) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GD1aα) = Glycan((Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT1a) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT1aα) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GD3) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GD2) = Glycan((Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GD1b) = Glycan((Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT1b) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT1bα) = Glycan((Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GQ1b) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GQ1bα) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT3) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT2) = Glycan((Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GT1c) = Glycan((Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GQ1c) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GQ1cα) = Glycan((Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GP1c) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])
generic_glycan(::GP1cα) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((NacetylneuraminicAcid(), ), [α(2) => lk(6)]), Nacetylgalactosamine(), Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), NacetylneuraminicAcid()), [α(2) => lk(8), α(2) => lk(8), α(2) => lk(3)]), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(2) => lk(6), β(1) => lk(4), α(2) => lk(3), β(1) => lk(4)])

generic_glycan(::SM1) = GlyComp([Nacetylgalactosamine() => 0x01, Galactose([0x03 => Sulfo()]) => 0x01, Galactose() => 0x01, Glucose() => 0x01])
generic_glycan(::GM1) = GlyComp([NacetylneuraminicAcid() => 0x01, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
generic_glycan(::GD1) = GlyComp([NacetylneuraminicAcid() => 0x02, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
generic_glycan(::GT1) = GlyComp([NacetylneuraminicAcid() => 0x03, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
generic_glycan(::GQ1) = GlyComp([NacetylneuraminicAcid() => 0x04, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])
generic_glycan(::GP1) = GlyComp([NacetylneuraminicAcid() => 0x05, Nacetylgalactosamine() => 0x01, Galactose() => 0x02, Glucose() => 0x01])

generic_glycan(::Gb3) = Glycan((Galactose(), Galactose(), Glucose()), [α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::Gb4) = Glycan((Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::Gb5) = Glycan((Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::ForssmanAntigen) = Glycan((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::ParaForssmanAntigen) = Glycan((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::BranchedForssmanAntigen) = Glycan((Nacetylgalactosamine(), Nacetylgalactosamine(), Glycan((Galactose(), Nacetylgalactosamine()), [β(1) => lk(3), β(1) => lk(4)]), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::SSEA4Antigen) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::SSEA4α6Antigen) = Glycan((NacetylneuraminicAcid(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(2) => lk(6), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::TypeIVHAntigen) = Glycan((Fucose(), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::TypeIVAAntigen) = Glycan((Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::TypeIVBAntigen) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])
generic_glycan(::GloboLex9) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Glycan((Galactose(), ), [β(1) => lk(3)]), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(6), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(4)])

generic_glycan(::iGb3) = Glycan((Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(4)])
generic_glycan(::iGb4) = Glycan((Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])
generic_glycan(::iGb5) = Glycan((Galactose(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])
generic_glycan(::ForssmanlikeiGb4) = Glycan((Nacetylgalactosamine(), Nacetylgalactosamine(), Galactose(), Galactose(), Glucose()), [α(1) => lk(3), β(1) => lk(3), α(1) => lk(3), β(1) => lk(4)])

generic_glycan(::Lc3) = Glycan((Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Lc4) = Glycan((Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Lea) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIHAntigen) = Glycan((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Leb) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIAAntigen) = Glycan((Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Aleb) = Glycan((Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIBAntigen) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Bleb) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LexA) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LeyA) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LM1) = Glycan((NacetylneuraminicAcid(), Nacetylgalactosamine(), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::sLea) = Glycan((NacetylneuraminicAcid(), Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::dsLea) = Glycan((NacetylneuraminicAcid(), NacetylneuraminicAcid(), Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(8), α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])

generic_glycan(::nLc4) = Glycan((Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::nLc5) = Glycan((Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIHAntigen) = Glycan((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIAAntigen) = Glycan((Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIBAntigen) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Lex5) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Ley6) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Lex7) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LexX8) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Ley8) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LeyX9) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::Lex9) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LexX10) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LexXX11) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::nLc6) = Glycan((Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIH7Antigen) = Glycan((Fucose(), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIA8Antigen) = Glycan((Nacetylgalactosamine(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::TypeIIB8Antigen) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(2)]), Galactose(), Nacetylglucosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(3), α(1) => lk(2), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::IAntigen) = Glycan((Galactose(), Nacetylglucosamine(), Glycan((Galactose(), Nacetylglucosamine()), [β(1) => lk(4), β(1) => lk(6)]), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(4), β(1) => lk(3), β(1) => lk(6), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LeaX) = Glycan((Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::LebX) = Glycan((Fucose(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(1) => lk(2), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::sLex6) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4), α(1) => lk(3), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::sLeaX) = Glycan((NacetylneuraminicAcid(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(4)]), Nacetylglucosamine(), Galactose(), Glycan((Fucose(), ), [α(1) => lk(3)]), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(3), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4), α(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::PAntigen) = Glycan((Nacetylgalactosamine(), Galactose(), Nacetylglucosamine(), Galactose(), Glucose()), [β(1) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
generic_glycan(::nLM1) = Glycan((NacetylneuraminicAcid(), Nacetylgalactosamine(), Nacetylglucosamine(), Galactose(), Glucose()), [α(2) => lk(3), β(1) => lk(4), β(1) => lk(3), β(1) => lk(4)])
