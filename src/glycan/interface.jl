# MassSpecChemicals
parse_chemical(::Type{<: Saccharide}, x) = parse_glycan(x)

getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::Glucose, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "Glc"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexosamine} = add_sub(mono, add_ldform(mono, string(uppercasefirst(first(string(T.name.name), 3)), "N")))
getchemicalattr(mono::Glucosamine, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "GlcN"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractNacetylhexosamine} = add_sub(mono, add_ldform(mono, string(uppercasefirst(string(T.name.name)[begin + 7:begin + 9]), "NAc")))
getchemicalattr(mono::Nacetylglucosamine, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "GlcNAc"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexuronicAcid} = add_sub(mono, add_ldform(mono, string(uppercasefirst(first(string(T.name.name), 3)), "A")))
getchemicalattr(mono::GlucuronicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "GlcA"))
getchemicalattr(mono::IduronicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "IdoA"))

getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxyhexose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::Deoxyhexose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "dHex"))
getchemicalattr(mono::Sixdeoxyaltose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "6dAlt"))
getchemicalattr(mono::Sixdeoxytalose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "6dTal"))
getchemicalattr(mono::Sixdeoxygulose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "6dGul"))

getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractNacetyldeoxyhexosamine} = add_sub(mono, add_ldform(mono, string(uppercasefirst(string(T.name.name)[begin + 7:begin + 9]), "NAc")))
getchemicalattr(mono::Nacetyldeoxyhexosamine, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "dHexNAc"))
getchemicalattr(mono::Nacetylsixdeoxyaltosamine, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "6dAltNAc"))
getchemicalattr(mono::Nacetylsixdeoxytalosamine, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "6dTalNAc"))

getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDideoxyhexose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::Dideoxyhexose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "ddHex"))

getchemicalattr(mono::NeuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu"))
getchemicalattr(mono::NacetylneuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu"), [0x05 => Acetyl()])
getchemicalattr(mono::NglycolylneuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu"), [0x05 => Glycolyl()])
getchemicalattr(mono::Kdn, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Kdn"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: DiaminodideoxynonulosonicAcid} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::FourepilegionaminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "4eLeg"))
getchemicalattr(mono::EightepilegionaminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "8eLeg"))
getchemicalattr(mono::EightepiacinetaminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "8eAci"))

getchemicalattr(mono::Inositol, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Ino"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractPentose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxypentose} = add_sub(mono, add_ldform(mono, string("d", uppercasefirst(first(string(T.name.name), 3)))))

getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: Monosaccharide} = add_name(mono, add_ldform(mono, string(T.name.name)))
getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: AbstractNacetylhexosamine} = add_name(mono, string("N-acetyl-", add_ldform(mono, replace(string(T.name.name), "Nacetyl" => ""))))
getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: AbstractHexuronicAcid} = add_name(mono, add_ldform(mono, replace(string(T.name.name), "Acid" => " acid")))
getchemicalattr(mono::Sixdeoxyaltose, ::Val{:name}; kwargs...) = add_name(mono, string("6-deoxy-", add_ldform(mono, "altose")))
getchemicalattr(mono::Sixdeoxytalose, ::Val{:name}; kwargs...) = add_name(mono, string("6-deoxy-", add_ldform(mono, "talose")))
getchemicalattr(mono::Sixdeoxygulose, ::Val{:name}; kwargs...) = add_name(mono, string("6-deoxy-", add_ldform(mono, "gulose")))
getchemicalattr(mono::Nacetylsixdeoxyaltosamine, ::Val{:name}; kwargs...) = add_name(mono, string("N-acetyl-6-deoxy-", add_ldform(mono, "altosamine")))
getchemicalattr(mono::Nacetylsixdeoxytalosamine, ::Val{:name}; kwargs...) = add_name(mono, string("N-acetyl-6-deoxy-", add_ldform(mono, "talosamine")))
getchemicalattr(mono::NacetylneuraminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "N-acetylneuraminic acid"))
getchemicalattr(mono::NglycolylneuraminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "N-glycolylneuraminic acid"))
getchemicalattr(mono::Kdn, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "Keto-deoxy-glycero-galactonononic acid"))
getchemicalattr(mono::FourepilegionaminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "4-epi-legionaminic acid"))
getchemicalattr(mono::EightepilegionaminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "8-epi-legionaminic acid"))
getchemicalattr(mono::EightepiacinetaminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "8-epi-acinetaminic acid"))

getchemicalattr(x::SeriesGlycan{T, Nothing}, ::Val{:name}; kwargs...) where T = glycanabbr(x.type)
getchemicalattr(x::SeriesGlycan{T, Nothing}, ::Val{:abbreviation}; kwargs...) where T = glycanabbr(x.type)
getchemicalattr(x::SeriesGlycan, ::Val{:name}; kwargs...) = string(glycanabbr(x.type), "(", glysubname(x.substituent), ")")
getchemicalattr(x::SeriesGlycan, ::Val{:abbreviation}; kwargs...) = string(glycanabbr(x.type), "(", glysubabbr(x.substituent), ")")
function getchemicalattr(glycan::GlyComp, ::Val{:abbreviation}; kwargs...) 
    v = map(glycan.comp) do (x, n)
        c = chemicalabbr(x)
        if n > 1
            if endswith(c, r"\d")
                c = string("(", c, ")", Int(n))
            else
                c = string(c, Int(n))
            end
        elseif endswith(c, r"\d")
            c = string("(", c, ")")
        end
        c
    end
    i = findfirst(x -> first(x) isa AbstractFunctionalGroup, glycan.comp)
    isnothing(i) ? join(v, "") : string(join(v[begin:i - 1], ""), "_", join(v[i:end], ""))
end

function getchemicalattr(glycan::Glycan, ::Val{:abbreviation}; kwargs...)
    s = ""
    for (i, x) in enumerate(glycan.linkage)
        c = glycan.chain[i]
        xs = repr_linkage(x)
        cs = chemicalabbr(c)
        if c isa Glycan
            m = match(r"\(?((?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?$", cs)
            sp, si, sj = m
            mp, mi, mj = match(r"([αβab]*)(\d*)-(\d*)$", xs)
            if isempty(sp)
                sp = mp
            end
            if isempty(si)
                si = mi
            end
            if isempty(mj)
                mj = sj
            end
            x = string(sp, si, "-", mj)
            if startswith(x, "-") || startswith(x, "α") || startswith(x, "β")
                s = string(s, "[", cs[begin:m.match.offset], sp, si, "-", mj, "]")
            else
                s = string(s, "[", cs[begin:m.match.offset], "(", sp, si, "-", mj, ")]")
            end
        else
            if startswith(xs, "-") || startswith(xs, "α") || startswith(xs, "β")
                s = string(s, cs, xs)
            else
                s = string(s, cs, "(", xs, ")")
            end
        end
    end
    if length(glycan.linkage) == length(glycan.chain) - 1
        s = string(s, chemicalabbr(last(glycan.chain)))
    end
    s
end

getchemicalattr(glycan::GlyComp, ::Val{:name}; kwargs...) = join(map(glycan.comp) do (x, n)
    c = chemicalname(x)
    if n > 1
        if endswith(c, r"\d")
            c = string("(", c, ")", Int(n))
        else
            c = string(c, Int(n))
        end
    elseif endswith(c, r"\d")
        c = string("(", c, ")")
    end
    c
end, "_")

function getchemicalattr(glycan::Glycan, ::Val{:name}; kwargs...)
    s = ""
    for (i, x) in enumerate(glycan.linkage)
        c = glycan.chain[i]
        xs = repr_linkage(x)
        cs = chemicalname(c)
        if c isa Glycan
            m = match(r"\(?((?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?$", cs)
            sp, si, sj = m
            mp, mi, mj = match(r"([αβab]*)(\d*)-(\d*)$", xs)
            if isempty(sp)
                sp = mp
            end
            if isempty(si)
                si = mi
            end
            if isempty(mj)
                mj = sj
            end
            x = string(sp, si, "-", mj)
            if startswith(x, "-") || startswith(x, "α") || startswith(x, "β")
                s = string(s, "[", cs[begin:m.match.offset], sp, si, "-", mj, "]")
            else
                s = string(s, "[", cs[begin:m.match.offset], "(", sp, si, "-", mj, ")]")
            end
        else
            if startswith(xs, "-") || startswith(xs, "α") || startswith(xs, "β")
                if startswith(cs, r"\d") && endswith(s, r"\d")
                    s = string(s, "-", cs, xs)
                else
                    s = string(s, cs, xs)
                end
            else
                s = string(s, cs, "(", xs, ")")
            end
        end
    end
    if length(glycan.linkage) == length(glycan.chain) - 1
        c = chemicalname(last(glycan.chain))
        if startswith(c, r"\d") && endswith(s, r"\d")
            s = string(s, "-", c)
        else
            s = string(s, c)
        end
    end
    s
end

# specify num -> num
# default 6 -> 1 === dehydrogen -> dehydroxy 
# other than 6, => reverse smiles of sub
function getchemicalattr(sugar::Hexose, ::Val{:SMILES}; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}C{5}(O1)C{4}(O)C{3}(O)C{2}(O)C{1}1O"
end
function getchemicalattr(sugar::Glucose, ::Val{:SMILES}; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end
function getchemicalattr(sugar::Galactose, ::Val{:SMILES}; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end

getchemicalattr(sugar::SeriesGlycan, ::Val{:elements}; kwargs...) = chemicalelements(generic_glycan(sugar); kwargs...)
function getchemicalattr(sugar::Glycan, ::Val{:elements}; kwargs...)
    es = [losswaterelements(m, i, i % length(sugar)) for (i, m) in enumerate(getchaincomponent(sugar))]
    vcat(es...)
end

function getchemicalattr(sugar::GlyComp, ::Val{:elements}; kwargs...)
    es = mapreduce(vcat, getchaincomponent(sugar)) do x 
        repeat([first(x)], last(x))
    end
    es = [losswaterelements(m, i, i % length(sugar)) for (i, m) in enumerate(es)]
    vcat(es...)
end
function losswaterelements(sugar::Monosaccharide, p, q)
    e = chemicalelements(sugar)
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
monosacchrideelements(::AbstractHexose) = ["C" => 6, "H" => 12, "O" => 6]
monosacchrideelements(::AbstractHexosamine) = ["C" => 6, "H" => 13, "N" => 1, "O" => 5]
monosacchrideelements(::AbstractNacetylhexosamine) = ["C" => 8, "H" => 15, "N" => 1, "O" => 6]
monosacchrideelements(::AbstractHexuronicAcid) = ["C" => 7, "H" => 12, "O" => 7]
monosacchrideelements(::AbstractDeoxyhexose) = ["C" => 6, "H" => 12, "O" => 5]
monosacchrideelements(::AbstractDideoxyhexose) = ["C" => 6, "H" => 12, "O" => 4]
monosacchrideelements(::AbstractNacetyldeoxyhexosamine) = ["C" => 8, "H" => 15, "N" => 1, "O" => 5]
monosacchrideelements(::AbstractPentose) = ["C" => 5, "H" => 10, "O" => 5]
monosacchrideelements(::AbstractDeoxypentose) = ["C" => 5, "H" => 10, "O" => 4]
monosacchrideelements(::NeuraminicAcid) = ["C" => 9, "H" => 17, "N" => 1, "O" => 8]
monosacchrideelements(::NacetylneuraminicAcid) = ["C" => 11, "H" => 19, "N" => 1, "O" => 9]
monosacchrideelements(::NglycolylneuraminicAcid) = ["C" => 11, "H" => 19, "N" => 1, "O" => 10]
monosacchrideelements(::Kdn) = ["C" => 11, "H" => 18, "O" => 11]
monosacchrideelements(::DiaminodideoxynonulosonicAcid) = ["C" => 9, "H" => 18, "N" => 2, "O" => 6]

function getchemicalattr(sugar::Monosaccharide{D, P, T}, ::Val{:elements}; kwargs...) where {D, P, T}
    es = monosacchrideelements(sugar)
    (isnothing(sugar.substituent) || isempty(sugar.substituent)) && return es
    i = findfirst(x -> first(x) == "H", es)
    if mono.substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        m = mapreduce(+, sugar.substituent) do x 
            ntotallinkage(first(x)) * Int(last(x))
        end
        es[i] = "H" => (leat(es[i]) - m)
        es = vcat(es, map(sugar.substituent) do x 
            repeat(chemicalelements(first(x)); outer = Int(last(x)))
        end...)
    else
        es[i] = "H" => (last(es[i]) - sum(ntotallinkage ∘ last, sugar.substituent))
        es = vcat(es, chemicalelements.(last.(sugar.substituent))...)
    end
    es
end
function getchemicalattr(x::Substituent{<: Glycan}, ::Val{:elements}; kwargs...) 
	ec = [losswaterelements(m, i, i % length(parentchemical(x))) for (i, m) in enumerate(getchaincomponent(parentchemical(x)))]
    es = last(ec)
	ls = leavinggroupelements(leavinggroup(x))
	for (e, n) in ls 
		if n > 0
			i = findfirst(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = findfirst(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	vcat(ec...)
end
function getchemicalattr(x::Substituent{<: GlyComp}, ::Val{:elements}; kwargs...) 
	es = mapreduce(vcat, getchaincomponent(sugar)) do x 
        repeat([first(x)], last(x))
    end
    ec = [losswaterelements(m, i, i % length(sugar)) for (i, m) in enumerate(es)]
    es = last(ec)
	ls = leavinggroupelements(leavinggroup(x))
	for (e, n) in ls 
		if n > 0
			i = findfirst(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = findfirst(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	vcat(ec...)
end
# named glycan

getchemicalattr(sugar::Saccharide, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(sugar); unique, kwargs...)

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

getchaincomponent(sugar::SeriesGlycan) = getchaincomponent(generic_glycan(sugar))
composition(sugar::AbstractGlycan) = Dict(getchaincomponent(GlyComp(sugar))...)
composition(sugar::GlyComp) = Dict(getchaincomponent(sugar)...)
length(sugar::AbstractGlycan) = sum(length, getchaincomponent(sugar))
length(sugar::GlyComp) = length(getchaincomponent(sugar))
# concatchemical, makechemical for Glycomp