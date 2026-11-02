function make_monosaccharide(::Type{M}, substituent = nothing; new_attr = true, kwargs...) where {M <: Monosaccharide}
    kw = (; pairs(defaultchemicalkwargs(M, substituent))..., kwargs...)
    attr = Base.structdiff(kw, NamedTuple{(:name, :elements, :abbreviation)})
    M(substituent, String(kw.name), collect(Pair{String, Int}, kw.elements), String(kw.abbreviation),
      buildattribute(new_attr, pairs(attr)))
end

defaultchemicalkwargs(::Type{M}, substituent) where {M <: Monosaccharide} =
    (; name = chemicaldefaultname(M, substituent), elements = chemicaldefaultelements(M, substituent),
       abbreviation = chemicaldefaultabbr(M, substituent))
defaultchemicalkwargs(::Type{<: SeriesGlycan}, type, substituent) =
    (; name = _serieglycanname(type, substituent), elements = _seriesglycanelements(type, substituent),
       abbreviation = _serieglycanabbr(type, substituent))
defaultchemicalkwargs(::Type{<: Glycan}, chain, linkage) =
    (; name = _glycanname(chain, linkage), elements = _glycanelements(chain),
       abbreviation = _glycanabbr(chain, linkage))
defaultchemicalkwargs(::Type{<: GlyComp}, comp) =
    (; name = _glycompname(comp), elements = _glycompelements(comp), abbreviation = _glycompabbr(comp))

# Monosaccharide abbreviation
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Monosaccharide} = chemicaldefaultname(T, substituent)
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractHexose} = add_sub(substituent, add_ldform(T, uppercasefirst(first(string(T.name.name), 3))))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Glucose} = add_ldform(T, add_sub(substituent, "Glc"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractHexosamine} = add_sub(substituent, add_ldform(T, string(uppercasefirst(first(string(T.name.name), 3)), "N")))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Glucosamine} = add_ldform(T, add_sub(substituent, "GlcN"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractNacetylhexosamine} = add_sub(substituent, add_ldform(T, string(uppercasefirst(string(T.name.name)[begin + 7:begin + 9]), "NAc")))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Nacetylglucosamine} = add_ldform(T, add_sub(substituent, "GlcNAc"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractHexuronicAcid} = add_sub(substituent, add_ldform(T, string(uppercasefirst(first(string(T.name.name), 3)), "A")))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: GlucuronicAcid} = add_sub(substituent, add_ldform(T, "GlcA"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: IduronicAcid} = add_sub(substituent, add_ldform(T, "IdoA"))

chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractDeoxyhexose} = add_sub(substituent, add_ldform(T, uppercasefirst(first(string(T.name.name), 3))))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Deoxyhexose} = add_sub(substituent, add_ldform(T, "dHex"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Sixdeoxyaltose} = add_sub(substituent, add_ldform(T, "6dAlt"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Sixdeoxytalose} = add_sub(substituent, add_ldform(T, "6dTal"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Sixdeoxygulose} = add_sub(substituent, add_ldform(T, "6dGul"))

chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractNacetyldeoxyhexosamine} = add_sub(substituent, add_ldform(T, string(uppercasefirst(string(T.name.name)[begin + 7:begin + 9]), "NAc")))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Nacetyldeoxyhexosamine} = add_sub(substituent, add_ldform(T, "dHexNAc"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Nacetylsixdeoxyaltosamine} = add_sub(substituent, add_ldform(T, "6dAltNAc"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Nacetylsixdeoxytalosamine} = add_sub(substituent, add_ldform(T, "6dTalNAc"))

chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractDideoxyhexose} = add_sub(substituent, add_ldform(T, uppercasefirst(first(string(T.name.name), 3))))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Dideoxyhexose} = add_sub(substituent, add_ldform(T, "ddHex"))

chemicaldefaultabbr(::Type{T}, substituent) where {T <: NeuraminicAcid} = add_sub(substituent, add_ldform(T, "Neu"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: NacetylneuraminicAcid} = add_sub(substituent, add_ldform(T, "Neu"), [0x05 => Acetyl()])
chemicaldefaultabbr(::Type{T}, substituent) where {T <: NglycolylneuraminicAcid} = add_sub(substituent, add_ldform(T, "Neu"), [0x05 => Glycolyl()])
chemicaldefaultabbr(::Type{T}, substituent) where {T <: Kdn} = add_sub(substituent, add_ldform(T, "Kdn"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractdiaminodideoxynonulosonicAcid} = add_sub(substituent, add_ldform(T, uppercasefirst(first(string(T.name.name), 3))))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: FourepilegionaminicAcid} = add_sub(substituent, add_ldform(T, "4eLeg"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: EightepilegionaminicAcid} = add_sub(substituent, add_ldform(T, "8eLeg"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: EightepiacinetaminicAcid} = add_sub(substituent, add_ldform(T, "8eAci"))

chemicaldefaultabbr(::Type{T}, substituent) where {T <: Inositol} = add_sub(substituent, add_ldform(T, "Ino"))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractPentose} = add_sub(substituent, add_ldform(T, uppercasefirst(first(string(T.name.name), 3))))
chemicaldefaultabbr(::Type{T}, substituent) where {T <: AbstractDeoxypentose} = add_sub(substituent, add_ldform(T, string("d", uppercasefirst(first(string(T.name.name), 3)))))

# Monosaccharide name
chemicaldefaultname(::Type{T}, substituent) where {T <: Monosaccharide} = add_name(T, substituent, add_ldform(T, string(T.name.name)))
chemicaldefaultname(::Type{T}, substituent) where {T <: AbstractNacetylhexosamine} = add_name(T, substituent, string("N-acetyl-", add_ldform(T, replace(string(T.name.name), "Nacetyl" => ""))))
chemicaldefaultname(::Type{T}, substituent) where {T <: AbstractHexuronicAcid} = add_name(T, substituent, add_ldform(T, replace(string(T.name.name), "Acid" => " acid")))
chemicaldefaultname(::Type{T}, substituent) where {T <: Sixdeoxyaltose} = add_name(T, substituent, string("6-deoxy-", add_ldform(T, "altose")))
chemicaldefaultname(::Type{T}, substituent) where {T <: Sixdeoxytalose} = add_name(T, substituent, string("6-deoxy-", add_ldform(T, "talose")))
chemicaldefaultname(::Type{T}, substituent) where {T <: Sixdeoxygulose} = add_name(T, substituent, string("6-deoxy-", add_ldform(T, "gulose")))
chemicaldefaultname(::Type{T}, substituent) where {T <: Nacetylsixdeoxyaltosamine} = add_name(T, substituent, string("N-acetyl-6-deoxy-", add_ldform(T, "altosamine")))
chemicaldefaultname(::Type{T}, substituent) where {T <: Nacetylsixdeoxytalosamine} = add_name(T, substituent, string("N-acetyl-6-deoxy-", add_ldform(T, "talosamine")))
chemicaldefaultname(::Type{T}, substituent) where {T <: NacetylneuraminicAcid} = add_name(T, substituent, add_ldform(T, "N-acetylneuraminic acid"))
chemicaldefaultname(::Type{T}, substituent) where {T <: NglycolylneuraminicAcid} = add_name(T, substituent, add_ldform(T, "N-glycolylneuraminic acid"))
chemicaldefaultname(::Type{T}, substituent) where {T <: Kdn} = add_name(T, substituent, add_ldform(T, "Keto-deoxy-glycero-galactonononic acid"))
chemicaldefaultname(::Type{T}, substituent) where {T <: FourepilegionaminicAcid} = add_name(T, substituent, add_ldform(T, "4-epi-legionaminic acid"))
chemicaldefaultname(::Type{T}, substituent) where {T <: EightepilegionaminicAcid} = add_name(T, substituent, add_ldform(T, "8-epi-legionaminic acid"))
chemicaldefaultname(::Type{T}, substituent) where {T <: EightepiacinetaminicAcid} = add_name(T, substituent, add_ldform(T, "8-epi-acinetaminic acid"))

# Monosaccharide elements
monosacchrideelements(::Type{<: AbstractHexose}) = ["C" => 6, "H" => 12, "O" => 6]
monosacchrideelements(::Type{<: AbstractHexosamine}) = ["C" => 6, "H" => 13, "N" => 1, "O" => 5]
monosacchrideelements(::Type{<: AbstractNacetylhexosamine}) = ["C" => 8, "H" => 15, "N" => 1, "O" => 6]
monosacchrideelements(::Type{<: AbstractHexuronicAcid}) = ["C" => 7, "H" => 12, "O" => 7]
monosacchrideelements(::Type{<: AbstractDeoxyhexose}) = ["C" => 6, "H" => 12, "O" => 5]
monosacchrideelements(::Type{<: AbstractDideoxyhexose}) = ["C" => 6, "H" => 12, "O" => 4]
monosacchrideelements(::Type{<: AbstractNacetyldeoxyhexosamine}) = ["C" => 8, "H" => 15, "N" => 1, "O" => 5]
monosacchrideelements(::Type{<: AbstractPentose}) = ["C" => 5, "H" => 10, "O" => 5]
monosacchrideelements(::Type{<: AbstractDeoxypentose}) = ["C" => 5, "H" => 10, "O" => 4]
monosacchrideelements(::Type{<: NeuraminicAcid}) = ["C" => 9, "H" => 17, "N" => 1, "O" => 8]
monosacchrideelements(::Type{<: NacetylneuraminicAcid}) = ["C" => 11, "H" => 19, "N" => 1, "O" => 9]
monosacchrideelements(::Type{<: NglycolylneuraminicAcid}) = ["C" => 11, "H" => 19, "N" => 1, "O" => 10]
monosacchrideelements(::Type{<: Kdn}) = ["C" => 11, "H" => 18, "O" => 11]
monosacchrideelements(::Type{<: AbstractdiaminodideoxynonulosonicAcid}) = ["C" => 9, "H" => 18, "N" => 2, "O" => 6]
monosacchrideelements(::Type{<: Inositol}) = ["C" => 6, "H" => 12, "O" => 6]

function chemicaldefaultelements(::Type{T}, substituent) where {T <: Monosaccharide}
    es = monosacchrideelements(T)
    (isnothing(substituent) || isempty(substituent)) && return es
    i = findfirst(x -> first(x) == "H", es)
    if substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        m = mapreduce(+, substituent) do x 
            ntotallinkage(first(x)) * Int(last(x))
        end
        es[i] = "H" => (last(es[i]) - m)
        es = vcat(es, map(substituent) do x 
            repeat(chemicalelements(first(x)); outer = Int(last(x)))
        end...)
    else
        es[i] = "H" => (last(es[i]) - sum(ntotallinkage ∘ last, substituent))
        es = vcat(es, chemicalelements.(last.(substituent))...)
    end
    es
end

# Glycan
_serieglycanname(type, substituent) = isnothing(substituent) ? glycanabbr(type) : string(glycanabbr(type), "(", glysubname(substituent), ")")
_serieglycanabbr(type, substituent) = isnothing(substituent) ? glycanabbr(type) : string(glycanabbr(type), "(", glysubabbr(substituent), ")")
# Elements are computed from the composition of the generic glycan; neither `Glycan` nor `GlyComp` is constructed
_seriesglycanelements(type, ::Nothing) = _glycompelements(_generic_composition(type))
function _seriesglycanelements(type, substituent::Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}})
    _glycompelements(vcat(_generic_composition(type), substituent))
end
function _seriesglycanelements(type, substituent::Vector{<: Pair{<: Monosaccharide, UInt8}})
    comp = _generic_composition(type)
    issa = x -> first(x) isa SialicAcid
    i = findall(issa, substituent)
    if !isempty(i) && !(type isa Union{SM1, GM1, GD1, GT1, GQ1, GP1})
        sa = first(first(substituent[i]))
        if allequal(first, substituent[i]) && sum(last, filter(issa, comp); init = 0x00) == sum(last, substituent[i])
            # all sialic acids are replaced
            return _glycompelements(vcat(filter(!issa, comp), [sa => sum(last, substituent[i])]))
        end
    end
    d = Dict{Monosaccharide, UInt8}(comp)
    for (m, n) in substituent
        get!(d, m, 0x00)
        if m isa SialicAcid && !(type isa Union{SM1, GM1, GD1, GT1, GQ1, GP1})
            d[NacetylneuraminicAcid()] -= n
        end
        d[m] += n
    end
    filter!(x -> last(x) > 0, d)
    _glycompelements(collect(pairs(d)))
end
function _seriesglycanelements(type, substituent::Vector{<: Pair{<: Tuple{UInt8, <: Vector}}})
    spec = replace_monosaccharide!(_glycanspec(type), [[deepcopy(collect(x)), y] for (x, y) in sort(substituent; by = x -> first(first(x)))])
    d = Dict{Monosaccharide, UInt8}()
    _push_composition!(d, spec)
    _glycompelements(collect(pairs(d)))
end

function _glycompabbr(comp)
    v = map(comp) do (x, n)
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
    i = findfirst(x -> first(x) isa AbstractFunctionalGroup, comp)
    isnothing(i) ? join(v, "") : string(join(v[begin:i - 1], ""), "_", join(v[i:end], ""))
end

function _glycanabbr(chain, linkage)
    s = ""
    for (i, x) in enumerate(linkage)
        c = chain[i]
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
    if length(linkage) == length(chain) - 1
        s = string(s, chemicalabbr(last(chain)))
    end
    s
end

_glycompname(comp) = join(map(comp) do (x, n)
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

function _glycanname(chain, linkage)
    s = ""
    for (i, x) in enumerate(linkage)
        c = chain[i]
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
    if length(linkage) == length(chain) - 1
        c = chemicalname(last(chain))
        if startswith(c, r"\d") && endswith(s, r"\d")
            s = string(s, "-", c)
        else
            s = string(s, c)
        end
    end
    s
end

function _glycanelements(chain)
    n = length(chain)
    reduce(vcat, (losswaterelements(m, i, i % n) for (i, m) in enumerate(chain)); init = Pair{String, Int}[])
end

function _glycompelements(comp)
    es = Pair{String, Int}[]
    nmono = 0
    for (x, n) in comp
        k = Int(n)
        if x isa AbstractFunctionalGroup
            # a substituent replaces hydrogens of the saccharides
            append!(es, ("H" => -ntotallinkage(x) * k, ))
            for (e, c) in chemicalelements(x)
                push!(es, e => c * k)
            end
        else
            nmono += k
            for (e, c) in chemicalelements(x)
                push!(es, e => c * k)
            end
        end
    end
    # each glycosidic bond loses one water
    nmono > 1 && push!(es, "H" => -2 * (nmono - 1), "O" => -(nmono - 1))
    merged = Pair{String, Int}[]
    for (e, c) in es
        i = findfirst(y -> first(y) == e, merged)
        isnothing(i) ? push!(merged, e => c) : (merged[i] = e => last(merged[i]) + c)
    end
    filter!(x -> last(x) != 0, merged)
end
