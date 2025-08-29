# MassSpecChemicals
parse_chemical(::Type{<: Saccharide}, x) = parse_glycan(x)

getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractPentose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
getchemicalattr(mono::Glucose, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "Glc"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxyhexose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
# getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxyhexose} = add_sub(mono, string("d", uppercasefirst(string(T.name.name)[begin + 5:begin + 7])))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxypentose} = add_sub(mono, add_ldform(mono, uppercasefirst(first(string(T.name.name), 3))))
# getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractDeoxypentose} = add_sub(mono, string("d", uppercasefirst(string(T.name.name)[begin + 5:begin + 7])))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexosamine} = add_sub(mono, add_ldform(mono, string(uppercasefirst(first(string(T.name.name), 3)), "N")))
getchemicalattr(mono::Glucosamine, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "GlcN"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractNacetylhexosamine} = add_sub(mono, add_ldform(mono, string(uppercasefirst(string(T.name.name)[begin + 7:begin + 9]), "NAc")))
getchemicalattr(mono::Nacetylglucosamine, ::Val{:abbreviation}; kwargs...) = add_ldform(mono, add_sub(mono, "GlcNAc"))
getchemicalattr(mono::T, ::Val{:abbreviation}; kwargs...) where {T <: AbstractHexuronicAcid} = add_sub(mono, add_ldform(mono, string(uppercasefirst(first(string(T.name.name), 3)), "A")))
getchemicalattr(mono::GlucuronicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "GlcA"))
getchemicalattr(mono::NeuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu"))
getchemicalattr(mono::NacetylneuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu5Ac"))
getchemicalattr(mono::NglycolylneuraminicAcid, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Neu5Gc"))
getchemicalattr(mono::Kdn, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Kdn"))
getchemicalattr(mono::Inositol, ::Val{:abbreviation}; kwargs...) = add_sub(mono, add_ldform(mono, "Ino"))
# getchemicalattr(mono::Sulfoquinovose, ::Val{:abbreviation}; kwargs...) = add_sub(mono, "Sulfoquinovose")

getchemicalattr(::GM4, ::Val{:abbreviation}; kwargs...) = "GM4"
getchemicalattr(::SM4, ::Val{:abbreviation}; kwargs...) = "SM4"
getchemicalattr(::Lac, ::Val{:abbreviation}; kwargs...) = "Lac"
getchemicalattr(::T, ::Val{:abbreviation}; kwargs...) where {T <: Ganglioseries} = string(T)
getchemicalattr(::T, ::Val{:abbreviation}; kwargs...) where {T <: Globoseries} = string(T)
getchemicalattr(::T, ::Val{:abbreviation}; kwargs...) where {T <: Isogloboseries} = string(T)
getchemicalattr(::T, ::Val{:abbreviation}; kwargs...) where {T <: Lactoseries} = string(T)
getchemicalattr(::T, ::Val{:abbreviation}; kwargs...) where {T <: Neolactoseries} = string(T)
getchemicalattr(x::SM1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "SM1", 2)
getchemicalattr(x::GM1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "GM1", 3)
getchemicalattr(x::GD1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "GD1", 5)
getchemicalattr(x::GT1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "GT1", 5)
getchemicalattr(x::GQ1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "GQ1", 4)
getchemicalattr(x::GP1, ::Val{:abbreviation}; kwargs...) = isomername(chemicalabbr.(x.isomer), "GP1", 2)

getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: Monosaccharide} = add_name(mono, add_ldform(mono, string(T.name.name)))
getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: AbstractNacetylhexosamine} = add_name(mono, add_ldform(mono, replace(string(T.name.name), "Nacetyl" => "N-acetyl")))
getchemicalattr(mono::T, ::Val{:name}; kwargs...) where {T <: AbstractHexuronicAcid} = add_name(mono, add_ldform(mono, replace(string(T.name.name), "Acid" => " acid")))
getchemicalattr(mono::NacetylneuraminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "N-acetylneuraminic acid"))
getchemicalattr(mono::NglycolylneuraminicAcid, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "N-glycolylneuraminic acid"))
getchemicalattr(mono::Kdn, ::Val{:name}; kwargs...) = add_name(mono, add_ldform(mono, "Keto-deoxy-glycero-galactonononic acid"))

getchemicalattr(::GM4, ::Val{:name}; kwargs...) = "GM4"
getchemicalattr(::SM4, ::Val{:name}; kwargs...) = "SM4"
getchemicalattr(::Lac, ::Val{:name}; kwargs...) = "Lac"
getchemicalattr(::T, ::Val{:name}; kwargs...) where {T <: Ganglioseries} = string(T)
getchemicalattr(::T, ::Val{:name}; kwargs...) where {T <: Globoseries} = string(T)
getchemicalattr(::T, ::Val{:name}; kwargs...) where {T <: Isogloboseries} = string(T)
getchemicalattr(::T, ::Val{:name}; kwargs...) where {T <: Lactoseries} = string(T)
getchemicalattr(::T, ::Val{:name}; kwargs...) where {T <: Neolactoseries} = string(T)
getchemicalattr(x::SM1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "SM1", 2)
getchemicalattr(x::GM1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "GM1", 3)
getchemicalattr(x::GD1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "GD1", 5)
getchemicalattr(x::GT1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "GT1", 5)
getchemicalattr(x::GQ1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "GQ1", 4)
getchemicalattr(x::GP1, ::Val{:name}; kwargs...) = isomername(chemicalname.(x.isomer), "GP1", 2)

getchemicalattr(glycan::GlyComp, ::Val{:abbreviation}; kwargs...) = join(map(glycan.comp) do (x, n)
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
end, "")

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

# concatchemical, makechemical for Glycomp