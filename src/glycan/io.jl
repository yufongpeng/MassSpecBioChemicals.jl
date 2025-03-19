const SNFG_AA = [    
    "Ala", 
    "Arg", 
    "Asn", 
    "Asp", 
    "Cys", 
    "Gln", 
    "Glu", 
    "Gly", 
    "His", 
    "Ile", 
    "Leu", 
    "Lys", 
    "Met", 
    "Phe", 
    "Pro", 
    "Ser", 
    "Thr", 
    "Trp", 
    "Tyr", 
    "Val", 
    "Orn"
]

# SNFG SPECIFIC functionalgroup ?
const SNFG_SUB = Dict{String, FunctionalGroup}(
    "Me"    => Methyl(),
    "Et"    => Ethyl(),
    "Fo"    => Formyl(),
    "Ac"    => Acetyl(),
    (aa => parse_aa_fg(aa) for aa in SNFG_AA)...,
    # Am, AmMe, AmMe2, AA2Ac, 5Glu
    "Gc"    => Glycoyl(),
    "Gr"    => Glyceryl(),
    # Gr2,3Me2
    # Hb 
    "Lt"    => Lactyl(),
    "NAc"   => Nacetyl(), # -OH
    "NFo"   => Nformyl(), # -OH
    "N"     => Amino(), # -OH 
    "P"     => Phosphate(),
    "Py"    => Pyruvyl(),
    # Pyr 
    "S"     => Sulfate(),
    "Tau"   => Tauryl()
)

"""
    snfg_abbr(x)

SNFG abbreviation
"""
snfg_abbr(x) = chemicalabbr(x)
for x in SNFG_AA
    T = typeof(parse_aa_fg(x))
    @eval snfg_abbr(aa::$T) = letter3_abbr(originalmolecule(aa))
end
snfg_abbr(x::Amino) = "N"

const MONO_STRUCT = Dict{String, Type{<: Monosaccharide}}(
    "Hex"   => Hexose,
    "Glc"   => Glucose,
    "Gal"   => Galactose,
    "Man"   => Mannose,
    "dHex"  => Deoxyhexose,
    "Fuc"   => Fucose,
    "HexN"  => Hexosamine,
    "GlcN"  => Glucosamine,
    "GalN"  => Galactosamine,
    "ManN"  => Mannosamine,
    "HexNAc"    => Nacetylhexosamine,
    "GlcNAc"    => Nacetylgluosamine,
    "GalNAc"    => Nacetylgalactosamine,
    "ManNAc"    => Nacetylmannosamine,
    "HexA"  => HexuronicAcid,
    "GlcA"  => GlucuronicAcid,
    "GalA"  => GalacturonicAcid,
    "ManA"  => MannuronicAcid,
    "Neu"   => NeuraminicAcid,
    "Neu5Ac"    => NacetylneuraminicAcid,
    "NeuAc"     => NacetylneuraminicAcid,
    "Neu5Gc"    => NglycolylneuraminicAcid,
    "NeuGc"     => NglycolylneuraminicAcid,
    "Kdn"   => Kdn,
    "Pen"   => Pentose,
    "Rib"   => Ribose,
    "Ara"   => Arabinose,
    "Xyl"   => Xylose,
    "dPen"  => Deoxypentose,
    "dRib"  => Deoxyribose,
    "Ino"  => Inositol,
    "Sulfoquinovose"    => Sulfoquinovose
)

const MONOSACCHRIDE = collect(keys(MONO_STRUCT))

function push_new_glycan!(cs, ls, x, c, d)
    mono, mp, mi, mj = x
    push!(cs, c)
    i = isempty(mi) ? first(d) : parse(UInt8, mi)
    lf = (isempty(mp) || mp == "(") ? Anomerposition(i) : (mp == "α" || mp == "(a") ? Alphaposition(i) : (mp == "β" || mp == "(b") ? Betaposition(i) : throw(ArgumentError("Invalid glycosidic linkage."))
    rt = isempty(mj) ? Linkageposition(last(d)) : Linkageposition(parse(UInt8, mj))
    push!(ls, lf => rt)
    isempty(mj)
end

function parse_glycan(s::AbstractString)
    p = 1
    cs = Saccharide[]
    ls = Pair{AbstractAnomerposition, Linkageposition}[]
    linktocheck = nothing
    @inbounds for m in eachmatch(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])", s)
        n = m.match.offset
        if n > p
            rs = collect(eachmatch(r"([^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:n]))
            ms = get_monosaccharide.(rs)
            if !isnothing(linktocheck)
                dp = dehydrogenposition(first(ms))
                ls[linktocheck] = first(ls[linktocheck]) => dp
            end 
            dmp = [dehydroxyposition(a) => dehydrogenposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
            push!(dmp, dehydroxyposition(last(ms)) => 0x00)
            trace = false
            for (x, c, d) in zip(rs, ms, dmp)
                trace = push_new_glycan!(cs, ls, x, c, d)
            end
            linktocheck = trace ? lastindex(ls) : nothing
        end
        push!(cs, parse_glycan(first(match(r"^\[?(.*?)\]?$", m.match))))
        push!(ls, last(last(cs).linkage))
        p = n + m.match.ncodeunits + 1
    end
    li = p
    rs = collect(eachmatch(r"([^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:end]))
    ms = get_monosaccharide.(rs)
    if !isnothing(linktocheck)
        dp = dehydrogenposition(first(ms))
        ls[linktocheck] = first(ls[linktocheck]) => Linkageposition(dp)
    end 
    dmp = [dehydroxyposition(a) => dehydrogenposition(b) for (a, b) in IterTools.partition(ms, 2, 1)]
    push!(dmp, dehydroxyposition(last(ms)) => 0x00)
    trace = false
    for (x, c, d) in zip(rs, ms, dmp)
        trace = push_new_glycan!(cs, ls, x, c, d)
        li = p + x.match.offset + x.match.ncodeunits
    end
    linktocheck = trace ? lastindex(ls) : nothing
    if li < lastindex(s)
        push!(cs, parse_monosaccharide(first(match(r"^\[?(.*?)\]?$", s[li:end]))))
        if !isnothing(linktocheck)
            ls[linktocheck] = first(ls[linktocheck]) => dehydrogenposition(last(cs))
        end 
    end
    Glycan((cs..., ), ls)
end

function split_snfg_sub(mod)
    s = Int[]
    e = Int[]
    dp = 0
    start = false
    for i in eachindex(mod)
        @inbounds x = mod[i:i]
        if dp == 0 && (x == "N" || !isnothing(match(r"[\d,?]", x)))
            if start 
                e[end] = i
                continue
            else
                push!(s, i)
                push!(e, i)
                start = true
                continue
            end
        end
        start = false
        @match x begin
            "(" => (dp += 1)
            ")" => (dp -= 1)
            _   => () 
        end
        e[end] = i
    end
    @inbounds [mod[i:j] for (i, j) in zip(s, e)]
end

function parse_monosaccharide(s)
    m = ""
    mono = nothing
    #return Hex()
    s = string(s)
    @inbounds for i in Iterators.reverse(eachindex(s))
        mono = s[begin:i]
        if haskey(MONO_STRUCT, mono)
            if i == lastindex(s)
                return MONO_STRUCT[mono]()
            else
                m = s[nextind(s, i):end]
                break
            end
        end
    end
    isnothing(mono) && throw(ArgumentError("Invalid monosacchride, \"$s\""))
    subs = split_snfg_sub(m)
    allsub = mapreduce(vcat, subs) do sub
        p, a = match(r"^([\d,?]*)(.*)$", sub)
        if startswith(a, "(") && endswith(a, r"\d\)")
            a, n = match(r"^\((.*[^\d])(\d+)\)$", a)
            # ignore p
            ps = zeros(UInt8, parse(Int, n))
        else
            ps = split(p, ",")
            replace!(ps, "?" => "0")
            ps = parse.(UInt8, ps)
        end
        ps .=> a
    end
    if all(iszero ∘ first, allsub)
        sk = last.(allsub)
        a = sort(unique(sk))
        substituent = [SNFG_SUB[x] => count(==(x), sk) for x in a]
    else
        sort!(allsub; by = reverse)
        substituent = [first(x) => SNFG_SUB[last(x)] for x in allsub]
    end
    MONO_STRUCT[mono](substituent)
end

function get_monosaccharide(x)
    m, mp, mi, mj = x
    parse_monosaccharide(first(match(r"^\[?(.*?)\]?$", m)))
end

function add_sub(mono::T, head) where {T <: Monosaccharide}
    (isnothing(mono.substituent) || isempty(mono.substituent)) && return head
    if mono.substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        string(head, join([n == 1 ? string("?", snfg_abbr(m)) : string("?(", snfg_abbr(m), Int(n), ")") for (m, n) in mono.substituent], ""))
    else
        ps = ""
        pp = String[]
        lm = nothing 
        for (p, m) in mono.substituent
            if isnothing(lm)
                lm = m
                push!(pp, p == 0 ? "?" : string(Int(p)))
            elseif m == lm
                push!(pp, p == 0 ? "?" : string(Int(p)))
            else
                ps = string(ps, join(pp, ","), snfg_abbr(lm))
                lm = m 
                pp = String[p == 0 ? "?" : string(Int(p))]
            end
        end
        ps = string(ps, join(pp, ","), snfg_abbr(lm))
        string(head, ps)
    end
end

function add_name(mono::T, head) where {T <: Monosaccharide}
    (isnothing(mono.substituent) || isempty(mono.substituent)) && return head
    if mono.substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        string(join([n == 1 ? string("?", chemicalname(m)) : string("?(", chemicalname(m), Int(n), ")") for (m, n) in mono.substituent], "-"), "-", string(T.name.name))
    else
        ps = ""
        pp = String[]
        lm = nothing 
        for (p, m) in mono.substituent
            if isnothing(lm)
                lm = m
                push!(pp, p == 0 ? "?" : string(Int(p)))
            elseif m == lm
                push!(pp, p == 0 ? "?" : string(Int(p)))
            else
                ps = string(ps, join(pp, ","), chemicalname(lm), "-")
                lm = m 
                pp = String[p == 0 ? "?" : string(Int(p))]
            end
        end
        ps = string(ps, join(pp, ","), chemicalname(lm), "-")
        string(ps, head)
    end
end

function isomername(isomer, c, n)
    if length(isomer) >= n
        c
    else
        string(c, "?(", join([replace(i, c => "") for i in isomer], ","), ")")
    end
end
