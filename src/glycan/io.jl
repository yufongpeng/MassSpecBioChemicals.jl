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

for x in SNFG_AA
    T = typeof(parse_aa_fg(x))
    @eval snfg_abbr(aa::$T) = letter3_abbr(originalmolecule(aa))
end
snfg_abbr(x) = chemicalabbr(x)
snfg_abbr(x::Amino) = "N"

const MONO_STRUCT = Dict{String, Type{<: Monosaccharide}}(
    "Hex"   => Hex,
    "Glc"   => Glc,
    "Gal"   => Gal,
    "Man"   => Man,
    "dHex"  => dHex,
    "Fuc"   => Fuc,
    "HexN"  => HexN,
    "GlcN"  => GlcN,
    "GalN"  => GalN,
    "ManN"  => ManN,
    "HexNAc"    => HexNAc,
    "GlcNAc"    => GlcNAc,
    "GalNAc"    => GalNAc,
    "ManNAc"    => ManNAc,
    "HexA"  => HexA,
    "GlcA"  => GlcA,
    "GalA"  => GalA,
    "ManA"  => ManA,
    "Neu"   => Neu,
    "Neu5Ac"    => Neu5Ac,
    "NeuAc"     => NeuAc,
    "Neu5Gc"    => Neu5Gc,
    "NeuGc"     => NeuGc,
    "Kdb"   => Kdb,
    "Pen"   => Pen,
    "Rib"   => Rib,
    "Ara"   => Ara,
    "Xyl"   => Xyl,
    "dPen"  => dPen,
    "dRib"  => dRib,
    "Ino"  => Ino,
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
parse_chemical(::Type{<: Saccharide}, x) = parse_glycan(x)
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

repr_linkage(l::Anomerposition) = l.position > 0 ? string(Int(l.position)) : ""
repr_linkage(l::Alphaposition) = l.position > 0 ? string("α", Int(l.position)) : "α"
repr_linkage(l::Betaposition) = l.position > 0 ? string("β", Int(l.position)) : "β"

function chemicalname(mono::T) where {T <: Monosaccharide}
    head = string(T.name.name)
    isnothing(mono.substituent) && return head 
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
chemicalname(::GM4) = "GM4"
chemicalname(::SM4) = "SM4"
chemicalname(::Lac) = "Lac"
chemicalname(::T) where {T <: Ganglioseries} = string(T)
chemicalname(::T) where {T <: Globoseries} = string(T)
chemicalname(::T) where {T <: Isogloboseries} = string(T)
chemicalname(::T) where {T <: Lactoseries} = string(T)
chemicalname(x::SM1) = isomername(x, "SM1", 2)
chemicalname(x::GM1) = isomername(x, "GM1", 3)
chemicalname(x::GD1) = isomername(x, "GD1", 5)
chemicalname(x::GT1) = isomername(x, "GT1", 5)
chemicalname(x::GQ1) = isomername(x, "GQ1", 4)
chemicalname(x::GP1) = isomername(x, "GP1", 2)
function isomername(x, c, n)
    if length(x.isomer) >= n
        c
    else
        string(c, "?(", join([replace(chemicalname(i), c => "") for i in x.isomer], ","), ")")
    end
end
chemicalname(glycan::GlyComp) = join(map(glycan.comp) do (x, n)
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
end, "")

function chemicalname(glycan::Glycan)
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
                s = string(s, cs, xs)
            else
                s = string(s, cs, "(", xs, ")")
            end
        end
    end
    if length(glycan.linkage) == length(glycan.chain) - 1
        s = string(s, chemicalname(last(glycan.chain)))
    end
    s
end

# specify num -> num
# default 6 -> 1 === dehydrogen -> dehydroxy 
# other than 6, => reverse smiles of sub
function repr_smiles(sugar::Hex)
    isnothing(sugar.substituent) && return "OC{6}C{5}(O1)C{4}(O)C{3}(O)C{2}(O)C{1}1O"
end
function repr_smiles(sugar::Glc)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end
function repr_smiles(sugar::Gal)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end