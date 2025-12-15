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
    "Gc"    => Glycolyl(),
    "Gr"    => Glyceryl{RSChirality}(),
    # Gr2,3Me2
    # Hb 
    "Lt"    => Lactyl{RSChirality}(),
    "NAc"   => Nacetyl(), # -OH
    "NFo"   => Nformyl(), # -OH
    "N"     => Amino(), # -OH 
    "P"     => Phosphoryl(),
    "Py"    => Pyruvyl(),
    # Pyr 
    "S"     => Sulfo(),
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
    "Man"   => Mannose,
    "Gal"   => Galactose,
    "Gul"   => Gulose,
    "Alt"   => Altose,
    "All"   => Allose,
    "Tal"   => Talose,
    "Ido"   => Idose,
    "Api"   => Apiose,
    "Fru"   => Fructose,
    "Tag"   => Tagatose,
    "Sor"   => Sorbose,
    "Psi"   => Psicose,
    "HexN"  => Hexosamine,
    "GlcN"  => Glucosamine,
    "ManN"  => Mannosamine,
    "GalN"  => Galactosamine,
    "GulN"  => Gulosamine,
    "AltN"  => Altosamine,
    "AllN"  => Allosamine,
    "TalN"  => Talosamine,
    "IdoN"  => Idosamine,
    "HexNAc"    => Nacetylhexosamine,
    "GlcNAc"    => Nacetylglucosamine,
    "ManNAc"    => Nacetylmannosamine,
    "GalNAc"    => Nacetylgalactosamine,
    "GulNAc"    => Nacetylgulosamine,
    "AltNAc"    => Nacetylaltosamine,
    "AllNAc"    => Nacetylallosamine,
    "TalNAc"    => Nacetyltalosamine,
    "IdoNAc"    => Nacetylidosamine,
    "HexA"  => HexuronicAcid,
    "GlcA"  => GlucuronicAcid,
    "ManA"  => MannuronicAcid,
    "GalA"  => GalacturonicAcid,
    "GulA"  => GuluronicAcid,
    "AltA"  => AlturonicAcid,
    "AllA"  => AlluronicAcid,
    "TalA"  => TaluronicAcid,
    "IdoA"  => IduronicAcid,
    "dHex"  => Deoxyhexose,
    "Qui"   => Quinovose,
    "Rha"   => Rhamnose,
    "Fuc"   => Fucose,
    "6dGul" => Sixdeoxygulose,
    "6dAlt" => Sixdeoxyaltose,
    "6dTal" => Sixdeoxytalose,
    "dHexNAc"   => Nacetyldeoxyhexosamine,
    "QuiNAc"    => Nacetylquinovosamine,
    "RhaNAc"    => Nacetylrhamnosamine,
    "FucNAc"    => Nacetylfucosamine,
    "6dAltNAc"  => Nacetylsixdeoxyaltosamine,
    "6dTalNAc"  => Nacetylsixdeoxytalosamine,
    "ddHex" => Dideoxyhexose,
    "oli"   => Olivose,
    "Tyv"   => Tyvelose,
    "Abe"   => Abequose,
    "Par"   => Paratose,
    "Col"   => Colitose,
    "Dig"   => Digitoxose,
    "Neu"   => NeuraminicAcid,
    "Neu5Ac"    => NacetylneuraminicAcid,
    "NeuAc"     => NacetylneuraminicAcid,
    "Neu5Gc"    => NglycolylneuraminicAcid,
    "NeuGc"     => NglycolylneuraminicAcid,
    "Kdn"   => Kdn,
    "Leg"   => LegionaminicAcid,
    "4eLeg" => FourepilegionaminicAcid,
    "8eLeg" => EightepilegionaminicAcid,
    "Aci"   => AcinetaminicAcid,
    "8eAci" => EightepiacinetaminicAcid,
    "Pse"   => PseudaminicAcid,
    "Pen"   => Pentose,
    "Ara"   => Arabinose,
    "Lyx"   => Lyxose,
    "Xyl"   => Xylose,
    "Rib"   => Ribose,
    "dPen"  => Deoxypentose,
    "dRib"  => Deoxyribose,
    "Ino"   => Inositol
)

const GLYCAN_STRUCT = Dict{String, Type{<: GlycanType}}(
    "GM4"   => GM4,
    "SM4"   => SM4,
    "Lac"   => Lac,
    "SM3"   => SM3,
    "SM2"   => SM2,
    "SM1a"  => SM1a,
    "SM1b"  => SM1b,
    "SB1a"  => SB1a,
    "SB1"   => SB1a,
    "GA2"   => GA2,
    "GA1"   => GA1,
    "GM1b"  => GM1b,
    "GM1α"  => GM1α,
    "GD1c"  => GD1c,
    "GD1α"  => GD1α,
    "GD1e"  => GD1α,
    "GM3"   => GM3,
    "GM2"   => GM2,
    "GM1a"  => GM1a,
    "GD1a"  => GD1a,
    "GD1aα" => GD1aα,
    "GD1aa" => GD1aα,
    "GT1a"  => GT1a,
    "GT1aα" => GT1aα,
    "GT1aa" => GT1aα,
    "GD3"   => GD3,
    "GD2"   => GD2,
    "GD1b"  => GD1b,
    "GT1b"  => GT1b,
    "GT1bα" => GT1bα,
    "GT1ba" => GT1bα,
    "GQ1b"  => GQ1b,
    "GQ1bα" => GQ1bα,
    "GQ1ba" => GQ1bα,
    "GT3"   => GT3,
    "GT2"   => GT2,
    "GT1c"  => GT1c,
    "GQ1c"  => GQ1c,
    "GQ1cα" => GQ1cα,
    "GQ1ca" => GQ1cα,
    "GP1c"  => GP1c,
    "GP1cα" => GP1cα,
    "GP1ca" => GP1cα,
    "SM1"   => SM1,
    "GM1"   => GM1,
    "GD1"   => GD1,
    "GT1"   => GT1,
    "GQ1"   => GQ1,
    "GP1"   => GP1,
    # "GM1?"  => GM1,
    # "GD1?"  => GD1,
    # "GT1?"  => GT1,
    # "GQ1?"  => GQ1,
    # "GP1?"  => GP1,

    "Gb3"   => Gb3,
    "Gb4"   => Gb4,
    "Gb5"   => Gb5,
    "Cytolipin K"   => Gb4,
    "SSEA3 antigen" => Gb5,
    "SSEA4 antigen" => SSEA4Antigen,
    "SSEA4a6 antigen"   => SSEA4α6Antigen,
    "SSEA4α6 antigen"   => SSEA4α6Antigen,
    "Type IV H antigen" => TypeIVHAntigen,
    "Type IV A antigen" => TypeIVAAntigen,
    "Forssman antigen"  => ForssmanAntigen,
    "Forssman"  => ForssmanAntigen,
    "Para-Forssman x3b" => ParaForssmanAntigen,
    "Branched Forssman" => BranchedForssmanAntigen,
    "Globo-Lex-9"   => GloboLex9,
    "iGb3"  => iGb3,
    "iGb4"  => iGb4,
    "iGb5"  => iGb5,
    "Forssman-like iGb4" => ForssmanlikeiGb4,
    "Cytolipin R"   => iGb4,
    "Lc3"   => Lc3,
    "Lc4"   => Lc4,
    "Lea"   => Lea,
    "Leb"   => Leb,
    "Type I H antigen"  => TypeIHAntigen,
    "Type I A antigen"  => TypeIAAntigen,
    "Aleb"  => Aleb,
    "Type I B antigen"  => TypeIBAntigen,
    "Bleb"  => Bleb,
    "Lex-A-9"   => LexA9,
    "Lex-A" => LexA,
    "Ley-A-9"   => LeyA9,
    "Ley-A" => LeyA,
    "LM1"   => LM1,
    "sLea"  => sLea,
    # dsLea
    "nLc4"  => nLc4,
    "nLc5"  => nLc5,
    "nLc6"  => nLc6,
    "I antigen"     => IAntigen,
    "Type II H-7 antigen" => TypeIIH7Antigen,
    "Type II A-8 antigen" => TypeIIA8Antigen,
    "Type II B-8 antigen" => TypeIIB8Antigen,
    "nLM1"  => nLM1,
    "Type II H antigen"   => TypeIIHAntigen,
    "Type II A antigen" => TypeIIAAntigen,
    "Type II B antigen" => TypeIIBAntigen,
    "i antigen"     => nLc6,
    "Lex"   => Lex,
    "Ley"   => Ley,
    "Lex-5" => Lex5,
    "Ley-6" => Ley6,
    "SSEA1 antigen" => Lex5,
    "Lex-7" => Lex7,
    "Lex-X" => LexX,
    "Lex-X-8"   => LexX8,
    "Dimeric Lex-8" => DimericLex8,
    "Ley-8" => Ley8,
    "Ley-X" => LeyX,
    "Ley-X-9"   => LeyX9,
    "Trimeric Ley-9"    => TrimericLey9,
    "Lex-9" => Lex9,
    "Lex-X-10"  => LexX10,
    "Dimeric Lex-10"    => DimericLex10,
    "Lex-X-X"   => LexXX,
    "Lex-X-X-11"=> LexXX11,
    "Trimeric Lex-11"   => TrimericLex11, 
    "Lea-X" => LeaX,
    "Leb-X" => LebX,
    "sLex"  => sLex,
    "sLex-6"=> sLex,
    "sLea-X"=> sLeaX,
    "P antigen"     => PAntigen
)

const MONOSACCHRIDE = collect(keys(MONO_STRUCT))

glycanabbr(::GM4) = "GM4"
glycanabbr(::SM4) = "SM4"
glycanabbr(::Lac) = "Lac"
glycanabbr(::T) where {T <: Ganglioseries} = string(T)
glycanabbr(::T) where {T <: Globoseries} = string(T)
glycanabbr(::T) where {T <: Isogloboseries} = string(T)
glycanabbr(::T) where {T <: Lactoseries} = string(T)
glycanabbr(::T) where {T <: Neolactoseries} = string(T)
glycanabbr(x::T) where {T <: SM1} = isomername(glycanabbr.(x.isomer), "SM1", 2)
glycanabbr(x::T) where {T <: GM1} = isomername(glycanabbr.(x.isomer), "GM1", 3)
glycanabbr(x::T) where {T <: GD1} = isomername(glycanabbr.(x.isomer), "GD1", 5)
glycanabbr(x::T) where {T <: GT1} = isomername(glycanabbr.(x.isomer), "GT1", 5)
glycanabbr(x::T) where {T <: GQ1} = isomername(glycanabbr.(x.isomer), "GQ1", 4)
glycanabbr(x::T) where {T <: GP1} = isomername(glycanabbr.(x.isomer), "GP1", 2)
glycanabbr(::T) where T = string(T)
glycanabbr(::SSEA4Antigen) = "SSEA4 antigen"
glycanabbr(::SSEA4α6Antigen) = "SSEA4a6 antigen"
glycanabbr(::TypeIVHAntigen) = "Type IV H antigen"
glycanabbr(::TypeIVAAntigen) = "Type IV A antigen"
glycanabbr(::TypeIVBAntigen) = "Type IV B antigen"
glycanabbr(::ForssmanAntigen) = "Forssman"
glycanabbr(::ParaForssmanAntigen) = "Para-Forssman x3b"
glycanabbr(::BranchedForssmanAntigen) = "Branched Forssman"
glycanabbr(::GloboLex9) = "Globo-Lex-9"
glycanabbr(::ForssmanlikeiGb4) = "Forssman-like iGb4"
glycanabbr(::TypeIHAntigen) = "Type I H antigen"
glycanabbr(::TypeIAAntigen) = "Type I A antigen"
glycanabbr(::TypeIBAntigen) = "Type I B antigen"
glycanabbr(::LexA) = "Lex-A"
glycanabbr(::LeyA) = "Ley-A"
glycanabbr(::IAntigen) = "I antigen"
glycanabbr(::PAntigen) = "P antigen"
glycanabbr(::TypeIIH7Antigen) = "Type II H-7 antigen"
glycanabbr(::TypeIIA8Antigen) = "Type II A-8 antigen"
glycanabbr(::TypeIIB8Antigen) = "Type II B-8 antigen"
glycanabbr(::TypeIIHAntigen) = "Type II H antigen"
glycanabbr(::TypeIIAAntigen) = "Type II A antigen"
glycanabbr(::TypeIIBAntigen) = "Type II B antigen"
glycanabbr(::Lex7) = "Lex-7"
glycanabbr(::LexX) = "Lex-X"
glycanabbr(::Ley8) = "Ley-8"
glycanabbr(::LeyX) = "Ley-X"
glycanabbr(::Lex9) = "Lex-9"
glycanabbr(::LexX10) = "Lex-X-10"
glycanabbr(::LexXX) = "Lex-X-X"
glycanabbr(::LeaX) = "Lea-X"
glycanabbr(::LebX) = "Leb-X"
glycanabbr(::sLeaX) = "sLea-X"

glycanname(x::GM4) = string("Ganglioside ", glycanabbr(x))
glycanname(::SM4) = "Sulfatide"
glycanname(::Lac) = "Lactose"
glycanname(x::T) where {T <: Ganglioseries} = string("Ganglioside ", glycanabbr(x))
glycanname(x::T) where {T <: Sulfoganglioseries} = string("Sulfoganglioside ", glycanabbr(x))
# glycanname(x::T) where {T <: Globoseries} = string(T)
# glycanname(x::T) where {T <: Isogloboseries} = string(T)
# glycanname(x::T) where {T <: Lactoseries} = string(T)
# glycanname(x::T) where {T <: Neolactoseries} = string(T)
glycanname(x::T) where {T <: SM1} = string("Sulfoganglioseries ", glycanabbr(x))
glycanname(x::T) where {T <: GM1} = string("Ganglioside ", glycanabbr(x))
glycanname(x::T) where {T <: GD1} = string("Ganglioside ", glycanabbr(x))
glycanname(x::T) where {T <: GT1} = string("Ganglioside ", glycanabbr(x))
glycanname(x::T) where {T <: GQ1} = string("Ganglioside ", glycanabbr(x))
glycanname(x::T) where {T <: GP1} = string("Ganglioside ", glycanabbr(x))
glycanname(::Gb3) = "Globoside 3"
glycanname(::CytolipinK) = "Globoside 4 (Cytolipin K)"
glycanname(::SSEA3Antigen) = "Globoside 5 (SSEA3 antigen)"
glycanname(::SSEA4Antigen) = "SSEA4 antigen"
glycanname(::SSEA4α6Antigen) = "SSEA4(α1-6) antigen"
glycanname(::TypeIVHAntigen) = "Type IV H antigen"
glycanname(::TypeIVAAntigen) = "Type IV A antigen"
glycanname(::TypeIVBAntigen) = "Type IV B antigen"
glycanname(::ForssmanAntigen) = "Forssman antigen"
glycanname(::ParaForssmanAntigen) = "Para-Forssman x3b antigen"
glycanname(::BranchedForssmanAntigen) = "Branched Forssman antigen"
glycanname(::GloboLex9) = "Globo-Lex-9"
glycanname(::iGb3) = "Isogloboside 3"
glycanname(::iGb4) = "Isogloboside 4 (Cytolipin R)"
glycanname(::iGb5) = "Isogloboside 5"
glycanname(::ForssmanlikeiGb4) = "Forssman-like iGb4"
glycanname(::Lc3) = "Lacto series 3"
glycanname(::Lc4) = "Lacto series 4"
glycanname(::Lea) = "Lewis a"
glycanname(::Leb) = "Lewis b"
glycanname(::TypeIHAntigen) = "Type I H antigen"
glycanname(::TypeIAAntigen) = "Type I A antigen"
glycanname(::TypeIBAntigen) = "Type I B antigen"
glycanname(::Aleb) = "Type I A-Lewis b"
glycanname(::Bleb) = "Type I B-Lewis b"
glycanname(::LexA) = "Lewis x-a"
glycanname(::LeyA) = "Lewis y-a"
glycanname(::sLea) = "Sialyl Lewis a"
glycanname(::dsLea) = "Disialyl Lewis a"
glycanname(::LM1) = "Sialyl Lacto series 4"
glycanname(::nLc4) = "Neolacto series 4"
glycanname(::nLc5) = "Neolacto series 5"
glycanname(::nLc6) = "Neolacto series 6 (i antigen)"
glycanname(::IAntigen) = "I antigen"
glycanname(::PAntigen) = "P antigen"
glycanname(::TypeIIH7Antigen) = "Type II H-7 antigen"
glycanname(::TypeIIA8Antigen) = "Type II A-8 antigen"
glycanname(::TypeIIB8Antigen) = "Type II B-8 antigen"
glycanname(::nLM1) = "Sialyl Neolacto series 4"
glycanname(::TypeIIHAntigen) = "Type II H antigen"
glycanname(::TypeIIAAntigen) = "Type II A antigen"
glycanname(::TypeIIBAntigen) = "Type II B antigen"
glycanname(::Lex) = "Lewis x (SSEA1 antigen)"
glycanname(::Ley) = "Lewis y"
glycanname(::Lex7) = "Lewis x-7"
glycanname(::LexX) = "Lewis x-x"
glycanname(::Ley8) = "Lewis y-8"
glycanname(::LeyX) = "Lewis y-x"
glycanname(::Lex9) = "Lewis x-9"
glycanname(::LexX10) = "Lewis x-x-10"
glycanname(::LexXX) = "Lewis x-x-x"
glycanname(::LeaX) = "Lewis a-x"
glycanname(::LebX) = "Lewis b-x"
glycanname(::sLex) = "Sialyl Lewis x"
glycanname(::sLeaX) = "Sialyl Lewis a-x"

function glysubname(x::Vector{<: Pair{<: G, UInt8}}) where G 
    ms = map(chemicalname ∘ first, x)
    ns = last.(x)
    join([n == 1 ? m : endswith(m, r"\d") ? string("(", m ,")", n) : string(m, n) for (m, n) in zip(ms, ns)], "/")
end
function glysubabbr(x::Vector{<: Pair{<: G, UInt8}}) where G 
    ms = map(chemicalabbr ∘ first, x)
    ns = last.(x)
    join([n == 1 ? m : endswith(m, r"\d") ? string("(", m ,")", n) : string(m, n) for (m, n) in zip(ms, ns)], "/")
end
function glysubname(x::Vector{<: Pair{<: Tuple{UInt8, <: Vector}}})  
    ms = map(chemicalname ∘ last, x)
    ns = map(x) do n 
        p, l = first(n)
        if isempty(l)
            string(Int(p))
        else
            string(Int(p), "[", join([string(Int(i)) for i in l], ","), "]")
        end
    end
    join([string(n, m) for (m, n) in zip(ms, ns)], "/")
end
function glysubabbr(x::Vector{<: Pair{<: Tuple{UInt8, <: Vector}}})  
    ms = map(chemicalabbr ∘ last, x)
    ns = map(x) do n 
        p, l = first(n)
        if isempty(l)
            string(Int(p))
        else
            string(Int(p), "[", join([string(Int(i)) for i in l], ","), "]")
        end
    end
    join([string(n, m) for (m, n) in zip(ms, ns)], "/")
end

function push_new_glycan!(cs, ls, x, c, d)
    mono, mp, mi, mj = x
    push!(cs, c)
    i = isempty(mi) ? first(d) : parse(UInt8, mi)
    lf = (isempty(mp) || mp == "(") ? Anomerposition(i) : (mp == "α" || mp == "(a") ? Alphaposition(i) : (mp == "β" || mp == "(b") ? Betaposition(i) : throw(ArgumentError("Invalid glycosidic linkage."))
    rt = isempty(mj) ? Linkageposition(last(d)) : Linkageposition(parse(UInt8, mj))
    push!(ls, lf => rt)
    isempty(mj)
end

"""
    parse_saccharide(s::AbstractString)

Parse saccharide abbreviation into
* `Monosaccharide`: a single monosaccharide abbreviation, e.g. `Glc`.
* `Glycan`: multiple connected monosaccharides abbreviation connected, e.g. `Neu5Acα-3Galβ-3[Neu5Acα-6]GalNAcβ-4[Neu5Acα-8Neu5Acα-3]Galβ-4Glc`.
* `SeriesGlycan`: named glycan, e.g. `GQ1bα(3[3]NeuGc/4[3]Neu?(Ac2)/4[4,6]NeuGc8Ac)`.
* `GlyComp`: glycan composition, e.g. `HexNAcHex3`.
"""
function parse_saccharide(s::AbstractString)
    try 
        parse_monosaccharide(s)
    catch
        try 
            parse_glycan(s)
        catch
            try
                parse_glycomp(s)
            catch
                throw(ArgumentError("Invalid sacchride representation."))
            end
        end
    end
end

function parse_glycan(s::AbstractString)
    try 
        parse_series_glycan(s)
    catch 
        parse_generic_glycan(s)
    end
end

function parse_series_glycan(s::AbstractString)
    m, s = match(r"([^\(\)]*)(?:\((.*)\))?", s)
    m = replace(m, "alpha" => "α")
    if haskey(GLYCAN_STRUCT, m)
        m = GLYCAN_STRUCT[m]()
    else
        m, p = split(m, "?"; limit = 2)
        iso = split(p, r"\s*/\s*")
        sort!(iso)
        m = GLYCAN_STRUCT[m](ntuple(i -> eval(Meta.parse(string(m, iso[i])))(), length(iso)))
    end
    if isnothing(s)
        sub = nothing
    else
        sub = [parse_glycansub(x) for x in split(s, r"\s*/\s*")]
        sub = convert(Vector{Pair{promote_type(typeof.(first.(sub))...), promote_type(typeof.(last.(sub))...)}}, sub)
    end
    g = SeriesGlycan(m, sub)
    generic_glycan(g)
    g
end

function parse_generic_glycan(s::AbstractString)
    p = 1
    cs = Saccharide[]
    ls = Pair{AbstractAnomerposition, Linkageposition}[]
    linktocheck = nothing
    trace = false
    @inbounds for m in eachmatch(r"(\[[^\]\[]*(?:(?1)[^\]\[]*)*-\d*\)?\])", s)
        n = m.match.offset
        if n > p
            rs = collect(eachmatch(r"((?:[DL]+-)*[^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:n]))
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
            end
            linktocheck = trace ? lastindex(ls) : nothing
        end
        push!(cs, parse_glycan(first(match(r"^\[?(.*?)\]?$", m.match))))
        push!(ls, last(last(cs).linkage))
        p = n + m.match.ncodeunits + 1
    end
    li = p
    rs = collect(eachmatch(r"((?:[DL]+-)*[^-]*?[^\d])(\(?(?:\(a)?(?:\(b)?[αβ]?)(\d*)-(\d*)\)?", s[p:end]))
    if !isempty(rs) && first(last(rs)) in ["L", "D", "DL"]
        if firstindex(rs) == lastindex(rs)
            rs = empty(rs)
        else
            rs = rs[begin:end - 1]
        end
    end
    if !isempty(rs)
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
    end
    linktocheck = trace ? lastindex(ls) : nothing
    if li < lastindex(s)
        push!(cs, parse_monosaccharide(first(match(r"^\[?(.*?)\]?$", s[li:end]))))
        if !isnothing(linktocheck)
            ls[linktocheck] = first(ls[linktocheck]) => Linkageposition(dehydrogenposition(last(cs)))
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
    l, s = match(r"^((?:[DL]+-))*(.*)", string(s))
    @inbounds for i in Iterators.reverse(eachindex(s))
        mono = s[begin:i]
        if haskey(MONO_STRUCT, mono)
            if i == lastindex(s)
                M = MONO_STRUCT[mono]
                if implicit_config(M)
                    DL = NoDLForm
                else
                    DL = l == "L-" ? LForm : l == "D-" ? DForm : l == "DL-" ? DLForm : nothing
                end
                return isnothing(DL) ? M() : M(; DL)
            else
                m = s[nextind(s, i):end]
                break
            end
        end
    end
    isnothing(mono) && throw(ArgumentError("Invalid monosacchride, \"$s\""))
    subs = split_snfg_sub(m)
    M = MONO_STRUCT[mono]
    if implicit_config(M)
        DL = NoDLForm
    else
        DL = l == "L-" ? LForm : l == "D-" ? DForm : l == "DL-" ? DLForm : nothing
    end
    isempty(subs) && return isnothing(DL) ? M() : M(; DL)
    allsub = mapreduce(vcat, subs) do sub
        p, a = match(r"^([\d,?]*)(.*)$", sub)
        if startswith(a, "(") && endswith(a, r"\d\)") && p == "?"
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
        substituent = [SNFG_SUB[x] => UInt8(count(==(x), sk)) for x in a]
    else
        sort!(allsub; by = reverse)
        substituent = [UInt8(first(x)) => SNFG_SUB[last(x)] for x in allsub]
    end
    isnothing(DL) ? M(substituent) : M(substituent; DL)
end

function get_monosaccharide(x)
    m, mp, mi, mj = x
    parse_monosaccharide(first(match(r"^\[?(.*?)\]?$", m)))
end

function parse_glycansub(mod)
    if startswith(mod, r"\d")
        # position 
        # config
        m = match(r"^(\d*)(?:\[(.*)\])?(.*)$", mod)
        isnothing(m) && throw(ArgumentError("Invalid functional group modification, \"$mod\""))
        p, l, s = m 
        if isnothing(l)
            pos = (parse(UInt8, p), UInt8[])
        else
            pos = (parse(UInt8, p), [parse(UInt8, x) for x in split(l, ",")])
        end
        return pos => parse_monosaccharide(s)
    else
        for (k, v) in SNFG_SUB
            m = match(Regex(string("^\\(?", k, "\\)?(\\d*)\$")), mod)
            isnothing(m) ? continue : return v => (isempty(first(m)) ? 0x01 : parse(UInt8, first(m)))
        end
        e = lastindex(mod)
        for i in Iterators.reverse(eachindex(mod))
            if isnothing(match(r"\d", mod[i:i]))
                e = i 
                break
            else
                e = i 
            end
        end
        m = mod[firstindex(mod):e]
        n = e == lastindex(mod) ? "" : mod[nextind(mod, e):lastindex(mod)]
        return parse_monosaccharide(m) => (isempty(n) ? 0x01 : parse(UInt8, n))
    end
end

function add_ldform(mono::T, head) where {L, T <: Monosaccharide{L}} 
    L <: LForm ? string("L-", head) : L <: DForm ? string("D-", head) : head 
end

function add_sub(mono::T, head, sub = nothing) where {T <: Monosaccharide}
    if (isnothing(mono.substituent) || isempty(mono.substituent)) 
        isnothing(sub) && return head
        sub = [m => snfg_abbr(n) for (m, n) in sub]
    elseif mono.substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        if isnothing(sub)
            return string(head, join([n == 1 ? string("?", snfg_abbr(m)) : string("?(", snfg_abbr(m), Int(n), ")") for (m, n) in mono.substituent], ""))
        else
            sub = [snfg_abbr(n) => m for (m, n) in sub]
            for (m, n) in mono.substituent
                for _ in 1:n
                    push!(sub, 0x00 => snfg_abbr(m))
                end
            end
        end
    else
        sub = [m => snfg_abbr(n) for (m, n) in vcat(sub, mono.substituent)]
    end
    sort!(sub)
    ps = ""
    pp = String[]
    lm = "" 
    for (p, m) in sub
        if isnothing(lm)
            lm = m
            push!(pp, p == 0 ? "?" : string(Int(p)))
        elseif m == lm
            push!(pp, p == 0 ? "?" : string(Int(p)))
        else
            ps = string(ps, join(pp, ","), lm)
            lm = m 
            pp = String[p == 0 ? "?" : string(Int(p))]
        end
    end
    ps = string(ps, join(pp, ","), lm)
    string(head, ps)
end

function add_name(mono::T, head) where {T <: Monosaccharide}
    (isnothing(mono.substituent) || isempty(mono.substituent)) && return head
    if mono.substituent isa Vector{<: Pair{<: FunctionalGroup, <: UInt8}}
        string(join([n == 1 ? string("?-", lowercase(chemicalname(m; group = false))) : string(join(repeat(["?"], Int(n)), ","), "-", lowercase(chemicalname(m; group = false))) for (m, n) in mono.substituent], "-"), "-", string(T.name.name))
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
                ps = string(ps, join(pp, ","), "-", lowercase(chemicalname(lm; group = false)), "-")
                lm = m 
                pp = String[p == 0 ? "?" : string(Int(p))]
            end
        end
        ps = string(ps, join(pp, ","), "-", lowercase(chemicalname(lm; group = false)), "-")
        string(ps, head)
    end
end

function split_glycomp(glycan::AbstractString)
    s = Int[firstindex(glycan)]
    e = Int[firstindex(glycan)]
    dp = 0
    endp = true
    for i in eachindex(glycan)
        @inbounds x = glycan[i]
        if (dp == 0 && !endp && isnothing(match(r"\d", string(x)))) 
            push!(s, i)
            push!(e, i)
            endp = true
            continue
        end
        if x == '('
            dp += 1
        elseif x == ')'
            dp -= 1
        end
        if dp == 1 && endp
            push!(s, i)
            push!(e, i)
            endp = false
            continue
        end
        e[end] = i
    end
    @inbounds [glycan[i:j] for (i, j) in zip(s, e)]
end

function parse_glycomp(s::AbstractString)
    v = mapreduce(vcat, split_glycomp(s)) do x 
        result = Pair{Monosaccharide, UInt8}[]
        iter_glycomp!(result, x)
    end
    d = Dict{eltype(v).parameters...}()
    for (m, n) in v 
        get!(d, m, 0x00)
        d[m] += n 
    end
    v = collect(pairs(d))
    GlyComp(sort!(v; by = x -> (first(x) isa AbstractFunctionalGroup, (repr ∘ first)(x))))
end

function iter_glycomp!(result::Vector, s::AbstractString)
    if startswith(s, '(')
        m, n = match(r"^\((.*)\)(\d)*$", s)
        mono = parse_monosaccharide(m)
        push!(result, mono => isnothing(n) ? 0x01 : parse(UInt8, n))
        return result
    end
    @inbounds for i in Iterators.reverse(eachindex(s))
        if startswith(s, "_")
            result = convert(Vector{Pair{Any, UInt8}}, result)
            return iter_snfg!(result, s[begin + 1:end])
        end
        try
            mono = parse_monosaccharide(s[begin:i])
            if i == lastindex(s)
                push!(result, mono => 0x01)
                return result
            else
                s = string(s[nextind(s, i):end])
                n = match(r"^\d+", s)
                if isnothing(n)
                    push!(result, mono => 0x01)
                    return iter_glycomp!(result, s)
                else
                    i = n.match.offset + n.match.ncodeunits
                    n = parse(UInt8, n.match)
                    push!(result, mono => n)
                    return i == lastindex(s) ? result : iter_glycomp!(result, s[nextind(s, i):end])
                end
            end
        catch
            continue
        end
    end
    throw(ArgumentError("Invalid monosaccharide(s), \"$s\""))
end

function iter_snfg!(result::Vector, s::AbstractString)
    @inbounds for i in Iterators.reverse(eachindex(s))
        try
            mono = SNFG_SUB[s[begin:i]]
            if i == lastindex(s)
                push!(result, mono => 0x01)
                return result
            else
                s = string(s[nextind(s, i):end])
                n = match(r"^\d+", s)
                if isnothing(n)
                    push!(result, mono => 0x01)
                    return iter_snfg!(result, s)
                else
                    i = n.match.offset + n.match.ncodeunits
                    n = parse(UInt8, n.match)
                    push!(result, mono => n)
                    return i == lastindex(s) ? result : iter_snfg!(result, s[nextind(s, i):end])
                end
            end
        catch
            continue
        end
    end
    throw(ArgumentError("Invalid functionsl group(s), \"$s\""))
end

function isomername(isomer, c, n)
    if length(isomer) >= n
        c
    else
        string(c, "?", join([replace(i, c => "") for i in isomer], "/"))
    end
end
