const PROTEIN_3LETTER_AA = Dict{String, αAminoAcid}(
    "Ala"   => Alanine(),
    "Arg"   => Arginine(),
    "Asn"   => Aspargine(),
    "Asp"   => Aspartate(),
    "Cys"   => Cysteine(),
    "Gln"   => Glutamine(),
    "Glu"   => Glutamate(),
    "Gly"   => Glycine(),
    "His"   => Histidine(),
    "Ile"   => Isoleucine(),
    "Leu"   => Leucine(),
    "Lys"   => Lysine(),
    "Met"   => Methionine(),
    "Phe"   => Phenylalanine(),
    "Pro"   => Proline(),
    "Ser"   => Serine(),
    "Thr"   => Threonine(),
    "Trp"   => Tryptophan(),
    "Tyr"   => Tyrosine(),
    "Val"   => Valine(),
    "Sec"   => Selenocysteine(),
    "Pyl"   => Pyrrolysine(),
    "Orn"   => Ornithine()    
)
const PROTEIN_3LETTER_FG = Dict{String, FunctionalGroup}(
    "Ala"   => Alanyl(),
    "Arg"   => Arginyl(),
    "Asn"   => Asparginyl(),
    "Asp"   => Aspartyl(),
    "Cys"   => Cysteinyl(),
    "Gln"   => Glutaminyl(),
    "Glu"   => Glutamyl(),
    "Gly"   => Glycyl(),
    "His"   => Histidinyl(),
    "Ile"   => Isoleucyl(),
    "Leu"   => Leucyl(),
    "Lys"   => Lysyl(),
    "Met"   => Methionyl(),
    "Phe"   => Phenylalanyl(),
    "Pro"   => Prolyl(),
    "Ser"   => Seryl(),
    "Thr"   => Threonyl(),
    "Trp"   => Tryptophanyl(),
    "Tyr"   => Tyrosyl(),
    "Val"   => Valyl(),
    "Sec"   => Selenocysteinyl(),
    "Pyl"   => Pyrrolysyl(),
    "Orn"   => Ornithyl()    
)

const PROTEIN_1LETTER_3LETTER = Dict{String, String}(
    "A"   => "Ala",
    "R"   => "Arg",
    "N"   => "Asn",
    "D"   => "Asp",
    "C"   => "Cys",
    "Q"   => "Gln",
    "E"   => "Glu",
    "G"   => "Gly",
    "H"   => "His",
    "I"   => "Ile",
    "L"   => "Leu",
    "K"   => "Lys",
    "M"   => "Met",
    "F"   => "Phe",
    "P"   => "Pro",
    "S"   => "Ser",
    "T"   => "Thr",
    "W"   => "Trp",
    "Y"   => "Tyr",
    "V"   => "Val",
    "U"   => "Sec",
    "O"   => "Pyl",
    "Orn" => "Orn"
)

const PROTEIN_AA_1LETTER = Dict{αAminoAcid, String}(
    Alanine()       => "A",
    Arginine()      => "R",
    Aspargine()     => "N",
    Aspartate()     => "D",
    Cysteine()      => "C",
    Glutamine()     => "Q",
    Glutamate()     => "E",
    Glycine()       => "G",
    Histidine()     => "H",
    Isoleucine()    => "I",
    Leucine()       => "L",
    Lysine()        => "K",
    Methionine()    => "M",
    Phenylalanine() => "F",
    Proline()       => "P",
    Serine()        => "S",
    Threonine()     => "T",
    Tryptophan()    => "W",
    Tyrosine()      => "Y",
    Valine()        => "V",
    Selenocysteine() => "U",
    Pyrrolysine()    => "O",
    Ornithine()     => "Orn"
)

# isotope?
function parse_aa(s::AbstractString)
    length(s) == 3 ? PROTEIN_3LETTER_AA[string(s)] : PROTEIN_3LETTER_AA[PROTEIN_1LETTER_3LETTER[string(s)]]
end

function parse_aa_fg(s::AbstractString)
    length(s) == 3 ? PROTEIN_3LETTER_FG[string(s)] : PROTEIN_3LETTER_FG[PROTEIN_1LETTER_3LETTER[string(s)]]
end

function parse_aa3(s::AbstractString)
    PROTEIN_3LETTER_AA[string(s)]
end

function getchemicalattr(aa::T, ::Val{:name}; kwargs...) where {T <: αAminoAcid}
    replace(repr(T), "()" => "")
end

getchemicalattr(aa::T, ::Val{:abbreviation}; nletter = 3) where {T <: αAminoAcid} = nletter == 3 ? letter3_abbr(aa) : nletter == 1 ? letter1_abbr(aa) : nothing

function letter3_abbr(aa::T) where {T <: αAminoAcid}
    PROTEIN_1LETTER_3LETTER[PROTEIN_AA_1LETTER[aa]]
end

function letter1_abbr(aa::T) where {T <: αAminoAcid}
    PROTEIN_AA_1LETTER[aa]
end