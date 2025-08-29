const PROTEIN_3LETTER_AA = Dict{String, Type}(
    "Ala"   => Alanine,
    "Arg"   => Arginine,
    "Asn"   => Aspargine,
    "Asp"   => Aspartate,
    "Cys"   => Cysteine,
    "Gln"   => Glutamine,
    "Glu"   => Glutamate,
    "Gly"   => Glycine,
    "His"   => Histidine,
    "Ile"   => Isoleucine,
    "Leu"   => Leucine,
    "Lys"   => Lysine,
    "Met"   => Methionine,
    "Phe"   => Phenylalanine,
    "Pro"   => Proline,
    "Ser"   => Serine,
    "Thr"   => Threonine,
    "Trp"   => Tryptophan,
    "Tyr"   => Tyrosine,
    "Val"   => Valine,
    "Sec"   => Selenocysteine,
    "Pyl"   => Pyrrolysine,
    "Orn"   => Ornithine    
)
const PROTEIN_3LETTER_FG = Dict{String, Type}(
    "Ala"   => Alanyl,
    "Arg"   => Arginyl,
    "Asn"   => Asparginyl,
    "Asp"   => Aspartyl,
    "Cys"   => Cysteinyl,
    "Gln"   => Glutaminyl,
    "Glu"   => Glutamyl,
    "Gly"   => Glycyl,
    "His"   => Histidinyl,
    "Ile"   => Isoleucyl,
    "Leu"   => Leucyl,
    "Lys"   => Lysyl,
    "Met"   => Methionyl,
    "Phe"   => Phenylalanyl,
    "Pro"   => Prolyl,
    "Ser"   => Seryl,
    "Thr"   => Threonyl,
    "Trp"   => Tryptophanyl,
    "Tyr"   => Tyrosyl,
    "Val"   => Valyl,
    "Sec"   => Selenocysteinyl,
    "Pyl"   => Pyrrolysyl,
    "Orn"   => Ornithyl    
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

letter1_abbr(::Alanine) = "A"
letter1_abbr(::Arginine) = "R"
letter1_abbr(::Aspargine) = "N"
letter1_abbr(::Aspartate) = "D"
letter1_abbr(::Cysteine) = "C"
letter1_abbr(::Glutamine) = "Q"
letter1_abbr(::Glutamate) = "E"
letter1_abbr(::Glycine) = "G"
letter1_abbr(::Histidine) = "H"
letter1_abbr(::Isoleucine) = "I"
letter1_abbr(::Leucine) = "L"
letter1_abbr(::Lysine) = "K"
letter1_abbr(::Methionine) = "M"
letter1_abbr(::Phenylalanine) = "F"
letter1_abbr(::Proline) = "P"
letter1_abbr(::Serine) = "S"
letter1_abbr(::Threonine) = "T"
letter1_abbr(::Tryptophan) = "W"
letter1_abbr(::Tyrosine) = "Y"
letter1_abbr(::Valine) = "V"
letter1_abbr(::Selenocysteine) = "U"
letter1_abbr(::Pyrrolysine) = "O"
letter1_abbr(::Ornithine) = "Orn"

# isotope?
function parse_aa(s::AbstractString)
    l, r, s = match(r"([DL]-)*(\([RS]\)-)*([A-Z][a-z]*)", s)
    A = length(s) == 3 ? PROTEIN_3LETTER_AA[string(s)] : PROTEIN_3LETTER_AA[PROTEIN_1LETTER_3LETTER[string(s)]]
    l = isnothing(l) ? r : l 
    _parse_aa(A, l)
end

function parse_aa_fg(s::AbstractString)
    l, r, s = match(r"([DL]-)*(\([RS]\)-)*([A-Z][a-z]*)", s)
    A = length(s) == 3 ? PROTEIN_3LETTER_FG[string(s)] : PROTEIN_3LETTER_FG[PROTEIN_1LETTER_3LETTER[string(s)]]
    l = isnothing(l) ? r : l 
    _parse_aa(A, l)
end

function parse_aa3(s::AbstractString)
    l, r, s = match(r"([DL]-)*(\([RS]\)-)*([A-Z][a-z]*)", s)
    A = PROTEIN_3LETTER_AA[string(s)]
    l = isnothing(l) ? r : l 
    _parse_aa(A, l)
end

function parse_aa3_fg(s::AbstractString)
    l, r, s = match(r"([DL]-)*(\([RS]\)-)*([A-Z][a-z]*)", s)
    A = PROTEIN_3LETTER_AA_FG[string(s)]
    l = isnothing(l) ? r : l 
    _parse_aa(A, l)
end

function _parse_aa(A, l)
    if isnothing(l) 
        chiralchemical(A, DLForm)
    elseif l == "L-"
        A{LForm}()
    elseif l == "D-"
        A{DForm}()
    elseif l == "(S)-"
        chiralchemical(A, SChirality)
    elseif l == "(R)-"
        chiralchemical(A, RChirality)
    else
        @warn "Invalid Fisher projection or RS chirality"
        chiralchemical(A, DLForm)
    end
end

chiralchemical(::Type{A}, x::Type{I}) where {I <: GlyceraldehydeSystem, A <: Union{Glycine, Glycyl}} = A{NoDLForm}()
chiralchemical(::Type{A}, x::Type{I}) where {I <: RSSystem, A <: Union{Glycine, Glycyl}} = A{NoDLForm}()
chiralchemical(::Type{A}, ::Type{I}) where {I <: GlyceraldehydeSystem, A <: Union{<: αAminoAcid, <: FunctionalGroup{<: αAminoAcid}}} = A{I}()
chiralchemical(::Type{A}, ::Type{SChirality}) where {A <: Union{<: αAminoAcid, <: FunctionalGroup{<: αAminoAcid}}} = A{LForm}()
chiralchemical(::Type{A}, ::Type{RChirality}) where {A <: Union{<: αAminoAcid, <: FunctionalGroup{<: αAminoAcid}}} = A{DForm}()
chiralchemical(::Type{A}, ::Type{SChirality}) where {A <: Union{Cysteine, Cysteinyl, Selenocysteine, Selenocysteinyl}} = A{DForm}()
chiralchemical(::Type{A}, ::Type{RChirality}) where {A <: Union{Cysteine, Cysteinyl, Selenocysteine, Selenocysteinyl}} = A{LForm}()

function getchemicalattr(aa::T, ::Val{:name}; kwargs...) where {L, T <: αAminoAcid{L}}
    s = replace(repr(T), r"\{*[DLForm]*\}*$" => "")
    L <: LForm ? string("L-", s) : L <: DForm ? string("D-", s) : s
end

function getchemicalattr(aa::T, ::Val{:abbreviation}; nletter = 3) where {L, T <: αAminoAcid{L}}
    s = nletter == 3 ? letter3_abbr(aa) : nletter == 1 ? letter1_abbr(aa) : nothing
    L <: LForm ? string("L-", s) : L <: DForm ? string("D-", s) : s
end
getchemicalattr(aa::Peptide, ::Val{:name}; kwargs...) = join([getchemicalattr(x, :name; kwargs...) for x in aa.chain], "-")
getchemicalattr(aa::Peptide, ::Val{:abbreviation}; kwargs...) = join([getchemicalattr(x, :abbreviation; kwargs...) for x in aa.chain], "-")

letter3_abbr(aa::T) where {T <: αAminoAcid} = PROTEIN_1LETTER_3LETTER[letter1_abbr(aa)]
alpha_rs(aa::αAminoAcid{LForm}) = SChirality()
alpha_rs(aa::αAminoAcid{DForm}) = RChirality()
alpha_rs(aa::αAminoAcid{DLForm}) = RSChirality()
alpha_rs(aa::Glycine) = AChirality()
alpha_rs(aa::Cysteine{LForm}) = RChirality()
alpha_rs(aa::Selenocysteine{LForm}) = RChirality()
alpha_rs(aa::Cysteine{DForm}) = SChirality()
alpha_rs(aa::Selenocysteine{DForm}) = SChirality()