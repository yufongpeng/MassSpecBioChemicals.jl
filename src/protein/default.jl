defaultchemicalkwargs(::Type{T}) where {T <: αAminoAcid} =
    (; name = chemicaldefaultname(T), elements = chemicaldefaultelements(T),
       abbreviation = chemicaldefaultabbr(T))

defaultchemicalkwargs(::Type{<:Peptide}, chain) =
    (; name = join(chemicalname.(chain), "-"), elements = _peptideelements(chain),
       abbreviation = join(chemicalabbr.(chain), "-"))

function _peptideelements(chain)
    n = length(chain)
    n == 0 && return Pair{String, Int}[]
    reduce(vcat, (losswaterelements(m, i, i % n) for (i, m) in enumerate(chain)); init = Pair{String, Int}[])
end

chemicaldefaultname(::Type{T}) where {T <: αAminoAcid} = string(nameof(Base.typename(T).wrapper))
chemicaldefaultname(::Type{T}) where {T <: αAminoAcid{DForm}} = string("D-", nameof(Base.typename(T).wrapper))
chemicaldefaultname(::Type{T}) where {T <: αAminoAcid{LForm}} = string("L-", nameof(Base.typename(T).wrapper))
chemicaldefaultabbr(::Type{T}) where {T <: αAminoAcid} = first(k for (k, v) in PROTEIN_3LETTER_AA if v === Base.typename(T).wrapper)
chemicaldefaultabbr(::Type{T}) where {T <: αAminoAcid{DForm}} = string("D-", first(k for (k, v) in PROTEIN_3LETTER_AA if v === Base.typename(T).wrapper))
chemicaldefaultabbr(::Type{T}) where {T <: αAminoAcid{LForm}} = string("L-", first(k for (k, v) in PROTEIN_3LETTER_AA if v === Base.typename(T).wrapper))

chemicaldefaultelements(::Type{<:αAminoAcid}) = Pair{String, Int}[]
chemicaldefaultelements(::Type{<:Alanine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Arginine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 6, "N" => 1, "H" => 1, "C" => 1, "N" => 1, "H" => 1, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Aspargine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Aspartate}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Cysteine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "S" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Glutamate}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 2, "H" => 4, "C" => 1, "O" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Glycine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Histidine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 3, "H" => 3, "N" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Isoleucine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 9, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Leucine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 9, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Lysine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 8, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Methionine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "S" => 1, "C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Phenylalanine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 5, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Proline}) = ["H" => 1, "N" => 1, "C" => 4, "H" => 7, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Serine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Threonine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 2, "H" => 4, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Tryptophan}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 8, "H" => 6, "N" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Tyrosine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 4, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Valine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 7, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Selenocysteine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "Se" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Pyrrolysine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 8, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 4, "H" => 7, "N" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Ornithine}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 6, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:fiveHydroxytryptophan}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 8, "H" => 6, "N" => 1, "O" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:DOPA}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
