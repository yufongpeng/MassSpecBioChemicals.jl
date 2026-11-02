defaultchemicalkwargs(::Type{T}) where {T <: Metabolite} =
    (; name = chemicaldefaultname(T), elements = chemicaldefaultelements(T),
       abbreviation = chemicaldefaultabbr(T))
       
chemicaldefaultname(::Type{T}) where {T <: Metabolite} = string(nameof(Base.typename(T).wrapper))
chemicaldefaultelements(::Type{<:Metabolite}) = Pair{String, Int}[]
chemicaldefaultabbr(::Type{T}) where {T <: Metabolite} = chemicaldefaultname(T)
# chemicaldefaultname(::Type{Ethanolamine}) = "Ethanolamine"
# chemicaldefaultname(::Type{Nmethylethanolamine}) = "N-Methylethanolamine"
# chemicaldefaultname(::Type{NNdimethylethanolamine}) = "Dimethylethanolamine"
# chemicaldefaultname(::Type{GABA}) = "GABA"
# chemicaldefaultname(::Type{Dopamine}) = "Dopamine"
# chemicaldefaultname(::Type{Taurine}) = "Taurine"
# chemicaldefaultname(::Type{<:Carnitine}) = "Carnitine"
chemicaldefaultname(::Type{Carnitine{LForm}}) = "L-Carnitine"
chemicaldefaultname(::Type{Carnitine{DForm}}) = "D-Carnitine"
# chemicaldefaultname(::Type{Choline}) = "Choline"
# chemicaldefaultname(::Type{CoA}) = "CoA"
# chemicaldefaultname(::Type{GlycolicAcid}) = "Glycolic acid"
# chemicaldefaultname(::Type{<:LacticAcid}) = "Lactic acid"
# chemicaldefaultname(::Type{LacticAcid{LForm}}) = "L-Lactic acid"
# chemicaldefaultname(::Type{LacticAcid{DForm}}) = "D-Lactic acid"
# chemicaldefaultname(::Type{PyruvicAcid}) = "Pyruvic acid"
# chemicaldefaultname(::Type{<:Glycerol}) = "Glycerol"
chemicaldefaultname(::Type{Glycerol{RChirality}}) = "(R)-Glycerol"
chemicaldefaultname(::Type{Glycerol{SChirality}}) = "(S)-Glycerol"
# chemicaldefaultname(::Type{<:Serotonin}) = "Serotonin"
chemicaldefaultname(::Type{Serotonin{RChirality}}) = "(R)-Serotonin"
chemicaldefaultname(::Type{Serotonin{SChirality}}) = "(S)-Serotonin"
# chemicaldefaultname(::Type{<:Melatonin}) = "Melatonin"
chemicaldefaultname(::Type{Melatonin{RChirality}}) = "(R)-Melatonin"
chemicaldefaultname(::Type{Melatonin{SChirality}}) = "(S)-Melatonin"

chemicaldefaultelements(::Type{Ethanolamine}) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{Nmethylethanolamine}) = ["H" => 1, "C" => 1, "H" => 3, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{NNdimethylethanolamine}) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{GABA}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{Dopamine}) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "C" => 6, "H" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{Taurine}) = ["H" => 2, "N" => 1, "C" => 2, "H" => 4, "S" => 1, "O" => 2, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Carnitine}) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 3, "H" => 5, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1]
chemicaldefaultelements(::Type{Choline}) = ["C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "N" => 1, "C" => 2, "H" => 3, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{CoA}) = ["H" => 1, "S" => 1, "C" => 2, "H" => 4, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 2, "H" => 4, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 5, "H" => 9, "O" => 1, "H" => 1, "O" => 1, "P" => 2, "O" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 5, "H" => 8, "O" => 3, "O" => 1, "P" => 1, "O" => 1, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 5, "H" => 2, "N" => 4, "N" => 1, "H" => 2]
chemicaldefaultelements(::Type{GlycolicAcid}) = ["H" => 1, "O" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:LacticAcid}) = ["H" => 1, "O" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{PyruvicAcid}) = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Glycerol}) = ["H" => 1, "O" => 1, "C" => 3, "H" => 5, "O" => 1, "H" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Serotonin}) = ["H" => 2, "N" => 1, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 8, "H" => 5, "N" => 1, "O" => 1, "H" => 1]
chemicaldefaultelements(::Type{<:Melatonin}) = ["H" => 1, "N" => 1, "C" => 1, "O" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 2, "C" => 1, "H" => 2, "C" => 8, "H" => 5, "N" => 1, "O" => 1, "C" => 1, "H" => 3]
