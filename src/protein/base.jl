length(::αAminoAcid) = 1
length(pep::AbstractPeptide) = length(getchaincomponent(pep))
