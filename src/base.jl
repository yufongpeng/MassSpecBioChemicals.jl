length(cc::ChainedChemical) = length(getchaincomponent(cc))
length(cc::DehydratedChemical) = length(getchaincomponent(cc))