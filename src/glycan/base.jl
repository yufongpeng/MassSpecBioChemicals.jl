length(::Monosaccharide) = 1
length(sugar::AbstractGlycan) = sum(length, getchaincomponent(sugar))
length(sugar::GlyComp) = sum(last, getchaincomponent(sugar))
