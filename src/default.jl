_chainname(chain) = join(chemicalname.(chain), "-")
_chainabbr(chain) = join(chemicalabbr.(chain), "-")
_chainelements(chain) = reduce(vcat, chemicalelements.(chain); init = Pair{String, Int}[])
chemicaldefaultname(::Type{<:ChainedChemical}, chain) = _chainname(chain)
chemicaldefaultabbr(::Type{<:ChainedChemical}, chain) = _chainabbr(chain)
chemicaldefaultelements(::Type{<:ChainedChemical}, chain) = _chainelements(chain)
chemicaldefaultname(::Type{<:DehydratedChemical}, chain) = _chainname(chain)
chemicaldefaultabbr(::Type{<:DehydratedChemical}, chain) = _chainabbr(chain)
function _dehydratedelements(chain)
    n = length(chain)
    n == 0 && return Pair{String, Int}[]
    reduce(vcat, (losswaterelements(m, i, i % n) for (i, m) in enumerate(chain)); init = Pair{String, Int}[])
end
chemicaldefaultelements(::Type{<:DehydratedChemical}, chain) = _dehydratedelements(chain)