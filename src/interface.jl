function getchemicalattr(dm::DehydratedChemical, ::Val{:elements}; kwargs...)
    es = [losswaterelements(m, i, i % length(dm)) for (i, m) in enumerate(getchaincomponent(dm))]
    vcat(es...)
end
getchemicalattr(x::XLinkedFunctionalGroup, ::Val{:elements}; kwargs...) = vcat(chemicalelements(x.xlinkage), chemicalelements(x.functionalgroup))
function getchemicalattr(x::Substituent, ::Val{:elements}; kwargs...) 
	es = chemicalelements(parentchemical(x))
	ls = leavinggroupelements(leavinggroup(x))
	fn_f = length(ls) > 1 ? findlast : findfirst
	for (e, n) in ls 
		if n > 0
			i = fn_f(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = fn_f(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	es
end
getchemicalattr(x::NullSubstituent, ::Val{:elements}; kwargs...) = Pair{String, Int}[]
getchemicalattr(x::NullSubstituent, ::Val{:formula}; kwargs...) = ""

# functionalgroup
getchemicalattr(x::Substituent{<: T, <: Union{Dehydrogen, Odehydrogen, Ndehydrogen}}, ::Val{:SMILES}; kwargs...) where T = getchemicalattr(x.chemical, :SMILES; kwargs...)
getchemicalattr(x::Substituent{<: T, <: Dehydroxy}, ::Val{:SMILES}; kwargs...) where T = replace(getchemicalattr(x.chemical, :SMILES; kwargs...), r"^O" => "")
getchemicalattr(x::XLinkedFunctionalGroup, ::Val{:SMILES}; kwargs...) = string("(", getchemicalattr(x.xlinkage, :SMILES; kwargs...), replace(getchemicalattr(x.functionalgroup, :SMILES; kwargs...), r"^\(" => "", r"\)$" => ""), ")")

length(cc::ChainedChemical) = length(getchaincomponent(cc))
length(cc::DehydratedChemical) = length(getchaincomponent(cc))