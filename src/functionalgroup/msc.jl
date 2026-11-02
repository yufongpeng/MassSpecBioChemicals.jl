chemicalname(cc::AbstractFunctionalGroup; kwargs...) = string(typeof(cc))
chemicalelements(cc::AbstractFunctionalGroup; kwargs...) = Pair{String, Int}[]
chemicalformula(cc::AbstractFunctionalGroup; kwargs...) = chemicalformula(chemicalelements(cc); kwargs...)
chemicalabbr(cc::AbstractFunctionalGroup; kwargs...) = chemicalname(cc; kwargs...)
chemicalsmiles(cc::AbstractFunctionalGroup; kwargs...) = ""

chemicalelements(x::XLinkedFunctionalGroup; kwargs...) = vcat(chemicalelements(x.xlinkage; kwargs...), chemicalelements(x.functionalgroup; kwargs...))
function chemicalelements(x::Substituent; kwargs...) 
	es = copy(chemicalelements(parentchemical(x)))
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
chemicalelements(x::NullSubstituent; kwargs...) = Pair{String, Int}[]
chemicalformula(x::NullSubstituent; kwargs...) = ""
chemicalname(x::T; group = true, kwargs...) where {T <: FunctionalGroup{<: AttributedChemical}} = group ? string(T, " Group") : string(T)

function chemicalelements(x::T; kwargs...) where {T <: FunctionalGroup{<: AttributedChemical}} 
	es = copy(chemicalelements(parentchemical(x)))
	ls = leavinggroupelements(leavinggroup(x))
	fn_f = length(ls) > 1  ? findlast : findfirst
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
chemicalformula(x::T; unique = false, kwargs...) where {T <: FunctionalGroup{<: AttributedChemical}} = chemicalformula(chemicalelements(x); unique, kwargs...)

# functionalgroup
chemicalsmiles(x::Substituent{<: T, <: Union{Dehydrogen, Odehydrogen, Ndehydrogen}}; kwargs...) where T = chemicalsmiles(x.chemical; kwargs...)
chemicalsmiles(x::Substituent{<: T, <: Dehydroxy}; kwargs...) where T = replace(chemicalsmiles(x.chemical; kwargs...), r"^O" => "")
chemicalsmiles(x::XLinkedFunctionalGroup; kwargs...) = string("(", chemicalsmiles(x.xlinkage; kwargs...), replace(chemicalsmiles(x.functionalgroup; kwargs...), r"^\(" => "", r"\)$" => ""), ")")
