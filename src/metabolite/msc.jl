chemicalname(m::Metabolite; kwargs...) = m.name
chemicalelements(m::Metabolite; kwargs...) = m.elements
chemicalabbr(m::Metabolite; kwargs...) = m.abbreviation
chemicalsmiles(m::Metabolite; kwargs...) = ""
chemicalformula(m::Metabolite; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)

chemicalsmiles(::Ethanolamine; kwargs...) = "NCCO"
chemicalsmiles(::Nmethylethanolamine; kwargs...) = "CNCCO"
chemicalsmiles(::NNdimethylethanolamine; kwargs...) = "CN(C)CCO"
chemicalsmiles(::GABA; kwargs...) = "NCCCC(=O)[O-]"
chemicalsmiles(::Dopamine; kwargs...) = "NCCc1cc(O)c(O)cc1"
chemicalsmiles(::Taurine; kwargs...) = "NCCS(=O)(=O)O"
chemicalsmiles(::Carnitine; kwargs...) = "C[N+](C)(C)CC(O)CC(=O)[O-]"
chemicalsmiles(::Choline; kwargs...) = "C[N+](C)(C)CCO"
chemicalsmiles(::CoA; kwargs...) = "SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O"
chemicalsmiles(::GlycolicAcid; kwargs...) = "OCC(=O)[O-]"
chemicalsmiles(::LacticAcid; kwargs...) = "OC(C)C(=O)[O-]"
chemicalsmiles(::PyruvicAcid; kwargs...) = "CC(=O)C(=O)[O-]"
function chemicalsmiles(::Glycerol; chain = (:O1, :C3), kwargs...)
    if chain == (:O1, :C3)
        "OCC(O)C(O)"
    elseif chain == (:C3, :O1)
        "C(O)C(O)CO"
    elseif chain == (:C3, :C1)
        "C(O)C(O)C(O)"
    elseif chain == (:C1, :C3)
        "C(O)C(O)C(O)"
    elseif chain == (:O3, :C1)
        "OCC(O)C(O)"
    elseif chain == (:C1, :O3)
        "C(O)C(O)CO"
    elseif chain == (:O3, :O1)
        "OCC(O)CO"
    elseif chain == (:O1, :O3)
        "OCC(O)CO"
    end
end
chemicalsmiles(::Serotonin; kwargs...) = "NCCC1=CNC2=C1C=C(C=C2)O"
chemicalsmiles(::Melatonin; kwargs...) = "N(C(=O)C)CCC1=CNC2=C1C=C(C=C2)OC"

chemicalabbr(::Tauryl; kwargs...) = "Tau"
chemicalabbr(::Glycolyl; kwargs...) = "Gc"
chemicalabbr(::Lactyl; kwargs...) = "Lt"
chemicalabbr(::Pyruvyl; kwargs...) = "Py"
chemicalabbr(::Glyceryl; kwargs...) = "Gr"
function chemicalname(::FunctionalGroup{T}; kwargs...) where {T<:Metabolite} 
    s = nameof(Base.typename(T).wrapper)
    T isa UnionAll ? string(repr_chirality(first(T.parameters)), s) : s 
end

function chemicalelements(x::T; kwargs...) where {T <: FunctionalGroup{<: Metabolite}}
    es = copy(chemicalelements(parentchemical(x)))
	ls = leavinggroupelements(leavinggroup(x))
	for (e, n) in ls 
		if n > 0
			i = findfirst(x -> ==(first(x), e), es)
			es[i] = e => (last(es[i]) + n)
		else 
			while n < 0 
				i = findfirst(x -> ==(first(x), e), es)
				es[i] = e => (last(es[i]) - 1)
				n += 1
				filter!(x -> last(x) > 0, es)
			end
		end
	end
	es
end 
chemicalformula(x::T; unique = false, kwargs...) where {T <: FunctionalGroup{<: Metabolite}} = chemicalformula(chemicalelements(x); unique, kwargs...)

chemicalsmiles(::Tauryl) = "(NCCS(=O)(=O)O)"
chemicalsmiles(::Glycolyl) = "(C(=O)CO)"
chemicalsmiles(::Lactyl) = "(C(=O)C(O)C)"
chemicalsmiles(::Pyruvyl) = "(C(=O)C(=O)C)"
chemicalsmiles(::Glyceryl) = "(OCC(O)C(O))"
