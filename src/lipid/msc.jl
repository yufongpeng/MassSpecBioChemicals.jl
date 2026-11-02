struct LipidParser <: AbstractChemicalParser end

parse_chemical(::LipidParser, x; kwargs...) = parse_lipid(x; kwargs...)

chemicalname(m::Lipid; kwargs...) = m.name
chemicalelements(m::Lipid; kwargs...) = m.elements
chemicalformula(m::Lipid; unique = false, kwargs...) =
    chemicalformula(chemicalelements(m); unique, kwargs...)
chemicalabbr(m::Lipid; kwargs...) = chemicalname(m; kwargs...)
chemicalsmiles(m::Lipid; kwargs...) = ""

chemicalelements(chain::CarbonChain; kwargs...) = chain.elements
chemicalelements(chain::Union{CycloChain, IsoprenoidChain}; kwargs...) = chain.elements
