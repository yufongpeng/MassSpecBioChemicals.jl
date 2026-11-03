"""
    LipidParser(; explicit_config = false)

Lipid parser storing the saccharide parsers required for parsing headgroups, tailgroups and glycan-derived lipid classes.
"""
struct LipidParser <: AbstractChemicalParser
    monosaccharide::MonosaccharideParser
    seriesglycan::SeriesGlycanParser
    glycomp::GlycompParser
end
LipidParser(; explicit_config = false) = LipidParser(
    MonosaccharideParser(; explicit_config), 
    SeriesGlycanParser(; explicit_config), 
    GlycompParser(; explicit_config)
)

parse_chemical(p::LipidParser, x; kwargs...) = parse_lipid(p, x; kwargs...)

chemicalname(m::Lipid; kwargs...) = m.name
chemicalelements(m::Lipid; kwargs...) = m.elements
chemicalformula(m::Lipid; unique = false, kwargs...) =
    chemicalformula(chemicalelements(m); unique, kwargs...)
chemicalabbr(m::Lipid; kwargs...) = chemicalname(m; kwargs...)
chemicalsmiles(m::Lipid; kwargs...) = ""

chemicalelements(chain::CarbonChain; kwargs...) = chain.elements
chemicalelements(chain::Union{CycloChain, IsoprenoidChain}; kwargs...) = chain.elements
