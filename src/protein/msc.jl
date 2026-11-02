struct AminoAcidParser{Bool, Bool} <: AbstractChemicalParser end
AminoAcidParser(; threeletter = true, functionalgroup = false) = AminoAcidParser{threeletter,functionalgroup}() 

parse_chemical(::AminoAcidParser{true, false}, x; kwargs...) = parse_aa3(x)
parse_chemical(::AminoAcidParser{false, false}, x; kwargs...) = parse_aa(x)
parse_chemical(::AminoAcidParser{true, true}, x; kwargs...) = parse_aa3_fg(x)
parse_chemical(::AminoAcidParser{false, true}, x; kwargs...) = parse_aa_fg(x)

chemicalname(m::AbstractPeptide; kwargs...) = m.name
chemicalelements(m::AbstractPeptide; kwargs...) = m.elements
chemicalformula(m::AbstractPeptide; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)
chemicalabbr(m::AbstractPeptide; kwargs...) = m.abbreviation
chemicalsmiles(m::AbstractPeptide; kwargs...) = ""

function chemicalabbr(aa:: αAminoAcid; nletter = 3, kwargs...)
    nletter == 3 ? aa.abbreviation : letter1_abbr(aa)
    # L <: LForm ? string("L-", s) : L <: DForm ? string("D-", s) : s
end
chemicalabbr(pep::Peptide; nletter = 3) =
    nletter == 3 ? pep.abbreviation : join([chemicalabbr(aa; nletter) for aa in pep.chain], "-")

chemicalelements(aa::FunctionalGroup{<: αAminoAcid, Dehydroxy}; kwargs...) = chemicalelements(parentchemical(aa))[begin:end - 2]
chemicalformula(aa::FunctionalGroup{<: αAminoAcid, Dehydroxy}; unique = false, kwargs...) = chemicalformula(chemicalelements(aa); unique, kwargs...)
