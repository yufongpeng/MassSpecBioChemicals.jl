module Proteins
using Reexport
@reexport using ..MassSpecBioChemicals
using ..MassSpecBioChemicals: lk
using MassSpecChemicals: AbstractChemical
import ..MassSpecBioChemicals: originalmolecule, leavinggroup, conjugation, ischainedchemical, dehydroxyposition, dehydrogenposition
import MassSpecChemicals: chemicalname, chemicalformula, chemicalabbr, repr_smiles
export AbstractPeptide,
        αAminoAcid, 
       Alanine,  
       Alanyl,  
       Arginine,  
       Arginyl,  
       Aspargine,  
       Asparginyl,  
       Aspartate,  
       Aspartyl,  
       Cysteine,  
       Cysteinyl,  
       Glutamine,  
       Glutaminyl,  
       Glutamate,  
       Glutamyl,  
       Glycine,  
       Glycyl,  
       Histidine,  
       Histidinyl,  
       Isoleucine,  
       Isoleucyl,  
       Leucine,  
       Leucyl,  
       Lysine,  
       Lysyl,  
       Methionine,  
       Methionyl,  
       Phenylalanine,  
       Phenylalanyl,  
       Proline,  
       Prolyl,  
       Serine,  
       Seryl,  
       Threonine,  
       Threonyl,  
       Tryptophan,  
       Tryptophanyl,  
       Tyrosine,  
       Tyrosyl,  
       Valine,  
       Valyl,  
       Selenocysteine,  
       Selenocysteinyl,  
       Pyrrolysine,  
       Pyrrolysyl,  
       Ornithine,  
       Ornithyl,  
       Oxytriptan,  
       Levodopa,
       Peptide,

       letter3_abbr,
       letter1_abbr

abstract type AbstractPeptide <: AbstractChemical end
abstract type αAminoAcid <: AbstractPeptide end

struct Alanine <: αAminoAcid end # A
struct Alanyl <: FunctionalGroup{Alanine, Dehydroxy} end # Ala
struct Arginine <: αAminoAcid end # R
struct Arginyl <: FunctionalGroup{Arginine, Dehydroxy} end # Arg
struct Aspargine <: αAminoAcid end # N
struct Asparginyl <: FunctionalGroup{Aspargine, Dehydroxy} end # Asn
struct Aspartate <: αAminoAcid end # D
struct Aspartyl <: FunctionalGroup{Aspartate, Dehydroxy} end # Asp
struct Cysteine <: αAminoAcid end # C
struct Cysteinyl <: FunctionalGroup{Cysteine, Dehydroxy} end # Cys
struct Glutamine <: αAminoAcid end # Q
struct Glutaminyl <: FunctionalGroup{Glutamine, Dehydroxy} end # Gln
struct Glutamate <: αAminoAcid end # E
struct Glutamyl <: FunctionalGroup{Glutamate, Dehydroxy} end # Glu
struct Glycine <: αAminoAcid end # G
struct Glycyl <: FunctionalGroup{Glycine, Dehydroxy} end # Gly
struct Histidine <: αAminoAcid end # H
struct Histidinyl <: FunctionalGroup{Histidine, Dehydroxy} end # His
struct Isoleucine <: αAminoAcid end # I
struct Isoleucyl <: FunctionalGroup{Isoleucine, Dehydroxy} end # Ile
struct Leucine <: αAminoAcid end # L
struct Leucyl <: FunctionalGroup{Leucine, Dehydroxy} end # Leu
struct Lysine <: αAminoAcid end # K
struct Lysyl <: FunctionalGroup{Lysine, Dehydroxy} end # Lys
struct Methionine <: αAminoAcid end # M
struct Methionyl <: FunctionalGroup{Methionine, Dehydroxy} end # Met
struct Phenylalanine <: αAminoAcid end # F
struct Phenylalanyl <: FunctionalGroup{Phenylalanine, Dehydroxy} end # Phe
struct Proline <: αAminoAcid end # P
struct Prolyl <: FunctionalGroup{Proline, Dehydroxy} end # Pro
struct Serine <: αAminoAcid end # S
struct Seryl <: FunctionalGroup{Serine, Dehydroxy} end # Ser
struct Threonine <: αAminoAcid end # T
struct Threonyl <: FunctionalGroup{Threonine, Dehydroxy} end # Thr
struct Tryptophan <: αAminoAcid end # W
struct Tryptophanyl <: FunctionalGroup{Tryptophan, Dehydroxy} end # Trp
struct Tyrosine <: αAminoAcid end # Y
struct Tyrosyl <: FunctionalGroup{Tyrosine, Dehydroxy} end # Tyr
struct Valine <: αAminoAcid end # V
struct Valyl <: FunctionalGroup{Valine, Dehydroxy} end # Val
struct Selenocysteine <: αAminoAcid end # U
struct Selenocysteinyl <: FunctionalGroup{Selenocysteine, Dehydroxy} end # Sec
struct Pyrrolysine <: αAminoAcid end # O
struct Pyrrolysyl <: FunctionalGroup{Pyrrolysine, Dehydroxy} end # Pyl
struct Ornithine <: αAminoAcid end # Orn
struct Ornithyl <: FunctionalGroup{Ornithine, Dehydroxy} end # Orn
struct Oxytriptan <: αAminoAcid end # 5HTP
struct Levodopa <: αAminoAcid end # LDOPA

struct Peptide{T} <: AbstractPeptide
    chain::T
end
ischainedchemical(::Peptide) = true
dehydrogenposition(::AbstractPeptide) = nothing
dehydroxyposition(::AbstractPeptide) = nothing
# dehydroxyposition(::Glutamate) = 0x01

conjugation(aa::αAminoAcid) = Substituent(Ndehydrogen, aa, lk(dehydrogenposition(aa)))
include("io.jl")

end