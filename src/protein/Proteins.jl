module Proteins
using Reexport
using ..MassSpecBioChemicals
using ..MassSpecBioChemicals: lk, GlyceraldehydeSystem, RSSystem
using MassSpecChemicals: AbstractChemical
import ..MassSpecBioChemicals: parentchemical, leavinggroup, conjugation, chainedchemical, getchaincomponent, getchainlinkage, getchainconfig, ischainedchemical, requirelinkage, requireconfig, dehydroxyposition, dehydrogenposition, dehydrogengroup, dehydroxygroup, chiralchemical
import MassSpecChemicals: getchemicalattr
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
       fiveHydroxytryptophan,  
       DOPA,
       Peptide,

       letter3_abbr,
       letter1_abbr,
       alpha_rs

abstract type AbstractPeptide <: AbstractChemical end
abstract type αAminoAcid{DL} <: AbstractPeptide end

struct Alanine{DL} <: αAminoAcid{DL} end # A
struct Alanyl{DL} <: FunctionalGroup{Alanine{DL}, Dehydroxy} end # Ala
struct Arginine{DL} <: αAminoAcid{DL} end # R
struct Arginyl{DL} <: FunctionalGroup{Arginine{DL}, Dehydroxy} end # Arg
struct Aspargine{DL} <: αAminoAcid{DL} end # N
struct Asparginyl{DL} <: FunctionalGroup{Aspargine{DL}, Dehydroxy} end # Asn
struct Aspartate{DL} <: αAminoAcid{DL} end # D
struct Aspartyl{DL} <: FunctionalGroup{Aspartate{DL}, Dehydroxy} end # Asp
struct Cysteine{DL} <: αAminoAcid{DL} end # C
struct Cysteinyl{DL} <: FunctionalGroup{Cysteine{DL}, Dehydroxy} end # Cys
struct Glutamine{DL} <: αAminoAcid{DL} end # Q
struct Glutaminyl{DL} <: FunctionalGroup{Glutamine{DL}, Dehydroxy} end # Gln
struct Glutamate{DL} <: αAminoAcid{DL} end # E
struct Glutamyl{DL} <: FunctionalGroup{Glutamate{DL}, Dehydroxy} end # Glu
struct Glycine{DL} <: αAminoAcid{DL} end # G
struct Glycyl{DL} <: FunctionalGroup{Glycine{DL}, Dehydroxy} end # Gly
struct Histidine{DL} <: αAminoAcid{DL} end # H
struct Histidinyl{DL} <: FunctionalGroup{Histidine{DL}, Dehydroxy} end # His
struct Isoleucine{DL} <: αAminoAcid{DL} end # I
struct Isoleucyl{DL} <: FunctionalGroup{Isoleucine{DL}, Dehydroxy} end # Ile
struct Leucine{DL} <: αAminoAcid{DL} end # L
struct Leucyl{DL} <: FunctionalGroup{Leucine{DL}, Dehydroxy} end # Leu
struct Lysine{DL} <: αAminoAcid{DL} end # K
struct Lysyl{DL} <: FunctionalGroup{Lysine{DL}, Dehydroxy} end # Lys
struct Methionine{DL} <: αAminoAcid{DL} end # M
struct Methionyl{DL} <: FunctionalGroup{Methionine{DL}, Dehydroxy} end # Met
struct Phenylalanine{DL} <: αAminoAcid{DL} end # F
struct Phenylalanyl{DL} <: FunctionalGroup{Phenylalanine{DL}, Dehydroxy} end # Phe
struct Proline{DL} <: αAminoAcid{DL} end # P
struct Prolyl{DL} <: FunctionalGroup{Proline{DL}, Dehydroxy} end # Pro
struct Serine{DL} <: αAminoAcid{DL} end # S
struct Seryl{DL} <: FunctionalGroup{Serine{DL}, Dehydroxy} end # Ser
struct Threonine{DL} <: αAminoAcid{DL} end # T
struct Threonyl{DL} <: FunctionalGroup{Threonine{DL}, Dehydroxy} end # Thr
struct Tryptophan{DL} <: αAminoAcid{DL} end # W
struct Tryptophanyl{DL} <: FunctionalGroup{Tryptophan{DL}, Dehydroxy} end # Trp
struct Tyrosine{DL} <: αAminoAcid{DL} end # Y
struct Tyrosyl{DL} <: FunctionalGroup{Tyrosine{DL}, Dehydroxy} end # Tyr
struct Valine{DL} <: αAminoAcid{DL} end # V
struct Valyl{DL} <: FunctionalGroup{Valine{DL}, Dehydroxy} end # Val
struct Selenocysteine{DL} <: αAminoAcid{DL} end # U
struct Selenocysteinyl{DL} <: FunctionalGroup{Selenocysteine{DL}, Dehydroxy} end # Sec
struct Pyrrolysine{DL} <: αAminoAcid{DL} end # O
struct Pyrrolysyl{DL} <: FunctionalGroup{Pyrrolysine{DL}, Dehydroxy} end # Pyl
struct Ornithine{DL} <: αAminoAcid{DL} end # Orn
struct Ornithyl{DL} <: FunctionalGroup{Ornithine{DL}, Dehydroxy} end # Orn
struct fiveHydroxytryptophan{DL} <: αAminoAcid{DL} end # 5HTP
struct DOPA{DL} <: αAminoAcid{DL} end # DOPA

struct Peptide{T} <: AbstractPeptide
    chain::T
end
getchaincomponent(m::Peptide) = m.chain
getchainlinkage(::Peptide) = missing
ischainedchemical(::Peptide) = true
ischainedchemical(::Type{<: Peptide}) = true
requirelinkage(::Type{<: Peptide}) = false
requireconfig(::Type{<: Peptide}) = false
chainedchemical(::Type{<: Peptide}, chemicals; kwargs...) = Peptide(chemicals)

getchainconfig(::Peptide) = missing
getchainconfig(::αAminoAcid) = missing
dehydrogenposition(::AbstractPeptide) = nothing
dehydroxyposition(::AbstractPeptide) = nothing
# dehydroxyposition(::Glutamate) = 0x01

dehydroxygroup(::Alanine{L}; position = nothing) where L = Alanyl{L}()
dehydroxygroup(::Arginine{L}; position = nothing) where L = Arginyl{L}()
dehydroxygroup(::Aspargine{L}; position = nothing) where L = Asparginyl{L}()
dehydroxygroup(::Aspartate{L}; position = nothing) where L = Aspartyl{L}()
dehydroxygroup(::Cysteine{L}; position = nothing) where L = Cysteinyl{L}()
dehydroxygroup(::Glutamine{L}; position = nothing) where L = Glutaminyl{L}()
dehydroxygroup(::Glutamate{L}; position = nothing) where L = Glutamyl{L}()
dehydroxygroup(::Glycine{L}; position = nothing) where L = Glycyl{L}()
dehydroxygroup(::Histidine{L}; position = nothing) where L = Histidinyl{L}()
dehydroxygroup(::Isoleucine{L}; position = nothing) where L = Isoleucyl{L}()
dehydroxygroup(::Leucine{L}; position = nothing) where L = Leucyl{L}()
dehydroxygroup(::Lysine{L}; position = nothing) where L = Lysyl{L}()
dehydroxygroup(::Methionine{L}; position = nothing) where L = Methionyl{L}()
dehydroxygroup(::Phenylalanine{L}; position = nothing) where L = Phenylalanyl{L}()
dehydroxygroup(::Proline{L}; position = nothing) where L = Prolyl{L}()
dehydroxygroup(::Serine{L}; position = nothing) where L = Seryl{L}()
dehydroxygroup(::Threonine{L}; position = nothing) where L = Threonyl{L}()
dehydroxygroup(::Tryptophan{L}; position = nothing) where L = Tryptophanyl{L}()
dehydroxygroup(::Tyrosine{L}; position = nothing) where L = Tyrosyl{L}()
dehydroxygroup(::Valine{L}; position = nothing) where L = Valyl{L}()
dehydroxygroup(::Selenocysteine{L}; position = nothing) where L = Selenocysteinyl{L}()
dehydroxygroup(::Pyrrolysine{L}; position = nothing) where L = Pyrrolysyl{L}()
dehydroxygroup(::Ornithine{L}; position = nothing) where L = Ornithyl{L}() 
dehydrogengroup(aa::αAminoAcid; position = nothing) = Substituent(Ndehydrogen, aa, isnothing(position) ? lk(dehydrogenposition(aa)) : position)
conjugation(aa::αAminoAcid) = dehydroxygroup(aa)
include("io.jl")

end