export isdissociated,
       Dihydrogen,
       Hydrogen,
       Ethane, 
       Ethyl, # Et
       Ethanol, 
       Ethoxy, # OEt
       Methane, 
       Methyl, # Me
       Methanol, 
       Methoxy, # OMe
       FormicAcid, 
       Formyl, # Fo
       Oformyl, # OFo
       Formamide, 
       Nformyl, # NFo
       AceticAcid, 
       Acetyl, # Ac
       Oacetyl, # OAc
       Acetamide, 
       Nacetyl, # NAc

       HydrogenBromide, 
       Bromo, # Br
       HydrogenChloride, 
       Chloro, # Cl
       HydrogenFluoride, 
       Fluoro, # F
       HydrogenIodide, 
       Iodo, # I
       HydrogenOxide, 
       OxygenAtom, 
       Hydroxy, # OH
       OLinkage, # O...
       CarboxylicAcidGroup, # COOH
       CarboxylicLinkage, # CO...
       Oxo, # oxo 2
       Ether, # oxy 2
       HydrogenPeroxide, 
       Hydroperoxyl, # OOH
       Epoxy, # Ep 2
       Peroxy, # OO 2

       Ammonia, 
       Amino, # NH2/N
       MethylAmine, 
       MethylAmino,
       DimethylAmine,
       DimethylAmino,
       TrimethylAmine,
       TrimethylAmino,
       NLinkage, # N...

       HydrogenSulfide, 
       Sulfanyl, # SH
       SLinkage, 
	   Methanethiol,
	   MethanethiolGroup,
       HydrogenCyanide, 
       Cyano, # CN
       PhosphoricAcid, # P
       Phosphoryl, # P/MP
       DiphosphoricAcid, 
       Diphosphoryl, # DP
       TriphosphoricAcid, 
       Triphosphoryl, # TP
       SulfuricAcid, # S
       Sulfo, # S
       # Sulfate, 
       NitricAcid, 
       Nitro, # NO2
       Phenolhydroxy,
       Amide,
       Ester,
       Thioester

using MassSpecChemicals: make_attributedchemical

struct Dihydrogen <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Dihydrogen}) = (; name = "H2", elements = ["H" => 2], abbreviation = "H2")
Dihydrogen(; kwargs...) = make_attributedchemical(Dihydrogen; kwargs...)
struct Hydrogen <: FunctionalGroup{Dihydrogen, Dehydrogen} end # H
struct Ethane <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Ethane}) = (; name = "Ethane", elements = ["C" => 1, "H" => 3, "C" => 1, "H" => 3], abbreviation = "Ethane")
Ethane(; kwargs...) = make_attributedchemical(Ethane; kwargs...)
struct Ethyl <: FunctionalGroup{Ethane, Dehydrogen} end # Et
struct Ethanol <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Ethanol}) = (; name = "Ethanol", elements = ["C" => 1, "H" => 3, "C" => 1, "H" => 2, "O" => 1, "H" => 1], abbreviation = "EtOH")
Ethanol(; kwargs...) = make_attributedchemical(Ethanol; kwargs...)
struct Ethoxy <: FunctionalGroup{Ethanol, Dehydrogen} end # OEt
struct Methane <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Methane}) = (; name = "Methane", elements = ["C" => 1, "H" => 4], abbreviation = "Methane")
Methane(; kwargs...) = make_attributedchemical(Methane; kwargs...)
struct Methyl <: FunctionalGroup{Methane, Dehydrogen} end # Me
struct Methanol <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Methanol}) = (; name = "Methanol", elements = ["C" => 1, "H" => 3, "O" => 1, "H" => 1], abbreviation = "MeOH")
Methanol(; kwargs...) = make_attributedchemical(Methanol; kwargs...)
struct Methoxy <: FunctionalGroup{Methanol, Dehydrogen} end # OMe
struct Formyl <: FunctionalGroup{FormicAcid, Dehydroxy} end # Fo
struct Oformyl <: FunctionalGroup{FormicAcid, Dehydrogen} end # OFo
struct Formamide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Formamide}) = (; name = "Formamide", elements = ["H" => 1, "C" => 1, "O" => 1, "N" => 1, "H" => 2], abbreviation = "Formamide")
Formamide(; kwargs...) = make_attributedchemical(Formamide; kwargs...)
struct Nformyl <: FunctionalGroup{Formamide, Dehydrogen} end # NFo
struct Acetyl <: FunctionalGroup{AceticAcid, Dehydroxy} end # Ac
struct Oacetyl <: FunctionalGroup{AceticAcid, Dehydrogen} end # OAc
struct Acetamide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Acetamide}) = (; name = "Acetamide", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "N" => 1, "H" => 2], abbreviation = "Acetamide")
Acetamide(; kwargs...) = make_attributedchemical(Acetamide; kwargs...)
struct Nacetyl <: FunctionalGroup{Acetamide, Dehydrogen} end # NAc
struct HydrogenBromide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenBromide}) = (; name = "Hydrogen bromide", elements = ["H" => 1, "Br" => 1], abbreviation = "Hydrogen bromide")
HydrogenBromide(; kwargs...) = make_attributedchemical(HydrogenBromide; kwargs...)
struct Bromo <: FunctionalGroup{HydrogenBromide, Dehydrogen} end # Br
struct HydrogenChloride <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenChloride}) = (; name = "Hydrogen chloride", elements = ["H" => 1, "Cl" => 1], abbreviation = "Hydrogen chloride")
HydrogenChloride(; kwargs...) = make_attributedchemical(HydrogenChloride; kwargs...)
struct Chloro <: FunctionalGroup{HydrogenChloride, Dehydrogen} end # Cl
struct HydrogenFluoride <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenFluoride}) = (; name = "Hydrogen fluoride", elements = ["H" => 1, "F" => 1], abbreviation = "Hydrogen fluoride")
HydrogenFluoride(; kwargs...) = make_attributedchemical(HydrogenFluoride; kwargs...)
struct Fluoro <: FunctionalGroup{HydrogenFluoride, Dehydrogen} end # F
struct HydrogenIodide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenIodide}) = (; name = "Hydrogen iodide", elements = ["H" => 1, "I" => 1], abbreviation = "Hydrogen iodide")
HydrogenIodide(; kwargs...) = make_attributedchemical(HydrogenIodide; kwargs...)
struct Iodo <: FunctionalGroup{HydrogenIodide, Dehydrogen} end # I
struct HydrogenOxide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenOxide}) = (; name = "Hydrogen oxide", elements = ["H" => 1, "O" => 1, "H" => 1], abbreviation = "Hydrogen oxide")
HydrogenOxide(; kwargs...) = make_attributedchemical(HydrogenOxide; kwargs...)
struct OxygenAtom <: UnknownGroup{HydrogenOxide, Dehydrogen} end 
struct Hydroxy <: FunctionalGroup{HydrogenOxide, Dehydrogen} end # OH
struct OLinkage <: FunctionalGroup{Hydroxy, Dehydrogen} end # O...
struct CarboxylicAcidGroup <: FunctionalGroup{FormicAcid, Demethine} end # COOH
struct CarboxylicLinkage <: FunctionalGroup{CarboxylicAcidGroup, Dehydroxy} end # CO...
struct Oxo <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxo 2
struct Ether <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxy 2
struct HydrogenPeroxide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenPeroxide}) = (; name = "Hydrogen peroxide", elements = ["H" => 1, "O" => 1, "O" => 1, "H" => 1], abbreviation = "Hydrogen peroxide")
HydrogenPeroxide(; kwargs...) = make_attributedchemical(HydrogenPeroxide; kwargs...)
struct Hydroperoxyl <: FunctionalGroup{HydrogenPeroxide, Dehydrogen} end # OOH
struct Epoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # Ep 2
struct Peroxy <: FunctionalGroup{HydrogenPeroxide, Didehydrogen} end # OO 2
struct Amino <: FunctionalGroup{Ammonia, Dehydrogen} end # NH2/N
struct NLinkage <: FunctionalGroup{Amino, Dehydrogen} end # N...
struct MethylAmine <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{MethylAmine}) = (; name = "Methylamine", elements = ["N" => 1, "H" => 2, "C" => 1, "H" => 3], abbreviation = "Methylamine")
MethylAmine(; kwargs...) = make_attributedchemical(MethylAmine; kwargs...)
struct MethylAmino <: FunctionalGroup{MethylAmine, Dehydrogen} end 
struct DimethylAmine <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{DimethylAmine}) = (; name = "Dimethylamine", elements = ["N" => 1, "H" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 3], abbreviation = "Dimethylamine")
DimethylAmine(; kwargs...) = make_attributedchemical(DimethylAmine; kwargs...)
struct DimethylAmino <: FunctionalGroup{DimethylAmine, Dehydrogen} end 
struct TrimethylAmine <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{TrimethylAmine}) = (; name = "Trimethylamine", elements = ["N" => 1, "C" => 1, "H" => 3, "C" => 1, "H" => 3, "C" => 1, "H" => 3], abbreviation = "Trimethylamine")
TrimethylAmine(; kwargs...) = make_attributedchemical(TrimethylAmine; kwargs...)
struct TrimethylAmino <: FunctionalGroup{TrimethylAmine, Dehydrogen} end 

struct HydrogenSulfide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenSulfide}) = (; name = "Hydrogen sulfide", elements = ["H" => 2, "S" => 1], abbreviation = "Hydrogen sulfide")
HydrogenSulfide(; kwargs...) = make_attributedchemical(HydrogenSulfide; kwargs...)
struct Sulfanyl <: FunctionalGroup{HydrogenSulfide, Dehydrogen} end # SH
struct Methanethiol <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{Methanethiol}) = (; name = "Methanethiol", elements = ["H" => 1, "S" => 1, "C" => 1, "H" => 3], abbreviation = "Methanethiol")
Methanethiol(; kwargs...) = make_attributedchemical(Methanethiol; kwargs...)
struct MethanethiolGroup <: FunctionalGroup{Methanethiol, Dehydrogen} end
struct HydrogenCyanide <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{HydrogenCyanide}) = (; name = "Hydrogen cyanide", elements = ["H" => 1, "C" => 1, "N" => 1], abbreviation = "Hydrogen cyanide")
HydrogenCyanide(; kwargs...) = make_attributedchemical(HydrogenCyanide; kwargs...)
struct Cyano <: FunctionalGroup{HydrogenCyanide, Dehydrogen} end # CN
struct PhosphoricAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{PhosphoricAcid}) = (; name = "Phosphoric acid", elements = ["H" => 3, "P" => 1, "O" => 4], abbreviation = "Phosphoric acid")
PhosphoricAcid(; kwargs...) = make_attributedchemical(PhosphoricAcid; kwargs...)
struct Phosphoryl <: FunctionalGroup{PhosphoricAcid, Dehydroxy} end # P/MP
struct DiphosphoricAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{DiphosphoricAcid}) = (; name = "Diphosphoric acid", elements = ["H" => 4, "P" => 2, "O" => 7], abbreviation = "Diphosphoric acid")
DiphosphoricAcid(; kwargs...) = make_attributedchemical(DiphosphoricAcid; kwargs...)
struct Diphosphoryl <: FunctionalGroup{DiphosphoricAcid, Dehydroxy} end # DP
struct TriphosphoricAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{TriphosphoricAcid}) = (; name = "Triphosphoric acid", elements = ["H" => 5, "P" => 3, "O" => 10], abbreviation = "Triphosphoric acid")
TriphosphoricAcid(; kwargs...) = make_attributedchemical(TriphosphoricAcid; kwargs...)
struct Triphosphoryl <: FunctionalGroup{TriphosphoricAcid, Dehydroxy} end # TP
struct SulfuricAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{SulfuricAcid}) = (; name = "Sulfuric acid", elements = ["H" => 2, "S" => 1, "O" => 4], abbreviation = "Sulfuric acid")
SulfuricAcid(; kwargs...) = make_attributedchemical(SulfuricAcid; kwargs...)
struct Sulfo <: FunctionalGroup{SulfuricAcid, Dehydroxy} end # S
struct SLinkage <: FunctionalGroup{Sulfo, Dehydrogen} end # S...

# struct Sulfate <: FunctionalGroup{SulfuricAcid, Dehydrogen} end 
struct NitricAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end
defaultchemicalkwargs(::Type{NitricAcid}) = (; name = "Nitric acid", elements = ["H" => 1, "N" => 1, "O" => 3], abbreviation = "Nitric acid")
NitricAcid(; kwargs...) = make_attributedchemical(NitricAcid; kwargs...)
struct Nitro <: FunctionalGroup{NitricAcid, Dehydroxy} end # NO2

struct Phenolhydroxy <: AbstractFunctionalGroup end
struct Amide <: AbstractFunctionalGroup end
struct Ester <: AbstractFunctionalGroup end
struct Thioester <: AbstractFunctionalGroup end

include("msc.jl")
include("interface.jl")