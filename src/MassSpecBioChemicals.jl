module MassSpecBioChemicals

using Reexport
using IterTools
@reexport using MassSpecChemicals
using MassSpecChemicals: 
    AbstractChemical, AttributedChemical, ChemicalAttribute, 
    pushchemicalattr_new!,
    buildattribute,
    safecopy, chemicalelements_equal, sort_unique_elements
    
import Base: show, length, ==, hash, copy

import MassSpecChemicals: 
    chemicalname, chemicalformula, chemicalelements, chemicalabbr, chemicalsmiles, 
    pushchemicalattr!, defaultchemicalkwargs

export FunctionalGroup, UnknownGroup, UnknownChemical, 
    LeavingGroup, Dehydrogen, Odehydrogen, Ndehydrogen, Didehydrogen, Odidehydrogen, Ndidehydrogen, Deamine, Dehydroxy, Demethine, Tridehydrogen,
    Substituent, NullSubstituent, XLinkedFunctionalGroup, 
    
    SubstitutedChemical, ChainedChemical, DehydratedChemical,
    
    AbstractLinkageposition, Linkageposition, 

    parentchemical, leavinggroup, conjugation,
    getchainlinkage, getchaincomponent, getchainconfig, composition, ischainedchemical, 

    AChirality, RSChirality, RChirality, SChirality, 
    NoOpticRotation, UnknownRotation, ClockwiseRotation, CounterclockwiseRotation, 
    NoDLForm, DLForm, LForm, DForm, 
    NoEZConfiguration, EZConfiguration, EConfiguration, ZConfiguration 
        
include(joinpath("functionalgroup", "functionalgroup.jl"))

"""
    SubstitutedChemical{M <: AbstractChemical, S, T, C} <: AbstractChemical

Chemical with substituents. 

* `chemical`: unsubstituted chemical (`M`).
* `substituted`: leaving group and position information of substituted chemical (`S`).
    * `LeavingGroup`: single substituent and unknown substituted position.
    * `Vector{Pair{LeavingGroup, UInt8}}`: leaving group and number of substitution pairs. 
    * `Vector{Pair{AbstractLinkageposition, LeavingGroup}}`: linkage position and leaving group pairs. 
* `substituent`: information of substutent(s) (`S`). 
    * `Vector{Pair{FunctionalGroup, UInt8}}`: functional group and number pairs. 
    * `Vector{Pair{AbstractLinkageposition, FunctionalGroup}}`: linkage position and functional group pairs. 
* `config`: configuration of chiral centers or cis-trans isomers related to the substituents (`C`). 
    * `Nothing`: unknown configuration information. 
    * `Vector{UInt8, AbstractConfiguration}`: position configuration pairs. 
"""
struct SubstitutedChemical{M <: AbstractChemical, S, T, C} <: AttributedChemical
    chemical::M
    substituted::S
    substituent::T
    config::C
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 

"""
    ChainedChemical{M, L, C} <: AbstractChemical

Linearly chained chemicals. 

* `chain`: a tuple of chemicals (`M`). The element is considered as a side branch if it is also a lineraly chained chemical, and the last chemical of side branch is connected to the main branch via the next non linearly chained element, e.g., `(A, B, (C, D), E)` has a main branch `(A, B, C)` and a side branch `(C, D)` connected to `E` via `D-E` bonding.  
* `linkage`: linkage information.
    * `Nothing`: unknown position.
    * `Vector{Pair{LeavingGroup, LeavingGroup}}`: leaving group pairs of each element of `chain` and the next connected element. 
    * `Vector{Vector{Pair{LeavingGroup, LeavingGroup}}, Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}}`: leaving group pairs and linkage position pairsof each element of `chain` and the next connected element. 
* `config`: configuration of chiral centers or cis-trans isomers related to the linkage (`C`). 
    * `Nothing`: unknown configuration information. 
    * `Vector{Union{Nothing, Vector{UInt8, AbstractConfiguration}}}`: position configuration pairs of each element of `chain`. 
"""
struct ChainedChemical{M, L, C} <: AttributedChemical
    chain::M
    linkage::L
    config::C
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end

"""
    DehydratedChemical{M, L, C} <: AbstractChemical

Linearly chained chemicals from dehydration reactions. 

* `chain`: a tuple of chemicals (`M`). The element is considered as a side branch if it is also a lineraly chained chemical, and the last chemical of side branch is connected to the main branch via the next non linearly chained element. The order of chemicals follows the rule such that the former one loses a hydroxy group and the latter one loses a hydrogen atom. E.g., `(A, B, (C, D), E)` has a main branch `(A, B, C)` and a side branch `(C, D)` which `B` and `D` both lose hydroxy groups to connect to `E` losing two hydrogen atoms. 
* `linkage`: linkage information. 
    * `Nothing`: unknown position.
    * `Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}`: linkage position pairs of each element of `chain` and the next connected element. 
* `config`: configuration of chiral centers or cis-trans isomers related to the linkage (`C`). 
    * `Nothing`: unknown configuration information. 
    * `Vector{Union{Nothing, Vector{UInt8, AbstractConfiguration}}}`: position configuration pairs of each element of `chain`. 
"""
struct DehydratedChemical{M, L, C} <: AttributedChemical
    chain::M
    linkage::L
    config::C
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end

"""
    UnknownChemical{F <: UnknownGroup} <: AbstractChemical

Parent chemical of undetermined functional group. 
"""
struct UnknownChemical{F <: UnknownGroup} <: AbstractChemical end

include("base.jl")
include("msc.jl")
include("chemicals.jl")
include("default.jl")
include("interface.jl")
# include("io.jl")

include(joinpath("basic", "basic.jl"))
include(joinpath("metabolite", "Metabolites.jl"))
include(joinpath("protein", "Proteins.jl"))
include(joinpath("glycan", "Glycans.jl"))
include(joinpath("nucleotide", "Nucleotides.jl"))
include(joinpath("lipid", "Lipids.jl"))

end