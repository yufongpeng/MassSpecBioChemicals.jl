module MassSpecBioChemicals

using Reexport
@reexport using MassSpecChemicals
using MassSpecChemicals: AbstractChemical
using IterTools
import Base: show
import MassSpecChemicals: getchemicalattr
export FunctionalGroup, UnknownGroup, UnknownChemical, LeavingGroup, Dehydrogen, Odehydrogen, Ndehydrogen, Didehydrogen, Odidehydrogen, Ndidehydrogen, Deamine, Dehydroxy, Demethyl,
        Substituent, SubstitutedChemical, XLinkedFunctionalGroup, ChainedChemical, DehydratedChemical, IsotopiclabeledChemical,
        AbstractLinkageposition, Linkageposition, 
        parentchemical, leavinggroup, conjugation, includeSIL, getchainlinkage, getchaincomponent, ischainedchemical

include("functionalgroup.jl")

"""
    SubstitutedChemical{M <: AbstractChemical, S, T} <: AbstractChemical

Chemical with substituents. 

* `chemical`: parent chemical (`M`)
* `substituted`: leaving group and position information (`S`).
    * `LeavingGroup`: single substituent and unknown substituted position.
    * `Vector{Pair{LeavingGroup, UInt8}}`: leaving group and number of substitution pairs. 
    * `Vector{Pair{AbstractLinkageposition, LeavingGroup}}`: linkage position and leaving group pairs. 
* `substituent`: information of substutent(s) (`S`). 
    * `Vector{Pair{FunctionalGroup, UInt8}}`: functional group and number pairs. 
    * `Vector{Pair{AbstractLinkageposition, FunctionalGroup}}`: linkage position and functional group pairs. 
"""
struct SubstitutedChemical{M <: AbstractChemical, S, T} <: AbstractChemical
    chemical::M
    substituted::S
    substituent::T
end 

"""
    ChainedChemical{M, L} <: AbstractChemical

Linearly chained chemicals. 

* `chain`: a tuple of chemicals (`M`). The element is considered as a side branch if it is also a lineraly chained chemical, and the last chemical of side branch is conncected to the main branch via the latter element, e.g., `(A, B, (C, D), E)` has a main branch `(A, B, C)` and a side branch `(C, D)` connected to `E`.  
* `linkage`: linkage information.
    * `Nothing`: unknown position.
    * `Vector{Pair{LeavingGroup, LeavingGroup}}`: leaving group pairs.
    * `Vector{Vector{Pair{LeavingGroup, LeavingGroup}}, Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}}`: leaving group pairs and linkage position pairs. 
"""
struct ChainedChemical{M, L} <: AbstractChemical
    chain::M
    linkage::L
end

"""
    DehydratedChemical{M, L} <: AbstractChemical

Linearly chained chemicals from dehydration reactions. 

* `chain`: a tuple of chemicals (`M`). The element is considered as a side branch if it is also a lineraly chained chemical, and the last chemical of side branch is conncected to the main branch via the latter element. The order of chemicals follows the rule such that the former one loses a hydroxy group and the latter one loses a hydrogen atom. E.g., `(A, B, (C, D), E)` has a main branch `(A, B, C)` and a side branch `(C, D)` which `B` and `D` both lose hydroxy groups to connect to `E` losing two hydrogen atoms. 
* `linkage`: linkage information. 
    * `Nothing`: unknown position.
    * `Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}`: linkage position pairs.
"""
struct DehydratedChemical{M, L} <: AbstractChemical
    chain::M
    linkage::L
end

"""
    IsotopiclabeledChemical{M, I} <: AbstractChemical

Chemical with known isotopic label position. For chemicals with unknown or undefined isotopic label position, please refer to `MassSpecChemicals.Isotopomers`.

* `chemical`: parent chemical (`M`).
* `isotopiclabel`: isotopic label elements and position (`I`)
"""
struct IsotopiclabeledChemical{M, I} <: AbstractChemical
    chemical::M
    isotopiclabel::I
end

"""
    UnknownChemical{F <: UnknownGroup} <: AbstractChemical

Parent chemical of undetermined functional group. 
"""
struct UnknownChemical{F <: UnknownGroup} <: AbstractChemical end

include("attr.jl")
include("chemicals.jl")
include("io.jl")

include("BasicCompounds.jl")
include("Metabolites.jl")
include(joinpath("protein", "Proteins.jl"))
include(joinpath("glycan", "Glycans.jl"))
include(joinpath("nucleotide", "Nucleotides.jl"))
include(joinpath("lipid", "Lipids.jl"))

end