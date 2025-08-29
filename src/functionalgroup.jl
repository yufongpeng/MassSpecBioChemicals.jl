"""
    LeavingGroup

Abstract type for all leaving groups.
"""
abstract type LeavingGroup end
"""
    Dehydrogen <: LeavingGroup

Leaving group of losing one hydrogen atom. 
"""
struct Dehydrogen <: LeavingGroup end # Dehydrogen
"""
    Odehydrogen <: LeavingGroup

Leaving group of losing one hydrogen atom from hydroxy group. 
"""
struct Odehydrogen <: LeavingGroup end # Dehydroxy
"""
    Ndehydrogen <: LeavingGroup

Leaving group of losing one hydrogen atom from amine group.
"""
struct Ndehydrogen <: LeavingGroup end # Dehydroxy
"""
    Didehydrogen <: LeavingGroup

Leaving group of losing two hydrogen atoms. 
"""
struct Didehydrogen <: LeavingGroup end # Didehydrogen/Dehydrogen,Dehydrogen
# struct Odidehydrogen <: LeavingGroup end
# struct Ndidehydrogen <: LeavingGroup end
"""
    Dehydroxy <: LeavingGroup

Leaving group of losing one hydroxy group.
"""
struct Dehydroxy <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen
"""
    Deamine <: LeavingGroup

Leaving group of losing one amine group.
"""
struct Deamine <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen
"""
    Demethyl <: LeavingGroup

Leaving group of losing one methyl group.
"""
struct Demethyl <: LeavingGroup end # Dehydrogen/Odehydrogen/Ndehydrogen

"""
    AbstractLinkageposition

Abstract type for all linkage position types.
"""
abstract type AbstractLinkageposition end
"""
    Linkageposition <: AbstractLinkageposition

Normal linkage position type. 

* `position`: `Union{Nothing, UInt8}`
"""
struct Linkageposition <: AbstractLinkageposition
    position::Union{Nothing, UInt8}
end
"""
    lk(x)

User freiendly function to construct `Linkageposition`. 
"""
lk(x) = Linkageposition(UInt8(x))
lk(x::Nothing) = Linkageposition(x)

abstract type AbstractConfiguration end 
abstract type AbstractChirality <: AbstractConfiguration end
abstract type AbsoluteConfiguration <: AbstractChirality end
abstract type RSSystem <: AbsoluteConfiguration end 
struct AChirality <: RSSystem end 
struct RSChirality <: RSSystem end 
struct RChirality <: RSSystem end 
struct SChirality <: RSSystem end 

abstract type OpticRotationSystem <: AbsoluteConfiguration end 
struct NoOpticRotation <: OpticRotationSystem end 
struct UnknownRotation <: OpticRotationSystem end 
struct ClockwiseRotation <: OpticRotationSystem end 
struct CounterclockwiseRotation <: OpticRotationSystem end 

abstract type GlyceraldehydeSystem <: AbstractConfiguration end 
struct NoDLForm <: GlyceraldehydeSystem end 
struct DLForm <: GlyceraldehydeSystem end 
struct LForm <: GlyceraldehydeSystem end 
struct DForm <: GlyceraldehydeSystem end 

abstract type GeometricConfiguration <: AbstractConfiguration end
struct NoEZConfiguration <: GeometricConfiguration end
struct EZConfiguration <: GeometricConfiguration end
struct EConfiguration <: GeometricConfiguration end
struct ZConfiguration <: GeometricConfiguration end

"""
    AbstractFunctionalGroup <: AbstractChemical

Abstract type for all functional groups.
"""
abstract type AbstractFunctionalGroup <: AbstractChemical end
"""
    FunctionalGroup{M, S} <: AbstractFunctionalGroup

Abstract type for all known functional groups.
"""
abstract type FunctionalGroup{M, S} <: AbstractFunctionalGroup end
"""
    UnknownGroup{M, S} <: AbstractFunctionalGroup

Abstract type for all unknown or undetermined functional groups.
"""
abstract type UnknownGroup{M, S} <: AbstractFunctionalGroup end

"""
    Substituent{M <: AbstractChemical, S <: LeavingGroup} <: FunctionalGroup{M, S}

Substituent originated from parent chemical `M` and leaving group `S` at certain position. 

* `chemical`: parent chemical (`M`)
* `position`: linkage position. 
"""
struct Substituent{M <: AbstractChemical, S <: LeavingGroup} <: FunctionalGroup{M, S}
    chemical::M
    position::AbstractLinkageposition
end
Substituent(::Type{S}, m::M, p = Linkageposition(0x00)) where {M, S} = Substituent{M, S}(m, p)
Substituent(::S, m::M, p = Linkageposition(0x00)) where {M, S} = Substituent{M, S}(m, p)

"""
    XLinkedFunctionalGroup{M, L} <: AbstractFunctionalGroup

Functional groups linked with other chemical using additional linkagers. 

* `xlinkage`: linker (`M`).
* `functionalgroup`: functional group (`L`).
"""
struct XLinkedFunctionalGroup{M, L} <: AbstractFunctionalGroup
    xlinkage::M
    functionalgroup::L
end
# isdissociated check adjacent C=O
