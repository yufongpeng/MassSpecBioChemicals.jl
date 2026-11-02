isdissociated(::Hydrogen) = false
isdissociated(::Ethyl) = false
isdissociated(::Ethoxy) = false
isdissociated(::Methyl) = false
isdissociated(::Methoxy) = false
isdissociated(::Formyl) = false
isdissociated(::Oformyl) = false
isdissociated(::Nformyl) = false
isdissociated(::Acetyl) = false
isdissociated(::Oacetyl) = false
isdissociated(::Nacetyl) = false
isdissociated(::Bromo) = false
isdissociated(::Chloro) = false
isdissociated(::Fluoro) = false
isdissociated(::Iodo) = false
isdissociated(::OxygenAtom) = false
isdissociated(::Hydroxy) = false
isdissociated(::OLinkage) = false
isdissociated(::CarboxylicAcidGroup) = true
isdissociated(::CarboxylicLinkage) = false
isdissociated(::Oxo) = false
isdissociated(::Ether) = false
isdissociated(::Hydroperoxyl) = false
isdissociated(::Epoxy) = false
isdissociated(::Peroxy) = false
isdissociated(::Amino) = true
isdissociated(::NLinkage) = true
isdissociated(::Sulfanyl) = true
isdissociated(::SLinkage) = false
isdissociated(::MethanethiolGroup) = false
isdissociated(::Cyano) = false
isdissociated(::Phosphoryl) = true
isdissociated(::Substituent{Phosphoryl, Dehydroxy}) = true
isdissociated(::Diphosphoryl) = true
isdissociated(::Triphosphoryl) = true
isdissociated(::Sulfo) = true
isdissociated(::Nitro) = false
isdissociated(::Phenolhydroxy) = true
isdissociated(::Amide) = false
isdissociated(::Ester) = false
isdissociated(::Thioester) = false

dehydrogenposition(::AttributedChemical) = nothing
dehydroxyposition(::AttributedChemical) = missing
dehydroxyposition(::Methanol) = nothing
dehydroxyposition(::Ethanol) = nothing
dehydroxyposition(::FormicAcid) = nothing
dehydroxyposition(::AceticAcid) = nothing
dehydroxyposition(::PhosphoricAcid) = nothing
dehydroxyposition(::DiphosphoricAcid) = nothing
dehydroxyposition(::TriphosphoricAcid) = nothing
dehydroxyposition(::SulfuricAcid) = nothing
dehydroxyposition(::NitricAcid) = nothing

dehydrogengroup(::Dihydrogen; position = nothing) = Hydrogen()
dehydrogengroup(::Ethane; position = nothing) = Ethyl()
dehydrogengroup(::Ethanol; position = nothing) = Ethoxy()
dehydrogengroup(::Methane; position = nothing) = Methyl()
dehydrogengroup(::Methanol; position = nothing) = Methoxy()
dehydroxygroup(::FormicAcid; position = nothing) = Formyl()
dehydrogengroup(::FormicAcid; position = nothing) = Oformyl()
dehydrogengroup(::Formamide; position = nothing) = Nformyl()
dehydroxygroup(::AceticAcid; position = nothing) = Acetyl()
dehydrogengroup(::AceticAcid; position = nothing) = Oacetyl()
dehydrogengroup(::Acetamide; position = nothing) = Nacetyl()
dehydrogengroup(::HydrogenBromide; position = nothing) = Bromo()
dehydrogengroup(::HydrogenChloride; position = nothing) = Chloro()
dehydrogengroup(::HydrogenFluoride; position = nothing) = Fluoro()
dehydrogengroup(::HydrogenIodide; position = nothing) = Iodo()
dehydrogengroup(::HydrogenOxide; position = nothing) = Hydroxy()
dehydrogengroup(::HydrogenPeroxide; position = nothing) = Hydroperoxyl()
dehydrogengroup(::Ammonia; position = nothing) = Amino()
dehydrogengroup(::MethylAmine; position = nothing) = MethylAmino()
dehydrogengroup(::DimethylAmine; position = nothing) = DimethylAmino()
dehydrogengroup(::TrimethylAmine; position = nothing) = TrimethylAmino()
dehydrogengroup(::HydrogenSulfide; position = nothing) = Sulfanyl()
dehydrogengroup(::Methanethiol; position = nothing) = MethanethiolGroup()
dehydrogengroup(::HydrogenCyanide; position = nothing) = Cyano()
dehydroxygroup(::PhosphoricAcid; position = nothing) = Phosphoryl()
dehydroxygroup(::DiphosphoricAcid; position = nothing) = Diphosphoryl()
dehydroxygroup(::TriphosphoricAcid; position = nothing) = Triphosphoryl()
dehydroxygroup(::SulfuricAcid; position = nothing) = Sulfo()
# dehydrogengroup(::SulfuricAcid; position = nothing) = Sulfate()
dehydroxygroup(::NitricAcid; position = nothing) = Nitro()

nlinkage(::CarboxylicAcidGroup) = 3
nlinkage(::OxygenAtom) = 0
nlinkage(::Ether) = 1
nbridge(::Ether) = 1
nlinkage(::Epoxy) = 1
nbridge(::Epoxy) = 1
nlinkage(::Peroxy) = 1
nbridge(::Peroxy) = 1
nlinkage(::XLinkedFunctionalGroup{CarboxylicLinkage}) = 3

dehydroxyfunctionalgroup(::FormicAcid; position = nothing) = CarboxylicAcidGroup()
dehydroxyfunctionalgroup(::AceticAcid; position = nothing) = CarboxylicAcidGroup()

function decompose(fg::XLinkedFunctionalGroup{CarboxylicLinkage}) 
    d = decompose(fg.functionalgroup)
    curr = dehydrogenfunctionalgroup(fg.functionalgroup)
    if curr == Hydroxy()
        d[Ester()] += 0x01
    elseif curr == Amino() || curr == MethylAmino() 
        d[Amide()] += 0x01
    elseif curr == Sulfo()
        d[Thioester()] += 0x01
    else
        d[Oxo()] += 0x01
    end
    d
end

linkage_condense(::Type{LeavingGroup}, prev::CarboxylicAcidGroup, curr::AbstractFunctionalGroup) = ([], [])
linkage_condense(::Type{Dehydroxy}) = ([Hydroxy()], [])
linkage_condense(::Type{Dehydroxy}, prev::AbstractFunctionalGroup, curr::AbstractFunctionalGroup) = ([prev], [])
linkage_condense(::Type{Dehydroxy}, prev::CarboxylicAcidGroup, curr::Union{Hydroxy, Phenolhydroxy}) = ([curr, prev], [Ester()])
linkage_condense(::Type{Dehydroxy}, prev::CarboxylicAcidGroup, curr::Union{Amino, MethylAmino}) = ([curr, prev], [Amide()])
linkage_condense(::Type{Dehydroxy}, prev::CarboxylicAcidGroup, curr::Sulfanyl) = ([curr, prev], [Thioester()])
linkage_condense(::Type{Dehydroxy}, prev::CarboxylicAcidGroup, curr::AbstractFunctionalGroup) = ([prev], [Oxo()])
linkage_condense(::Type{Dehydroxy}, prev::AbstractFunctionalGroup, curr::Union{Hydroxy, Phenolhydroxy}) = ([prev, curr], [OLinkage()])
linkage_condense(::Type{Dehydroxy}, prev::AbstractFunctionalGroup, curr::Union{Amino, MethylAmino}) = ([prev, curr], [NLinkage()])
linkage_condense(::Type{Dehydroxy}, prev::AbstractFunctionalGroup, curr::Sulfanyl) = ([prev, curr], [SLinkage()])

linkage_condense(::Type{Dehydroxy}, m::AbstractChemical; position = nothing) = ([Hydroxy()], [])

dehydrogenfunctionalgroup(a::AbstractChemical; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(a::AbstractChemical; position = nothing) = Hydroxy()
