function chiralitize(::Type{T}, args, kwargs) where {T<:Metabolite}
    if T <: UnionAll
        S = Core.apply_type(T, get(kwargs, :chirality, defaultchirality(T)))
        kwargs = Base.structdiff(NamedTuple(kwargs), NamedTuple{(:chirality,)})
    end
    S, args, kwargs
end

defaultchirality(::Type{Carnitine}) = DLForm
defaultchirality(::Type{LacticAcid}) = DLForm
defaultchirality(::Type{Glycerol}) = RSChirality
defaultchirality(::Type{Serotonin}) = RSChirality
defaultchirality(::Type{Melatonin}) = RSChirality

chiralchemical(::Type{Carnitine}, ::Type{I}) where {I <: GlyceraldehydeSystem} = Carnitine{I}()
chiralchemical(::Type{Carnitine}, ::Type{RSChirality}) = Carnitine{DLForm}()
chiralchemical(::Type{Carnitine}, ::Type{SChirality}) = Carnitine{DForm}()
chiralchemical(::Type{Carnitine}, ::Type{RChirality}) = Carnitine{LForm}()
chiralchemical(::Type{LacticAcid}, ::Type{I}) where {I <: GlyceraldehydeSystem} = LacticAcid{I}()
chiralchemical(::Type{LacticAcid}, ::Type{RSChirality}) = LacticAcid{DLForm}()
chiralchemical(::Type{LacticAcid}, ::Type{RChirality}) = LacticAcid{DForm}()
chiralchemical(::Type{LacticAcid}, ::Type{SChirality}) = LacticAcid{LForm}()
chiralchemical(::Type{Lactyl}, ::Type{I}) where {I <: GlyceraldehydeSystem} = Lactyl{I}()
chiralchemical(::Type{Lactyl}, ::Type{RSChirality}) = Lactyl{DLForm}()
chiralchemical(::Type{Lactyl}, ::Type{RChirality}) = Lactyl{DForm}()
chiralchemical(::Type{Lactyl}, ::Type{SChirality}) = Lactyl{LForm}()
chiralchemical(::Type{Glycerol}, ::Type{I}) where {I <: RSChirality} = Glycerol{I}()
chiralchemical(::Type{Glyceryl}, ::Type{I}) where {I <: RSChirality} = Glyceryl{I}()
chiralchemical(::Type{Serotonin}, ::Type{I}) where {I <: RSChirality} = Serotonin{I}()
chiralchemical(::Type{Melatonin}, ::Type{I}) where {I <: RSChirality} = Melatonin{I}()

dehydrogengroup(::Taurine; position = nothing) = Tauryl()
dehydroxygroup(::GlycolicAcid; position = nothing) = Glycolyl()
dehydroxygroup(::LacticAcid{T}; position = nothing) where T = Lactyl{T}()
dehydroxygroup(::PyruvicAcid; position = nothing) = Pyruvyl()
dehydrogengroup(::Glycerol{T}; position = nothing) where T = Glyceryl{T}()
conjugation(::Taurine) = Tauryl()

dehydrogenposition(::Ethanolamine) = nothing
dehydroxyposition(::Ethanolamine) = nothing
dehydrogenfunctionalgroup(::Ethanolamine; position = nothing) = Amino()
dehydroxyfunctionalgroup(::Ethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::Nmethylethanolamine) = nothing
dehydroxyposition(::Nmethylethanolamine) = nothing
dehydrogenfunctionalgroup(::Nmethylethanolamine; position = nothing) = MethylAmino()
dehydroxyfunctionalgroup(::Nmethylethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::NNdimethylethanolamine) = missing
dehydroxyposition(::NNdimethylethanolamine) = nothing
dehydroxyfunctionalgroup(::NNdimethylethanolamine; position = nothing) = Hydroxy()
dehydrogenposition(::GABA) = nothing
dehydroxyposition(::GABA) = nothing
dehydrogenfunctionalgroup(::GABA; position = nothing) = Amino()
dehydroxyfunctionalgroup(::GABA; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::Dopamine) = nothing
dehydroxyposition(::Dopamine) = missing
dehydrogenfunctionalgroup(::Dopamine; position = nothing) = Amino()
dehydrogenposition(::Taurine) = nothing
dehydrogenfunctionalgroup(::Taurine; position = nothing) = Sulfo()
dehydroxyposition(::Taurine) = missing
dehydrogenposition(::Choline) = missing
dehydroxyposition(::Choline) = nothing
dehydroxyfunctionalgroup(::Choline; position = nothing) = Hydroxy()
dehydrogenposition(::CoA) = nothing
dehydroxyposition(::CoA) = missing
dehydrogenfunctionalgroup(::CoA; position = nothing) = Sulfo()
dehydrogenposition(::GlycolicAcid) = nothing
dehydroxyposition(::GlycolicAcid) = nothing
dehydrogenfunctionalgroup(::GlycolicAcid; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::GlycolicAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::LacticAcid) = nothing
dehydroxyposition(::LacticAcid) = nothing
dehydrogenfunctionalgroup(::LacticAcid; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::LacticAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::PyruvicAcid) = missing
dehydroxyposition(::PyruvicAcid) = nothing
dehydroxyfunctionalgroup(::PyruvicAcid; position = nothing) = CarboxylicAcidGroup()
dehydrogenposition(::Glycerol) = 0x01
dehydroxyposition(::Glycerol) = 0x03
dehydrogenfunctionalgroup(::Glycerol; position = nothing) = Hydroxy()
dehydroxyfunctionalgroup(::Glycerol; position = nothing) = Hydroxy()

deisomerize(::Carnitine) = Carnitine{DLForm}()
deisomerize(::LacticAcid) = LacticAcid{DLForm}()
deisomerize(::Glycerol) = Glycerol{RSChirality}()

decompose(::Ethanolamine) = Dict(Hydroxy() => 0x01, Amino() => 0x01)
decompose(::Nmethylethanolamine) = Dict(Hydroxy() => 0x01, MethylAmino() => 0x01) # N-alky
decompose(::NNdimethylethanolamine) = Dict(Hydroxy() => 0x01, DimethylAmino() => 0x01)
decompose(::GABA) = Dict(CarboxylicAcidGroup() => 0x01, Amino() => 0x01)
decompose(::Dopamine) = Dict(Phenolhydroxy() => 0x02, Amino() => 0x01)
decompose(::Taurine) = Dict(Sulfo() => 0x01, Amino() => 0x01)
decompose(::Carnitine) = Dict(TrimethylAmino() => 0x01, Hydroxy() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Choline) = Dict(TrimethylAmino() => 0x01, Hydroxy() => 0x01)
decompose(::CoA) = Dict(CoA() => 0x01)
decompose(::GlycolicAcid) = Dict(CarboxylicAcidGroup() => 0x01, Hydroxy() => 0x01)
decompose(::LacticAcid) = Dict(CarboxylicAcidGroup() => 0x01, Hydroxy() => 0x01)
decompose(::PyruvicAcid) = Dict(CarboxylicAcidGroup() => 0x01, Oxo() => 0x01)
decompose(::Glycerol) = Dict(Hydroxy() => 0x03)
decompose(::Serotonin) = Dict(Phenolhydroxy() => 0x01, Amino() => 0x01, MethylAmino() => 0x01)
decompose(::Melatonin) = Dict(Methoxy() => 0x01, Nacetyl() => 0x01, MethylAmino() => 0x01)
