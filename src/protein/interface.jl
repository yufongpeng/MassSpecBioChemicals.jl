getchaincomponent(pp::Peptide) = pp.chain
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

function getchemicalattr(pp::Peptide, ::Val{:elements}; kwargs...)
    es = [losswaterelements(m, i, i % length(dm)) for (i, m) in enumerate(getchaincomponent(pp))]
    vcat(es...)
end

getchemicalattr(::Alanine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Arginine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 6, "N" => 1, "H" => 1, "C" => 1, "N" => 1, "H" => 1, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Aspargine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Aspartate, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Cysteine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "S" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Glutamate, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 2, "H" => 4, "C" => 1, "O" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Glycine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Histidine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 3, "H" => 3, "N" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Isoleucine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 9, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Leucine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 9, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Lysine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 8, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Methionine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "S" => 1, "C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Phenylalanine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 5, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Proline, ::Val{:elements}; kwargs...) = ["H" => 1, "N" => 1, "C" => 4, "H" => 7, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Serine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Threonine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 2, "H" => 4, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Tryptophan, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 8, "H" => 6, "N" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Tyrosine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 4, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Valine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 7, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Selenocysteine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "Se" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Pyrrolysine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 4, "H" => 8, "N" => 1, "H" => 1, "C" => 1, "O" => 1, "C" => 4, "H" => 7, "N" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::Ornithine, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 3, "H" => 6, "N" => 1, "H" => 2, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::fiveHydroxytryptophan, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 8, "H" => 6, "N" => 1, "O" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(::DOPA, ::Val{:elements}; kwargs...) = ["H" => 2, "N" => 1, "C" => 1, "H" => 1, "C" => 1, "H" => 2, "C" => 6, "H" => 3, "O" => 1, "H" => 1, "O" => 1, "H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]
getchemicalattr(aa::αAminoAcid, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(aa); unique, kwargs...)

getchemicalattr(aa::FunctionalGroup{<: αAminoAcid, Dehydroxy}, ::Val{:elements}; kwargs...) = chemicalelements(parentchemical(aa))[begin:end - 2]
getchemicalattr(aa::FunctionalGroup{<: αAminoAcid, Dehydroxy}, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(aa); unique, kwargs...)

deisomerize(pep::Peptide) = Peptide(map(deisomerize, getchaincomponent(pep)))
deisomerize(::Glycine) = Glycine{NoDLForm}()
deisomerize(::Alanine) = Alanine{DLForm}()
deisomerize(::Arginine) = Arginine{DLForm}()
deisomerize(::Aspargine) = Aspargine{DLForm}()
deisomerize(::Aspartate) = Aspartate{DLForm}()
deisomerize(::Cysteine) = Cysteine{DLForm}()
deisomerize(::Glutamine) = Glutamine{DLForm}()
deisomerize(::Glutamate) = Glutamate{DLForm}()
deisomerize(::Histidine) = Histine{DLForm}()
deisomerize(::Isoleucine) = Isoleucine{DLForm}()
deisomerize(::Leucine) = Leucine{DLForm}()
deisomerize(::Lysine) = Lysine{DLForm}()
deisomerize(::Methionine) = Methionine{DLForm}()
deisomerize(::Phenylalanine) = Phenylalanine{DLForm}()
deisomerize(::Proline) = Proline{DLForm}()
deisomerize(::Serine) = Serine{DLForm}()
deisomerize(::Threonine) = Threonine{DLForm}()
deisomerize(::Tryptophan) = Tryptophan{DLForm}()
deisomerize(::Tyrosine) = Tyrosine{DLForm}()
deisomerize(::Valine) = Valine{DLForm}()
deisomerize(::Selenocysteine) = Selenocysteine{DLForm}()
deisomerize(::Pyrrolysine) = Pyrrolysine{DLForm}()
deisomerize(::Ornithine) = Ornithine{DLForm}() 
deisomerize(::fiveHydroxytryptophan) = fiveHydroxytryptophan{DLForm}()
deisomerize(::DOPA) = DOPA{DLForm}()

composition(pep::Peptide) = Dict(getchaincomponent(pep))
composition(aa::αAminoAcid) = Dict(aa => 0x01)

function decompose(pep::Peptide) 
    p = map(decompose, getchaincomponent(pep))
    d = Dict{AbstractChemical, UInt8}(CarboxylicAcidGroup() => 0x01)
    for c in p 
        for (k, v) in c 
            get!(d, k, 0x00)
            d[k] += v 
        end
        d[Amino()] -= 0x01
        d[CarboxylicAcidGroup()] -= 0x01
        d[Amide()] += 0x01 
    end
    d[Amide()] -= 0x01
    d
end

decompose(::Glycine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Alanine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Arginine) = Dict(MethylAmino() => 0x01, Amino() => 0x03, CarboxylicAcidGroup() => 0x01)
decompose(::Aspargine) = Dict(Amide() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Aspartate) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x02)
decompose(::Cysteine) = Dict(Sulfanyl() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Glutamine) = Dict(Amide() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Glutamate) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x02)
decompose(::Histidine) = Dict(DimethylAmino() => 0x01, MethylAmino() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Isoleucine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Leucine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Lysine) = Dict(Amino() => 0x02, CarboxylicAcidGroup() => 0x01)
decompose(::Methionine) = Dict(Methanethiol() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Phenylalanine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Proline) = Dict(MethylAmino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Serine) = Dict(Hydroxy() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Threonine) = Dict(Hydroxy() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Tryptophan) = Dict(MethylAmino() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Tyrosine) = Dict(Phenolhydroxy() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Valine) = Dict(Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Selenocysteine) = Dict(Sulfanyl() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Pyrrolysine) = Dict(DimethylAmino() => 0x01, Amide() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::Ornithine) = Dict(Amino() => 0x02, CarboxylicAcidGroup() => 0x01)
decompose(::fiveHydroxytryptophan) = Dict(Hydroxy() => 0x01, MethylAmino() => 0x01, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)
decompose(::DOPA) = Dict(Phenolhydroxy() => 0x02, Amino() => 0x01, CarboxylicAcidGroup() => 0x01)

dehydrogenfunctionalgroup(pep::Peptide; position = nothing) = dehydrogenfunctionalgroup(first(getchaincomponent(pep)))
dehydrogenfunctionalgroup(::αAminoAcid; position = nothing) = Amino()
dehydrogenfunctionalgroup(::Proline; position = nothing) = MethylAmino()
dehydroxyfunctionalgroup(::Peptide; position = nothing) = CarboxylicAcidGroup()
dehydroxyfunctionalgroup(::αAminoAcid; position = nothing) = CarboxylicAcidGroup()