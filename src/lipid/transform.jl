function dissociate_carbonchain_group(lipid::Union{Hydrocarbon, FattyAlcohol, FattyAmide, FattyAldehyde}; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    MonoFattyAcyl(lipid.backbone, chain) => (fg, )
end
function dissociate_carbonchain_group(lipid::FattyAcylEstolid; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_group(lipid.backbone; presserve_fg)
    FattyAcylEster(chain2, chain1, lipid.position) => (presserve_fg ? (first(fg2), fg1) : (nothing, ))
end
function dissociate_carbonchain_group(lipid::WaxEster; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_group(lipid.backbone; presserve_fg)
    FattyAcylEster(chain2, chain1, lipid.position) => (presserve_fg ? (first(fg2), fg1) : (nothing, ))
end
function dissociate_carbonchain_group(lipid::FattyAcid; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(lipid.backbone, chain) => (nothing, )
    fg = ispropertyposition(fg) ? pushfirst!(fg, 0x01 => CarboxylicAcidGroup()) : sort_chainmodification!(pushfirst!(fg, CarboxylicAcidGroup() => 0x01))
    MonoFattyAcyl(lipid.backbone, chain) => (fg, )
end
function dissociate_carbonchain_group(lipid::FattyAmine; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(lipid.backbone, chain) => (nothing, )
    fg = ispropertyposition(fg) ? pushfirst!(fg, 0x01 => Amino()) : sort_chainmodification!(pushfirst!(fg, Amino() => 0x01))
    MonoFattyAcyl(lipid.backbone, chain) => (fg, )
end
function dissociate_carbonchain_group(lipid::MonoFattyAcyl; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fg = ispropertyposition(fg) ? pushfirst!(fg, 0x01 => dehydrogengroup(lipid.backbone)) : sort_chainmodification!(pushfirst!(fg, dehydrogengroup(lipid.backbone) => 0x01))
    MonoFattyAcyl(HydrogenOxide(), chain) => (fg, )
end
function dissociate_carbonchain_group(lipid::NacylAmine; presserve_fg = true)
    if lipid.backbone isa FattyAcyl
        chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
        MonoFattyAcyl(lipid.backbone, chain) => (fg, )
    else
        chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
        presserve_fg || return MonoFattyAcyl(HydrogenOxide(), chain) => (nothing, )
        fg = ispropertyposition(fg) ? pushfirst!(fg, 0x01 => dehydrogengroup(lipid.backbone)) : sort_chainmodification!(pushfirst!(fg, dehydrogengroup(lipid.backbone) => 0x01))
        MonoFattyAcyl(HydrogenOxide(), chain) => (fg, )
    end
end

function dissociate_carbonchain_group(lipid::Glycerolipid; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    n = ncarbonchain(lipid)
    cls = @match n begin
        1   => Monoradylglycerol
        2   => Diradylglycerol
        3   => Triradylglycerol
    end
    cls(lipid.backbone, chains, lipid.sn, lipid.chirality) => (presserve_fg ? tuplize(fgs) : (nothing, ))
end

dissociate_carbonchain_group(lipid::Union{<: Glycerophospholipid, <: Omodifiedradylglycerol}; presserve_fg = true) = _dissociate_carbonchain_group(lipid; presserve_fg)
function _dissociate_carbonchain_group(lipid::Union{<: Glycerophospholipid, <: Omodifiedradylglycerol}; presserve_fg = true)
    sn = decode_sn(lipid)
    n = ncarbonchain(lipid)
    cls = @match n begin
        1   => Monoradylglycerol
        2   => Diradylglycerol
        3   => Triradylglycerol
    end
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chirality = first(tuplize(lipid.chirality))
    chirality = chirality == RChirality() ? SChirality() : chirality == SChirality() ? RChirality() : chirality
    fgs = presserve_fg ? (tuplize(fgs)..., deletechemicalat(lipid.backbone, length(lipid.backbone))) : (nothing, )
    cls(last(getchaincomponent(lipid.backbone)), chains, UInt8(sum(sn .* ((nchainposition(cls) + 1) .^ ((n - 1):-1:0)))), chirality) => fgs
end

function dissociate_carbonchain_group(lipid::Union{LysophosphatidylNmodifiedethanolamine, LysophosphatidylNmodifiedserine}; presserve_fg = true)
    if first(getchaincomponent(lipid.backbone)) isa FattyAcyl
        sn = decode_sn(lipid)
        fa = first(getchaincomponent(lipid.backbone)).chain
        if ncarbon(fa) > 0
            chain = (lipid.chain, fa)
            push!(sn, 3)
        else
            chain = lipid.chain
        end
        chains, fgs = dissociate_carbonchain_xlinkage(chain; presserve_fg)
        fgs = presserve_fg ? (fgs..., deletechemicalat(lipid.backbone, [1, length(lipid.backbone)]) => FattyAcyl) : (nothing, )
        Diradylglycerol(last(getchaincomponent(lipid.backbone)), chains, UInt8(sum(sn .* (4 .^ (1:-1:0)))), first(tuplize(lipid.chirality))) => fgs
    else
        _dissociate_carbonchain_group(lipid; presserve_fg)
    end
end

function dissociate_carbonchain_group(lipid::Union{PhosphatidylNmodifiedethanolamine, PhosphatidylNmodifiedserine}; presserve_fg = true)
    if first(getchaincomponent(lipid.backbone)) isa FattyAcyl
        sn = decode_sn(lipid)
        fa = first(getchaincomponent(lipid.backbone)).chain
        if ncarbon(fa) > 0
            chain = (tuplize(lipid.chain)..., fa)
            push!(sn, 3)
        else
            chain = tuplize(lipid.chain)
        end
        chains, fgs = dissociate_carbonchain_xlinkage(chain; presserve_fg)
        fgs = presserve_fg ? (fgs..., deletechemicalat(lipid.backbone, [1, length(lipid.backbone)]) => FattyAcyl) : (nothing, )
        Triradylglycerol(last(getchaincomponent(lipid.backbone)), length(chains) > 1 ? chains : first(chains), UInt8(sum(sn .* (4 .^ (2:-1:0)))), first(tuplize(lipid.chirality))) => fgs
    else
        _dissociate_carbonchain_group(lipid)
    end
end

function dissociate_carbonchain_group(lipid::Bisradylglycerophosphate; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = presserve_fg ? (tuplize(fgs)..., lipid.backbone) : (nothing, )
    Radyldiglycerol(deletechemicalat(lipid.backbone, 2), chains, lipid.sn, lipid.chirality) => fgs
end
function dissociate_carbonchain_group(lipid::Bisradylglycerophosphoglycerol; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = presserve_fg ? (tuplize(fgs)..., lipid.backbone) : (nothing, )
    Radyltriglycerol(deletechemicalat(lipid.backbone, [2, 4]), chains, lipid.sn, lipid.chirality) => fgs
end
function dissociate_carbonchain_group(lipid::GlycerophosphoNacylethanolamine; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = presserve_fg ? (fgs, deletechemicalat(lipid.backbone, length(lipid.backbone)) => FattyAcyl) : (nothing, )
    Monoradylglycerol(last(getchaincomponent(lipid.backbone)), chains, 0x01, AChirality()) => fgs
end

function dissociate_carbonchain_group(lipid::CeramideBone; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    SphingoBone(nothing, chains, 0x00, lipid.chirality) => tuplize(fgs)
end
function dissociate_carbonchain_group(lipid::Sulfonolipid; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = tuplize(fgs)
    presserve_fg || return SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
    push!(first(fgs), 0x01 => Sulfo())
    SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
end
function dissociate_carbonchain_group(lipid::SphingoBone; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = tuplize(fgs)
    presserve_fg || return SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
    push!(first(fgs), 0x02 => Amino())
    SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
end
function dissociate_carbonchain_group(lipid::Lysosulfonolipid; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    fgs = tuplize(fgs)
    presserve_fg || return SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
    push!(first(fgs), 0x01 => Sulfo())
    push!(first(fgs), 0x02 => Amino())
    SphingoBone(nothing, chains, 0x00, lipid.chirality) => fgs
end
function dissociate_carbonchain_group(lipid::Acylceramide; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2 = dissociate_carbonchain(lipid.headgroup)
    SphingoBone(chain2, chain1, lipid.position, lipid.chirality) => tuplize(fg1)
end
function dissociate_carbonchain_group(lipid::Acylsphingomyelin; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2 = dissociate_carbonchain(first(lipid.headgroup))
    SphingoBone(chain2, chain1, first(lipid.position), first(lipid.chirality)) => tuplize(fg1)
end
function dissociate_carbonchain_group(lipid::MixSphingoBone; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    id = findall(x -> x isa MonoFattyAcyl, lipid.headgroup)
    chains = dissociate_carbonchain.(lipid.headgroup[id])
    if length(id) > 1
        MixSphingoBone(chains, chain1, lipid.position[id], lipid.chirality[id]) => tuplize(fg1)
    else
        SphingoBone(first(chains), chain1, lipid.position[first(id)], lipid.chirality[first(id)]) => tuplize(fg1)
    end
end

dissociate_carbonchain(lipid) = first(dissociate_carbonchain_group(lipid; presserve_fg = false))

function dissociate_carbonchain_xlinkage(chain::Tuple; presserve_fg = true) 
    x = ntuple(i -> dissociate_carbonchain_xlinkage(chain[i]; presserve_fg), length(chain))
    first.(x), last.(x)
end
function dissociate_carbonchain_xlinkage(chain::CarbonChain{T}; presserve_fg = true) where T
    x = dissociate_carbonchain_xlinkage(chain.substituent; presserve_fg) 
    make_carbonchain(T, chain.carbon, chain.doublebond, first(x), chain.chirality, chain.isotopiclabel), last(x)
end

# new fg for MonoFattyAcyl sub ?
function dissociate_carbonchain_xlinkage(sub::Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}}; presserve_fg = true)
    s = deepcopy(sub)
    fg = presserve_fg ? Pair{AbstractChemical, UInt8}[] : nothing
    for (i, v) in enumerate(sub)
        if first(v) isa XLinkedFunctionalGroup && !(first(v).functionalgroup isa Substituent{<: MonoFattyAcyl})
            f = transform_xlinkge(first(v))
            s[i] = f => last(v)
            presserve_fg && push!(fg, parentchemical(first(v)) => last(v))
        elseif first(v) isa XLinkedFunctionalGroup
            c, f = dissociate_carbonchain_xlinkage(first(v).functionalgroup.chemical.chain)
            s[i] = XLinkedFunctionalGroup(first(v).xlinkage, Substituent(leavinggroup(first(v).functionalgroup), MonoFattyAcyl(first(v).functionalgroup.chemical.backbone, c), first(v).functionalgroup.position)) => last(v)
            presserve_fg && (fg = vcat(fg, f...))
        end
        presserve_fg && isdissociated(last(s[i])) && push!(fg, s[i])
    end
    sort_chainmodification!(s), sort_chainmodification!(fg)
end

function dissociate_carbonchain_xlinkage(sub::Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}}; presserve_fg = true)
    s = deepcopy(sub)
    fg = presserve_fg ? Pair{UInt8, AbstractChemical}[] : nothing
    for (i, v) in enumerate(sub)
        if last(v) isa XLinkedFunctionalGroup && !(last(v).functionalgroup isa Substituent{<: MonoFattyAcyl})
            f = transform_xlinkge(last(v))
            s[i] = first(v) => f 
            presserve_fg && push!(fg, first(v) => parentchemical(last(v)))
        elseif last(v) isa XLinkedFunctionalGroup
            c, f = dissociate_carbonchain_xlinkage(last(v).functionalgroup.chemical.chain)
            s[i] = first(v) => XLinkedFunctionalGroup(last(v).xlinkage, Substituent(leavinggroup(last(v).functionalgroup), MonoFattyAcyl(last(v).functionalgroup.chemical.backbone, c), last(v).functionalgroup.position))
            f = map(f) do y
                (first(y) + first(v) + 0x01) => last(y)
            end
            presserve_fg && (fg = vcat(fg, f...))
        end
        presserve_fg && isdissociated(last(s[i])) && push!(fg, s[i])
    end
    sort_chainmodification!(s), sort_chainmodification!(fg)
end
dissociate_carbonchain_xlinkage(sub::Vector; presserve_fg = true) = sub, empty(sub)
dissociate_carbonchain_xlinkage(sub::Nothing; presserve_fg = true) = sub, []

transform_xlinkge(::XLinkedFunctionalGroup{OLinkage}) = Hydroxy()
transform_xlinkge(::XLinkedFunctionalGroup{NLinkage}) = Amino()
transform_xlinkge(::XLinkedFunctionalGroup{CarboxylicLinkage}) = CarboxylicAcidGroup()
# to INCHI, Search lipid map?
# fatty acyl to INCHI, Search lipid map?