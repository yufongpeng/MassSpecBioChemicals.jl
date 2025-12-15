function push_dissociate_fg!(d, fgs)
    for fg in tuplize(fgs)
        if ispropertyposition(fg) 
            for (p, f) in fg
                get!(d, f, 0x00)
                d[f] += 0x01
            end
        else
            for (f, p) in fg 
                get!(d, f, 0x00)
                d[f] += p
            end
        end
    end
    d
end
function dissociate_carbonchain_group(lipid::Union{Hydrocarbon, FattyAlcohol, FattyAmide, FattyAldehyde}; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(lipid.backbone, chain) => nothing
    if ispropertyposition(fg) 
        d = Dict{AbstractChemical, UInt8}()
        for (p, f) in fg 
            get!(d, f, 0x00)
            d[f] += 0x01
        end
        fg = collect(pairs(d))
    end
    MonoFattyAcyl(lipid.backbone, chain) => (nothing, fg)
end
function dissociate_carbonchain_group(lipid::FattyAcid; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(lipid.backbone, chain) => nothing
    if ispropertyposition(fg) 
        d = Dict{AbstractChemical, UInt8}()
        for (p, f) in fg 
            get!(d, f, 0x00)
            d[f] += 0x01
        end
        get!(d, CarboxylicAcidGroup(), 0x00)
        d[CarboxylicAcidGroup()] += 0x01
        fg = collect(pairs(d))
    else
        i = findfirst(x -> ==(first(x), CarboxylicAcidGroup()), fg)
        isnothing(i) ? pushfirst!(fg, CarboxylicAcidGroup() => 0x01) : (fg[i] = CarboxylicAcidGroup() => (last(fg[i] + 0x01)))
    end
    MonoFattyAcyl(lipid.backbone, chain) => (CarboxylicAcidGroup(), fg)
end
function dissociate_carbonchain_group(lipid::FattyAmine; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(lipid.backbone, chain) => nothing
    if ispropertyposition(fg) 
        d = Dict{AbstractChemical, UInt8}()
        for (p, f) in fg 
            get!(d, f, 0x00)
            d[f] += 0x01
        end
        get!(d, Amino(), 0x00)
        d[Amino()] += 0x01
        fg = collect(pairs(d))
    else
        i = findfirst(x -> ==(first(x), Amino()), fg)
        isnothing(i) ? pushfirst!(fg, Amino() => 0x01) : (fg[i] = Amino() => (last(fg[i] + 0x01)))
    end
    MonoFattyAcyl(lipid.backbone, chain) => (Amino(), fg)
end
function dissociate_carbonchain_group(lipid::FattyAcylEstolid; presserve_fg = true) 
    # parallel
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_xlinkage(lipid.backbone.chain; presserve_fg)
    dl = FattyAcylEster(MonoFattyAcyl(lipid.backbone.backbone, chain2), chain1, lipid.position)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, (fg1, fg2))
    get!(d, CarboxylicAcidGroup(), 0x00)
    d[CarboxylicAcidGroup()] += 0x01
    fg = collect(pairs(d))
    dl => (CarboxylicAcidGroup(), fg)
end
function dissociate_carbonchain_group(lipid::WaxEster; presserve_fg = true) 
    # linear 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_xlinkage(lipid.backbone.chain; presserve_fg)
    dl = FattyAcylEster(MonoFattyAcyl(lipid.backbone.backbone, chain2), chain1, lipid.position)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, (fg1, fg2))
    fg = collect(pairs(d))
    dl => (nothing, fg)
end

function dissociate_carbonchain_group(lipid::MonoFattyAcyl; presserve_fg = true) 
    chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    presserve_fg || return MonoFattyAcyl(HydrogenOxide(), chain) => nothing
    b = dehydrogengroup(lipid.backbone)
    if ispropertyposition(fg) 
        d = Dict{AbstractChemical, UInt8}()
        for (p, f) in fg 
            get!(d, f, 0x00)
            d[f] += 0x01
        end
        get!(d, b, 0x00)
        d[b] += 0x01
        fg = collect(pairs(d))
    else
        i = findfirst(x -> ==(first(x), b), fg)
        isnothing(i) ? pushfirst!(fg, b => 0x01) : (fg[i] = b => (last(fg[i] + 0x01)))
    end
    MonoFattyAcyl(HydrogenOxide(), chain) => (b, fg)
end
function dissociate_carbonchain_group(lipid::NacylAmine; presserve_fg = true)
    if lipid.backbone isa FattyAcyl
        # linear
        chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
        presserve_fg || return NacylAmine(lipid.backbone, chain) => nothing
        if ispropertyposition(fg) 
            d = Dict{AbstractChemical, UInt8}()
            for (p, f) in fg 
                get!(d, f, 0x00)
                d[f] += 0x01
            end
            fg = collect(pairs(d))
        end
        NacylAmine(lipid.backbone, chain) => (nothing, fg)
    else
        chain, fg = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
        presserve_fg || return MonoFattyAcyl(HydrogenOxide(), chain) => nothing
        b = dehydrogengroup(lipid.backbone)
        if ispropertyposition(fg) 
            d = Dict{AbstractChemical, UInt8}()
            for (p, f) in fg 
                get!(d, f, 0x00)
                d[f] += 0x01
            end
            get!(d, b, 0x00)
            d[b] += 0x01
            fg = collect(pairs(d))
        else
            i = findfirst(x -> ==(first(x), b), fg)
            isnothing(i) ? pushfirst!(fg, b => 0x01) : (fg[i] = b => (last(fg[i] + 0x01)))
        end
        MonoFattyAcyl(HydrogenOxide(), chain) => (b, fg)
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
    dl = cls(lipid.backbone, chains, lipid.sn, lipid.chirality) 
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    dl => (nothing, collect(pairs(d)))
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
    dl =  cls(last(getchaincomponent(lipid.backbone)), chains, UInt8(sum(sn .* ((nchainposition(cls) + 1) .^ ((n - 1):-1:0)))), chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    b = deletechemicalat(lipid.backbone, length(lipid.backbone))
    get!(d, b, 0x00)
    d[b] += 0x01
    dl => (b, collect(pairs(d)))
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
        dl = Diradylglycerol(last(getchaincomponent(lipid.backbone)), chains, UInt8(sum(sn .* (4 .^ (1:-1:0)))), first(tuplize(lipid.chirality)))
        presserve_fg || return dl => nothing
        d = Dict{AbstractChemical, UInt8}()
        push_dissociate_fg!(d, fgs)
        b = deletechemicalat(lipid.backbone, [1, length(lipid.backbone)])
        get!(d, b, 0x00)
        d[b] += 0x01
        dl => (b, collect(pairs(d)))
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
        dl = Triradylglycerol(last(getchaincomponent(lipid.backbone)), length(chains) > 1 ? chains : first(chains), UInt8(sum(sn .* (4 .^ (2:-1:0)))), first(tuplize(lipid.chirality))) 
        presserve_fg || return dl => nothing
        d = Dict{AbstractChemical, UInt8}()
        push_dissociate_fg!(d, fgs)
        b = deletechemicalat(lipid.backbone, [1, length(lipid.backbone)])
        get!(d, b, 0x00)
        d[b] += 0x01
        dl => (b, collect(pairs(d)))
    else
        _dissociate_carbonchain_group(lipid)
    end
end

function dissociate_carbonchain_group(lipid::Union{Bisradylglycerophosphate, Radylglycerophosphoglycerol}; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(tuplize(lipid.chain); presserve_fg)
    if lipid.sn isa UInt16 
        sn = lipid.sn 
    else
        sn = decode_sn(lipid)
        sn = UInt16(sum(sn .* (5 .^ ((length(sn) - 1):-1:0))))
    end
    todel = [length(lipid.backbone) - 2 * i + 1 for i in 1:(length(lipid.backbone) ÷ 2)]
    dl = Radyldiglycerol(deletechemicalat(lipid.backbone, todel), chains, sn, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    get!(d, lipid.backbone, 0x00)
    d[lipid.backbone] += 0x01
    dl => (lipid.backbone, collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::Bisradylglycerophosphoglycerol; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = Radyltriglycerol(deletechemicalat(lipid.backbone, [2, 4]), chains, lipid.sn, lipid.chirality) 
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    get!(d, lipid.backbone, 0x00)
    d[lipid.backbone] += 0x01
    dl => (lipid.backbone, collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::GlycerophosphoNacylethanolamine; presserve_fg = true)
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = Monoradylglycerol(last(getchaincomponent(lipid.backbone)), chains, 0x01, AChirality())
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    b = deletechemicalat(lipid.backbone, length(lipid.backbone)) 
    get!(d, b, 0x00)
    d[b] += 0x01
    dl => (b, collect(pairs(d)))
end

function dissociate_carbonchain_group(lipid::CeramideBone; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = SphingoBone(nothing, chains, 0x00, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    dl => (lipid.headgroup, collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::Sulfonolipid; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = SphingoBone(nothing, chains, 0x00, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    get!(d, Sulfo(), 0x00)
    d[Sulfo()] += 0x01
    dl => (Sulfo(), collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::SphingoBone; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = SphingoBone(nothing, chains, 0x00, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    get!(d, Amino(), 0x00)
    d[Amino()] += 0x01
    dl => ((lipid.headgroup, Amino()), collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::Lysosulfonolipid; presserve_fg = true) 
    chains, fgs = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    dl = SphingoBone(nothing, chains, 0x00, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, fgs)
    get!(d, Sulfo(), 0x00)
    d[Sulfo()] += 0x01
    get!(d, Amino(), 0x00)
    d[Amino()] += 0x01
    dl => ((Sulfo(), Amino()), collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::Acylceramide; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_xlinkage(lipid.headgroup.chain; presserve_fg)
    dl = SphingoBone(MonoFattyAcyl(lipid.headgroup.backbone, chain2), chain1, lipid.position, lipid.chirality)
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, (tuplize(fg1)..., fg2))
    dl => (nothing, collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::Acylsphingomyelin; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    chain2, fg2 = dissociate_carbonchain_xlinkage(first(lipid.headgroup).chain; presserve_fg)
    dl = SphingoBone(MonoFattyAcyl(first(lipid.headgroup).backbone, chain2), chain1, first(lipid.position), first(lipid.chirality))
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, (tuplize(fg1)..., fg2))
    dl => (last(lipid.headgroup), collect(pairs(d)))
end
function dissociate_carbonchain_group(lipid::MixSphingoBone; presserve_fg = true) 
    chain1, fg1 = dissociate_carbonchain_xlinkage(lipid.chain; presserve_fg)
    id = findall(x -> x isa MonoFattyAcyl, lipid.headgroup)
    chains, fgs = zip(dissociate_carbonchain_xlinkage.(lipid.headgroup[id]; presserve_fg)...)
    if length(id) > 1
        dl = MixSphingoBone(chains, chain1, lipid.position[id], lipid.chirality[id])
    else
        dl = SphingoBone(first(chains), chain1, lipid.position[first(id)], lipid.chirality[first(id)])
    end
    presserve_fg || return dl => nothing
    d = Dict{AbstractChemical, UInt8}()
    push_dissociate_fg!(d, (tuplize(fg1)..., tuplize(fgs)...))
    dl => (filter(x -> !isa(x, MonoFattyAcyl), lipid.headgroup), collect(pairs(d)))
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
            g = parentchemical(first(v))
            isnothing(g) && continue
            presserve_fg && push!(fg, g => last(v))
        elseif first(v) isa XLinkedFunctionalGroup
            c, f = dissociate_carbonchain_xlinkage(first(v).functionalgroup.chemical.chain)
            s[i] = XLinkedFunctionalGroup(first(v).xlinkage, Substituent(leavinggroup(first(v).functionalgroup), MonoFattyAcyl(first(v).functionalgroup.chemical.backbone, c), first(v).functionalgroup.position)) => last(v)
            presserve_fg && (fg = vcat(fg, f...))
        end
        presserve_fg && isdissociated(first(s[i])) && push!(fg, s[i])
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
            g = parentchemical(last(v))
            isnothing(g) && continue
            presserve_fg && push!(fg, first(v) => g)
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