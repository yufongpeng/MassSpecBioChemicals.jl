check_structure!(::Type{<: FattyAcyl}, bone, chain::CarbonChain) = 
    (bone, check_configuration!(bone, chain))
check_structure!(::Type{<: FattyAcyl}, bone, chain) = throw(ArgumentError("Invalid number of chains"))
check_structure!(::Type{<: FattyAcylEster}, bone, chain::CarbonChain, sn) = 
    (bone, check_configuration!(bone, chain), check_ester_position(chain, sn))
check_structure!(::Type{<: FattyAcylEster}, bone, chain, sn) = throw(ArgumentError("Invalid number of chains"))

function check_ester_position(chain, position)
    position == 0 && return position
    !(chain isa Tuple) || length(chain) != 2 || position <= ncarbon(first(chain)) ||
        throw(ArgumentError("Invalid position $position for the ester chain(s)"))
    position
end 

function check_structure!(::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, chain, sn)
    for c in tuplize(chain)
        check_configuration!(bone, c)
        if ispropertyposition(c.doublebond) && any(x -> first(x) == 0x00, c.doublebond)
            throw(ArgumentError("Position `1` cannot contain double bond."))
        end
    end
    (bone, chain, sn, glycerol_class_rs(bone))
end

check_structure!(::Type{<: GlycerophosphoNacylethanolamine}, bone, chain, sn) = 
    (bone, check_configuration!.(Ref(bone), chain), glycerol_class_rs(bone))

function glycerol_class_rs(backbone)
    s = getchaincomponent(backbone)
    rs = collect(map(glycerol_aa_rs, s))
    filter!(!isnothing, rs)
    if length(rs) == 1
        rs = first(rs)
    else
        rs = (reverse!(rs)..., )
    end
    rs
end

glycerol_aa_rs(x) = nothing
glycerol_aa_rs(::Glycerol{T}) where T = T()
glycerol_aa_rs(aa::αAminoAcid) = alpha_rs(aa)

function check_structure!(::Type{<: SphingoBone}, bone, chain, poss, bone_species)
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if bone_species && chain_species || isnothing(bone)
        check_configuration!(bone, fc) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    elseif chain_species
        all(iszero ∘ first, poss) || throw(ArgumentError("Headgroup position cannot be specified at (molecular)specieslevel or without headgroup(s)"))
        check_configuration!(bone, fc) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    elseif ispropertyposition(fc.substituent) || (isnothing(fc.substituent) && any(!iszero ∘ first, poss))
        newfc = false
        if isnothing(fc.substituent)
            sub = Pair{UInt8, AbstractFunctionalGroup}[]
            newfc = true
        elseif XLinkedFunctionalGroup <: eltype(fc.substituent).parameters[end]
            sub = fc.substituent
        else
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
            newfc = true
        end
        if isempty(sub) || all(x -> last(x) == Hydroxy() || last(x) == Hydrogen(), sub)
            if length(poss) > 1 
                push!(sub, first(first(poss)) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(dehydroxygroup(bone))))
                for p in @view poss[2:end]
                    push!(sub, first(p) => XLinkedFunctionalGroup(OLinkage(), NullSubstituent()))
                end
            else
                push!(sub, first(first(poss)) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone)))
            end
            chain_chi = vcat(fc.chirality, poss) 
        else
            chain_chi = fc.chirality
        end
        sort_chainmodification!(sub) 
        # check poss chirality
        if newfc
            fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, chain_chi, fc.isotopiclabel)
        end
        check_configuration!(bone, fc) 
        chi = [get_chainchirality(fc.chirality, p) for (p, _) in poss]
        chi = length(chi) == 1 ? first(chi) : (chi..., )
    elseif ispropertynumber(fc.substituent) || isnothing(fc.substituent)
        newfc = false
        if isnothing(fc.substituent)
            sub = Pair{XLinkedFunctionalGroup, UInt8}[]
            newfc = true
        elseif XLinkedFunctionalGroup <: eltype(fc.substituent).parameters[begin]
            sub = fc.substituent
        else
            sub = convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, fc.substituent)
            newfc = true
        end
        if isempty(sub) || all(x -> first(x) == Hydroxy() || first(x) == Hydrogen(), sub)
            if length(poss) > 1
                m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(dehydroxygroup(bone)))
                id = findfirst(x -> first(x) == m, sub)
                if isnothing(id)
                    push!(sub, m => 0x01)
                else
                    n = last(sub[id])
                    sub[id] = m => (n + 0x01)
                end
                m = XLinkedFunctionalGroup(OLinkage(), NullSubstituent())
                id = findfirst(x -> first(x) == m, sub)
                if isnothing(id)
                    push!(sub, m => UInt8(length(poss) - 1))
                else
                    n = last(sub[id])
                    sub[id] = m => (n + UInt8(length(poss) - 1))
                end
            else
                m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone))
                id = findfirst(x -> first(x) == m, sub)
                if isnothing(id)
                    push!(sub, m => 0x01)
                else
                    n = last(sub[id])
                    sub[id] = m => (n + 0x01)
                end
            end
        end
        sort_chainmodification!(sub)
        headposition = [first(p) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(bone)) for p in poss]
        if newfc
            fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, poss, fc.isotopiclabel)
        end
        check_configuration!(bone, fc; headposition) 
        chi = length(poss) == 1 ? RSChirality() : (RSChirality(), RSChirality())
    else
        throw(ArgumentError("Unknown headgroup position when \"OH\"'s positions are known"))
    end
    if length(poss) == 1
        enpos = first(first(poss))
    elseif length(poss) == 2
        enpos = first(first(poss)) * UInt8(32) + first(last(poss))
    else
        throw(ArgumentError("Invalid headgroup position, \"$poss\""))
    end
    (bone, length(chain) == 1 ? fc : (fc, check_configuration!.(nothing, chain[begin + 1:end]...)), enpos, chi)
end

function check_structure!(::Type{<: MixSphingoBone}, bone, chain, poss, bone_species)
    fc = first(chain)
    chain_species = fc isa CarbonChain{<: Tuple{<: AbstractSPB, <: Acyl}} || (!isnothing(fc.substituent) && any(x -> first(x) isa OxygenAtom, fc.substituent))
    if bone_species && chain_species || isnothing(bone)
        check_configuration!(bone, fc) 
        chi = last.(first.(poss))
    elseif chain_species
        all(iszero ∘ first, vcat(poss...)) || throw(ArgumentError("Headgroup position cannot be specified at specieslevel or without headgroup(s)"))
        check_configuration!(bone, fc) 
        chi = last.(first.(poss))
    elseif ispropertyposition(fc.substituent) || (isnothing(fc.substituent) && any(x -> any(!iszero ∘ first, x), poss))
        newfc = false
        if isnothing(fc.substituent)
            sub = Pair{UInt8, AbstractFunctionalGroup}[]
            newfc = true
        elseif XLinkedFunctionalGroup <: eltype(fc.substituent).parameters[end]
            sub = fc.substituent
        else
            sub = convert(Vector{Pair{UInt8, AbstractFunctionalGroup}}, fc.substituent)
            newfc = true
        end
        length(poss) == length(bone) || throw(ArgumentError("Number of headgroups does not equal number of positions."))
        if isempty(sub) || all(x -> last(x) == Hydroxy() || last(x) == Hydrogen(), sub)
            for (b, p) in zip(bone, poss)
                if length(p) > 1 
                    push!(sub, first(first(p)) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(dehydroxygroup(b))))
                    for x in @view poss[2:end]
                        push!(sub, first(x) => XLinkedFunctionalGroup(OLinkage(), NullSubstituent()))
                    end
                else
                    push!(sub, first(first(p)) => XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b)))
                end
            end
            chain_chi = vcat(fc.chirality, poss...) 
        else
            chain_chi = fc.chirality
        end
        sort_chainmodification!(sub) 
        if newfc
            fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, chain_chi, fc.isotopiclabel)
        end
        check_configuration!(bone, fc) 
        chi = map(poss) do pos 
            if length(pos) == 1
                get_chainchirality(fc.chirality, first(first(pos)))
            else 
                [get_chainchirality(fc.chirality, first(p)) for p in pos]
            end
        end
    elseif ispropertynumber(fc.substituent) || isnothing(fc.substituent)
        newfc = false
        if isnothing(fc.substituent)
            sub = Pair{AbstractFunctionalGroup, UInt8}[]
            newfc = true
        elseif XLinkedFunctionalGroup <: eltype(fc.substituent).parameters[begin]
            sub = fc.substituent
        else
            sub = convert(Vector{Pair{AbstractFunctionalGroup, UInt8}}, fc.substituent)
            newfc = true
        end
        length(poss) == length(bone) || throw(ArgumentError("Number of headgroups does not equal number of positions."))
        if isempty(sub) || all(x -> first(x) == Hydroxy() || first(x) == Hydrogen(), sub)
            for (b, p) in zip(bone, poss)
                if length(p) > 1
                    m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(dehydroxygroup(b)))
                    id = findfirst(x -> first(x) == m, sub)
                    if isnothing(id)
                        push!(sub, m => 0x01)
                    else
                        n = last(sub[id])
                        sub[id] = m => (n + 0x01)
                    end
                    m = XLinkedFunctionalGroup(OLinkage(), NullSubstituent())
                    id = findfirst(x -> first(x) == m, sub)
                    if isnothing(id)
                        push!(sub, m => UInt8(length(p) - 1))
                    else
                        n = last(sub[id])
                        sub[id] = m => (n + UInt8(length(p) - 1))
                    end
                else
                    m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b))
                    id = findfirst(x -> first(x) == m, sub)
                    if isnothing(id)
                        push!(sub, m => 0x01)
                    else
                        n = last(sub[id])
                        sub[id] = m => (n + 0x01)
                    end
                end
            end
        end
        sort_chainmodification!(sub)
        headposition = mapreduce(vcat, zip(bone, poss)) do (b, p) 
            m = XLinkedFunctionalGroup(OLinkage(), dehydroxygroup(b))
            [first(pp) => m for pp in p]
        end
        if newfc
            fc = make_carbonchain(typeof(fc).parameters[begin], fc.carbon, fc.doublebond, sub, vcat(poss...), fc.isotopiclabel)
        end
        check_configuration!(bone, fc; headposition) 
        chi = last.(first.(poss))
    else
        throw(ArgumentError("Unknown head group position when \"OH\"'s positions are known"))
    end
    if any(x -> isa(x, Vector), poss)
        enpos = map(poss) do p 
            p isa Int ? UInt8(p) : length(p) == 1 ? UInt8(first(first(p))) : length(p) == 2 ? (UInt8(first(first(p))) * UInt8(32) + UInt8(first(last(p)))) : throw(ArgumentError("Invalid head group position, \"$p\""))
        end
    else
        enpos = convert(Vector{UInt8}, poss)
    end
    if eltype(chi) <: RSSystem
        chi = convert(Vector{RSSystem}, chi)
    end
    (bone, length(chain) == 1 ? fc : (fc, check_configuration!.(nothing, chain[begin + 1:end]...)), enpos, chi)
end

function get_chainchirality(chirality, p)
    i = findfirst(x -> first(x) == p, chirality)
    isnothing(i) ? RSChirality() : last(chirality[i])
end