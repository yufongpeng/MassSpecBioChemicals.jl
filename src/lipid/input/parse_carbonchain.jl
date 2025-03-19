"""
    parse_carbonchain(Con, bone, echain, schain)

Parse carbon chain string into (`CarbonChain`, snposition code) or ((`CarbonChain`..., ), snposition code)
"""
function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcyl}
    mchain = split_carbonchain(schain)
    length(mchain) > 1 && throw(ArgumentError("Maximal number of chains is 1, got $(length(mchain))"))
    [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: NacylAmine}
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 # sum or single
        [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
    elseif length(mchain) == nchainposition(Con) # NA xx:x/xx:x
        map(echain, mchain) do e, m
            p = parse_position_carbonchain(e, m; force = true)
            first(p) == :start || first(p) == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            last(p)
        end, 0x00
    else
        throw(ArgumentError("Maximal number of chains is $(nchainposition(Con)), got $(length(mchain))"))
    end
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcylEster}
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 
        [last(parse_position_carbonchain(echain, first(mchain); force = true))], 0x00
    elseif echain == (Acyl, Acyl)
        fa = last(parse_position_carbonchain(Acyl, first(mchain); force = true))
        m = match(r"/(\d+)O\(FA[^\s]*(.*)\)", last(mchain))
        if isnothing(m)
            foh = parse_position_carbonchain(Acyl, last(mchain); force = true)
            first(foh) == :inorder || throw(ArgumentError("Invalid chain position, \"$(last(mchain))\""))
            [last(foh), fa], 0x00
        else
            p, a = m
            cchain, pos, ox, sil, mod, sn = parse_singlecarbonchain(a)
            isnothing(ox) || throw(ArgumentError("This fatty acid should have position of functional group, $(last(mchain))"))
            isnothing(sn) || throw(ArgumentError("This fatty acid should not have sn position, $(last(mchain))"))
            # [last(parse_position_carbonchain(Acyl, string(a, ";", p, "OH"); force = true)), fa], parse(UInt8, p)
            [last(parse_position_carbonchain(Acyl, a; force = true)), fa], parse(UInt8, p)
        end
    elseif length(mchain) == nchainposition(Con)
        map(echain, mchain) do e, m
            p = parse_position_carbonchain(e, m; force = true)
            first(p) == :start || first(p) == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            last(p)
        end, 0x00
    else
        throw(ArgumentError("Maximal number of chains is $(nchainposition(Con)), got $(length(mchain))"))
    end
end

function parse_carbonchain(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, echain, schain)
    mchain = split_carbonchain(schain)
    maxsn = nchainposition(Con)
    length(mchain) > maxsn && throw(ArgumentError("Maximal number of chains is $maxsn, got $(length(mchain))"))
    infos = [parse_singlecarbonchain(m) for m in mchain]
    pos = [parse_chainposition(Radyl, info.cchain, info.sn, nothing) for info in infos]
    isn = findall(x -> !isnothing(match(r"^\(sn.*\)$", string(x))), pos)
    if !isempty(isn)
        id = collect(eachindex(pos))
        id = vcat(setdiff!(id, isn), isn)
        infos = infos[id]
        pos = pos[id]
        mchain = mchain[id]
    end
    if length(pos) > 1 && any(==(:inorder), @view pos[begin + 1:end])
        allequal(pos[begin + 1:end]) || throw(ArgumentError("Chain separater should be all \"/\" or \"\\\""))
        pos[begin] = :inorder
    end
    position = map(enumerate(pos)) do (i, c)
        @match c begin 
            :start => 0x00
            :inorder => UInt8(i)
            :noorder => 0x00
            r"\(sn.*\)" => UInt8(findfirst(==(replace(c, "(" => "", ")" => "")), chainposition(Con)))
            _  => throw(ArgumentError("Invalid chain position, \"$(mchain[i])\""))
        end
    end
    if !isempty(isn)
        id = sortperm(position; alg = MergeSort)
        infos = infos[id]
        position = position[id]
        mchain = mchain[id]
    end
    carbonchain = if length(mchain) == maxsn
        [make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for info in infos] # ordered/unordered, ex TG
    elseif length(mchain) == length(echain) # unordered
        all(==(:inorder), pos) && throw(ArgumentError("Expected number of chains in order is $maxsn, got $(length(mchain)); change chain separater to \"_\""))
        [make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for info in infos]
    else
        # remove sn
        # sort 
        nr = length(mchain) - length(isn)
        ne = length(echain) - length(isn)
        snx = ne ÷ nr
        sns = repeat([snx], nr)
        δ = ne - nr * snx
        i = firstindex(sns) # large => small ? 
        while δ > 0
            δ -= 1
            sns[i] += 1
            i += 1
        end
        sns = vcat(sns, repeat([1], length(isn)))
        [sn == 1 ? make_carbonchain(Radyl, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) : 
                        make_carbonchain(ntuple(i -> Radyl, sn), info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod)) for (sn, info) in zip(sns, infos)]
    end

    allunique(filter(!=(0x00), position)) || throw(ArgumentError("Overlapped chain position, \"$schain\""))
    
    # check sn position distribution (sn)
    # :start => 0x01
    # :inorder => vector ID 
    # :noorder => 0x00
    # sn => chainposition ID 
    # base = maxsn + 1
    
    allunique(filter(>(0x00), position)) || throw(ArgumentError("sn position overlapped"))
    id = findall(x -> ncarbon(x) > 0, carbonchain)
    carbonchain = carbonchain[id]
    position = position[id]
    sn = maxsn > 3 ? 0x0000 : 0x00
    bs = convert(typeof(sn), maxsn + 1)
    for p in position
        sn *= bs
        sn += p
    end
    length(carbonchain) == 1 ? first(carbonchain) : ntuple(i -> carbonchain[i], length(carbonchain)), sn
    # prev = prev * base + next
end

function parse_carbonchain(Con::Type{<: Sphingolipid}, bone, echain, schain)
    mchain = split_carbonchain(schain)
    if length(echain) == 3 # ACer...
        # take chain from bone
        if bone isa FattyAcid && ncarbon(bone.chain) > 0 # FA x:y
            echain = first(echain, 2)
        elseif bone isa Tuple && ncarbon(first(bone).chain) > 0
            echain = first(echain, 2)
        else
            throw(ArgumentError("Invalid sphingolipid backbone, $bone"))
        end
        if length(mchain) == 1 # sum level
            pchain = [parse_position_carbonchain(echain, first(mchain))]
        elseif length(mchain) == 2 # species level
            pchain = [parse_position_carbonchain(C, m) for (C, m) in zip(first(echain, 2), mchain)]
        else
            throw(ArgumentError("Invalid number of chains, $(length(mchain))"))
        end
    elseif length(mchain) == length(echain)
        pchain = [parse_position_carbonchain(C, m) for (C, m) in zip(echain, mchain)]
    else # sum level
        pchain = [parse_position_carbonchain(echain, first(mchain))]
    end
    for (i, c) in enumerate(pchain)
        @match first(c) begin 
            :start => 0x00
            :inorder => UInt8(i)
            _ => throw(ArgumentError("Invalid chain position, \"$(mchain[i])\""))
        end
    end
    ntuple(i -> last(pchain[i]), length(pchain)), 0x00
end
function parse_carbonchain(Con::Type{<: Sterol}, bone, pos, chain, sn)
    throw(ArgumentError("`parse_carbonchain` not implemented for `Sterol`"))
end

function parse_carbonchain(Con::Type{<: Prenol}, bone, pos, chain, sn)
    throw(ArgumentError("`parse_carbonchain` not implemented for `Prenol`"))
end

"""
    split_carbonchain(schain)

Split carbon chain string by "/" and "_"; return a vector of carbon chain string
"""
function split_carbonchain(schain)
    s1 = split(schain, "/")
    fs = popfirst!(s1)
    chains = String[]
    s2 = split(fs, "_")
    push!(chains, popfirst!(s2))
    for s in s2
        push!(chains, string("_", s))
    end
    for s in s1
        s2 = split(string("/", s), "_")
        push!(chains, popfirst!(s2))
        for s in s2
            push!(chains, string("_", s))
        end
    end
    chains
end

"""
    parse_position_carbonchain(T, mchain; position = nothing, force = false)

Parse carbon chain string into snposition symbol => `CarbonChain`

`force`: force carbon chain to be `T`
"""
function parse_position_carbonchain(T, mchain; position = nothing, force = false)
    # pos => chain
    # pos:
    # :start start
    # :inorder ordered
    # :noorder unordered
    # number for ACer
    # sn-\d
    # cchain_sil = match(r"^([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?", mchain)
    # isnothing(cchain_sil) && throw(ArgumentError("Invalid fattyacyl chain, \"$mchain\""))
    # m = cchain_sil.match
    # cchain, pos, sil = cchain_sil
    # mchain = mchain[m.offset + m.ncodeunits + 1:end]
    # mod = collect(eachmatch(r"((?:;(([^)(;/_]*[\(\[][^)(]*+(?:(?3)[^)(]*)*+[\)\]])?[^)(;/_]*)))", mchain))
    #                           ((?:;[^)(;/_]*([\(\[][^)(]*+(?:(?2)[^)(]*)*+[\)\]])?[^)(;/_]*)*)
    # sn, = match(r"(\(sn-*\d*'*\))?$", mchain)
    # r"((?:;(([^)(\[\];/_]*\([^)(]*+(?:(?3)[^)(]*)*+\))?([^)(\[[^\[\]]*+(?:(?3)[^\[\]]*)*+\])?[^)(;/_]*)))"
    cchain, pos, ox, sil, mod, sn = parse_singlecarbonchain(mchain)
    parse_chainposition(T, cchain, sn, position) => make_carbonchain(T, cchain, pos, ox, sil, split_chainmodification(mod); force)
end

function parse_singlecarbonchain(mchain) 
    cchain, pos, ox, sil, mod, _, sn = (isnothing(x) ? nothing : isempty(x) ? nothing : x for x in match(REGEX[:chain], mchain))
    (; cchain, pos, ox, sil, mod, sn)
end

"""
    parse_chainposition(T, mc, sn, position)

Parse chain position into snposition symbol
"""
function parse_chainposition(::Type{<: Radyl}, mc, sn, position)
    !isnothing(position) ? position : 
    !isnothing(sn) ? sn : 
    startswith(mc, r"\s") ? :start : 
    startswith(mc, "/") ? :inorder :
    startswith(mc, "_") ? :noorder : throw(ArgumentError("Invalid fattyacyl chain, \"$mc\""))
end

function parse_chainposition(::Type{<: SPB}, mc, sn, position)
    isnothing(sn) || throw(ArgumentError("SPB does not have sn position"))
    !isnothing(position) ? position : 
    startswith(mc, r"\s") ? :start : 
    startswith(mc, "/") ? :inorder :
    startswith(mc, "_") ? :noorder : throw(ArgumentError("Invalid fattyacyl chain, \"$mc\""))
end

function parse_chainposition(::Type{T}, mc, sn, position) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_chainposition` not implemented for `$T"))
end

function parse_chainposition(T::Tuple, mc, sn, position)
    if length(T) == 1
        parse_chainposition(first(T), mc, sn, position)
    elseif SPB in T
        parse_chainposition(SPB, mc, sn, position)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        parse_chainposition(T[i], mc, sn, position)
    else
        parse_chainposition(Radyl, mc, sn, position)
    end
end

"""
    make_carbonchain(T, cchain::AbstractString, pos, ox, sil, mod; force = false)
    make_carbonchain(T, carbon::Number, doublebond, substituent, isotopiclabel)

Construct `CarbonChain` with parsed information
"""
function make_carbonchain(T::Tuple, cchain::AbstractString, pos, ox, sil, mod; force = false)
    if length(T) == 1
        make_carbonchain(first(T), cchain, pos, ox, sil, mod; force)
    elseif SPB in T
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(SPB, cchain, pos, ox, sil, mod)...)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(T[i], cchain, pos, ox, sil, mod)...)
    elseif force
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(Radyl, cchain, pos, ox, sil, mod)...)
    else
        n, rad, chain = match(r"([d,t,e]?)([P,O]?-?)(\d+:\d+.*)", cchain)
        if isnothing(rad) 
            Chain = Tuple{ntuple(i -> Acyl, length(T))...}
        else
            n = @match n begin
                ""  => 1
                "d" => 2
                "t" => 3
                "e" => 4
                _   => throw(ArgumentError("Invalid representation of numbers of alkyl or alkenyl chain"))
            end
            echain = [Acyl for i in eachindex(T)]
            rad = @match rad begin
                ""   => Acyl
                "O-" => Alkyl
                "P-" => Alkenyl
                _    => throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
            end
            echain[begin:begin + n - 1] .= rad
            Chain = Tuple{echain...}
        end
        CarbonChain{Chain}(parse_carbonchainbody(Radyl, chain, pos, ox, sil, mod)...)
    end
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: Radyl}
    force && return CarbonChain{T}(parse_carbonchainbody(T, cchain, pos, ox, sil, mod)...)
    n, rad, chain = match(r"([d,t,e]?)([P,O]?-?)(\d+:\d+.*)", cchain)
    isempty(n) || throw(ArgumentError("This should be a single fattyacyl chain, \"$cchain\""))
    Chain = @match rad begin
        ""   => Acyl
        "O-" => Alkyl
        "P-" => Alkenyl
        _    => throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
    end
    Chain <: T || throw(ArgumentError("Fattyacyl chain does not match to class"))
    CarbonChain{Chain}(parse_carbonchainbody(Chain, chain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: SPB}
    rad, chain = match(r"([d,t,e]?[P,O]?-?)(\d+:\d+.*)", cchain)
    isempty(rad) || throw(ArgumentError("Invalid fattyacyl chain, \"$cchain\""))
    CarbonChain{T}(parse_carbonchainbody(T, chain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchain` not implemented for `$T`"))
end

make_carbonchain(T, carbon::C, doublebond, substituent, isotopiclabel = nothing) where {C <: Number} = CarbonChain{T}(UInt8(carbon), uint8ize(doublebond), uint8ize(substituent), isotopiclabel)

"""
    parse_carbonchainbody(T, cchain, pos, ox, sil, mod; force = false) -> (carbon, doublebond, substituent, isotopiclabel)

Parse carbon chain info strings into correct type and format for `CarbonChain` construction
"""
function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: Radyl}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    # mod
    sub = parse_chainmodification(T, mod, ox)
    isnothing(pos) && return cb, db, sub, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid double bond configuration, \"$(string(ps))\""))
        @inbounds db[i] = x + e
    end
    return cb, db, sub, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: SPB}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    sub = parse_chainmodification(T, mod, ox)
    isnothing(pos) && return cb, db, sub, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid double bond configuration, \"$(string(ps))\""))
        @inbounds db[i] = x + e
    end
    return cb, db, sub, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchainbody` not implemented for `$T`"))
end
