"""
    parse_carbonchain(Con, bone, echain, schain)

Parse carbon chain string into (`CarbonChain`, snposition code, chirality(s)) or ((`CarbonChain`..., ), snposition code, chirality(s))
"""
function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcyl}
    mchain = split_carbonchain(schain)
    length(mchain) > 1 && throw(ArgumentError("Maximal number of chains is 1, got $(length(mchain))"))
    [make_carbonchain(echain, parse_singlecarbonchain(first(mchain)); force = true)], 0x00
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: NacylAmine}
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 # sum or single
        [make_carbonchain(echain, parse_singlecarbonchain(first(mchain)); force = true)], 0x00
    elseif length(mchain) == nchainposition(Con) # NA xx:x/xx:x
        map(echain, mchain) do e, m
            info = parse_singlecarbonchain(m)
            p = parse_chainposition(e, info)
            p == :start || p == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            make_carbonchain(e, info; force = true)
        end, 0x00
    else
        throw(ArgumentError("Maximal number of chains is $(nchainposition(Con)), got $(length(mchain))"))
    end
end

function parse_carbonchain(Con::Type{T}, bone, echain, schain) where {T <: FattyAcylEster}
    # chiral
    mchain = split_carbonchain(schain)
    if length(mchain) == 1 
        [make_carbonchain(echain, parse_singlecarbonchain(first(mchain)); force = true)], 0x00
    elseif echain == (Acyl, Acyl)
        fa = make_carbonchain(Acyl, parse_singlecarbonchain(first(mchain)); force = true)
        m = match(r"^/(\d+)O\(FA[^\s]*(.*)\)(\[[RS]\])?$", last(mchain))
        if isnothing(m)
            # /cb:db;O
            hfa = parse_singlecarbonchain(last(mchain))
            parse_chainposition(Acyl, hfa) == :inorder || throw(ArgumentError("Invalid chain position, \"$(last(mchain))\""))
            [make_carbonchain(Acyl, hfa; force = true), fa], 0x00
        else
            p, a, c = m 
            # isnothing(c) || throw(ArgumentError("Hydroxy fatty acid is not chiral."))
            info = parse_singlecarbonchain(a)
            isnothing(info.ox) || throw(ArgumentError("This fatty acid should have position of functional group, $(last(mchain))"))
            isnothing(info.sn) || throw(ArgumentError("This fatty acid should not have sn position, $(last(mchain))"))
            hfa = parse_singlecarbonchain(string(a, ";", p, "OH", isnothing(c) ? "" : c))
            # delete pOH in output 
            # hfa = last(parse_carbonchain_info(Acyl, a; force = true))
            [make_carbonchain(Acyl, hfa; force = true), fa], parse(UInt8, p)
        end
    elseif length(mchain) == nchainposition(Con)
        map(echain, mchain) do e, m
            info = parse_singlecarbonchain(m)
            p = parse_chainposition(e, info)
            p == :start || p == :inorder || throw(ArgumentError("Invalid chain position, \"$m\""))
            make_carbonchain(e, info; force = true)
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
    pos = [parse_chainposition(Radyl, info) for info in infos]
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
        [make_carbonchain(Radyl, info) for info in infos] # ordered/unordered, ex TG
    elseif length(mchain) == length(echain) # unordered
        all(==(:inorder), pos) && throw(ArgumentError("Expected number of chains in order is $maxsn, got $(length(mchain)); change chain separater to \"_\""))
        [make_carbonchain(Radyl, info) for info in infos]
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
        [sn == 1 ? make_carbonchain(Radyl, info) : 
                        make_carbonchain(ntuple(i -> Radyl, sn), info) for (sn, info) in zip(sns, infos)]
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
    # check chiral for parse_chainchirality
    # ch_pos = config_chain(Con)
    # if all(==(:inorder), pos)
    #     ch = map(enumerate(infos)) do (i, info) 
    #         parse_chainchirality(Radyl, info; chiral = UInt8(i) in ch_pos)
    #     end
    #     filter!(!=(AChirality()), ch)
    # else
    #     ch = map(ch_pos) do _ 
    #         RSChirality()
    #     end
    # end
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
            pchain = (make_carbonchain(echain, parse_singlecarbonchain(first(mchain))), )
        elseif length(mchain) == 2 # species level
            pchain = (make_carbonchain(first(echain), parse_singlecarbonchain(first(mchain))), 
                    make_carbonchain(echain[begin + 1], parse_singlecarbonchain(last(mchain))))
        else
            throw(ArgumentError("Invalid number of chains, $(length(mchain))"))
        end
    elseif length(mchain) == length(echain) == 1 
        pchain = (make_carbonchain(first(echain), parse_singlecarbonchain(first(mchain))), )
    elseif length(mchain) == length(echain) 
        pchain = (make_carbonchain(first(echain), parse_singlecarbonchain(first(mchain))), 
                    make_carbonchain(last(echain), parse_singlecarbonchain(last(mchain))))
    else # sum level
        pchain = (make_carbonchain(echain, parse_singlecarbonchain(first(mchain))), )
    end
    # for (i, c) in enumerate(pchain)
    #     @match first(c) begin 
    #         :start => 0x00
    #         :inorder => UInt8(i)
    #         _ => throw(ArgumentError("Invalid chain position, \"$(mchain[i])\""))
    #     end
    # end
    pchain, 0x00
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
    parse_carbonchain_info(T, mchain; position = nothing, force = false)

Parse carbon chain string into snposition symbol => `CarbonChain`

`force`: force carbon chain to be `T`
"""
function parse_carbonchain_info(T, mchain; position = nothing, force = false )
    # pos => chain
    # pos:
    # :start start
    # :inorder ordered
    # :noorder unordered
    # number for ACer
    # sn-\d
    # cchain_sil = match(r"^([\s,/,_](\[[RS]\])?[d,t,e]?[P,O,E]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?", mchain)
    # isnothing(cchain_sil) && throw(ArgumentError("Invalid fattyacyl chain, \"$mchain\""))
    # m = cchain_sil.match
    # cchain, pos, sil = cchain_sil
    # mchain = mchain[m.offset + m.ncodeunits + 1:end]
    # mod = collect(eachmatch(r"((?:;(([^)(;/_]*[\(\[][^)(]*+(?:(?3)[^)(]*)*+[\)\]])?[^)(;/_]*)))", mchain))
    #                           ((?:;[^)(;/_]*([\(\[][^)(]*+(?:(?2)[^)(]*)*+[\)\]])?[^)(;/_]*)*)
    # sn, = match(r"(\(sn-*\d*'*\))?$", mchain)
    # r"((?:;(([^)(\[\];/_]*\([^)(]*+(?:(?3)[^)(]*)*+\))?([^)(\[[^\[\]]*+(?:(?3)[^\[\]]*)*+\])?[^)(;/_]*)))"
    info = parse_singlecarbonchain(mchain)
    parse_chainposition(T, info; position) => make_carbonchain(T, info; force)
end

function parse_singlecarbonchain(mchain) 
    sep, rad, cchain, pos, ox, sil, mod, _, sn = (isnothing(x) ? nothing : isempty(x) ? nothing : x for x in match(REGEX[:chain], mchain))
    (; sep, rad, cchain, pos, ox, sil, mod, sn)
end   

"""
    parse_chainposition(T, sep, sn; position)

Parse chain position into snposition symbol
"""
parse_chainposition(T, info; position = nothing) = parse_chainposition(T, info.sep, info.sn; position)
function parse_chainposition(::Type{<: Radyl}, sep, sn; position = nothing)
    !isnothing(position) ? position : 
    !isnothing(sn) ? sn : 
    @match sep begin
        r"\s"   => :start
        "/"     => :inorder 
        "_"     => :noorder
        _       => throw(ArgumentError("Invalid fattyacyl chain, \"$sep\""))
    end
end

function parse_chainposition(::Type{<: AbstractSPB}, sep, sn; position = nothing)
    isnothing(sn) || throw(ArgumentError("SPB does not have sn position"))
    !isnothing(position) ? position : 
        @match sep begin
        r"\s"   => :start
        "/"     => :inorder 
        "_"     => :noorder
        _       => throw(ArgumentError("Invalid fattyacyl chain, \"$sep\""))
    end
end

function parse_chainposition(::Type{T}, sep, sn; position = nothing) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_chainposition` not implemented for `$T"))
end

function parse_chainposition(T::Tuple, sep, sn; position = nothing)
    if length(T) == 1
        parse_chainposition(first(T), sep, sn; position)
    elseif SPB in T
        parse_chainposition(SPB, sep, sn; position)
    elseif SulfoSPB in T
        parse_chainposition(SulfoSPB, sep, sn; position)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        parse_chainposition(T[i], sep, sn; position)
    else
        parse_chainposition(Radyl, sep, sn; position)
    end
end

# """
# """
# parse_chainchirality(T, info; chiral = false) = parse_chainchirality(T, info.ch; chiral)
# function parse_chainchirality(::Tuple, ch; chiral = false) 
#     isnothing(ch) || @warn "Sum carbon chain is not chiral"
#     AChirality() 
# end

# function parse_chainchirality(::Type{T}, ch; chiral = false) where {T <: Radyl}
#     if chiral 
#         parse_rschirality(ch)
#     else
#         isnothing(ch) || @warn "This carbon chain is not chiral"
#         AChirality()
#     end
# end

# function parse_chainchirality(::Type{T}, ch; chiral = false) where {T <: AbstractSPB}
#     isnothing(ch) || @warn "Sphingoid base is not chiral"
#     AChirality()
# end

function parse_rschirality(ch)
    @match ch begin
        "[R]"   => RChirality()
        "[S]"   => SChirality()
        _       => RSChirality()
    end
end
"""
    make_carbonchain(T, cchain::AbstractString, pos, ox, sil, mod; force = false)
    make_carbonchain(T, carbon::Number, doublebond, substituent, isotopiclabel)

Construct `CarbonChain` with parsed information
"""
make_carbonchain(T, info; force = false) = make_carbonchain(T, info.rad, info.cchain, info.pos, info.ox, info.sil, split_chainmodification(info.mod); force)
function make_carbonchain(T::Tuple, rad, cchain::AbstractString, pos, ox, sil, mod; force = false)
    length(T) == 1 && return make_carbonchain(first(T), rad, cchain, pos, ox, sil, mod; force)
    if SPB in T
        isnothing(rad) || throw(ArgumentError("Invalid fattyacyl chain, \"$rad$cchain\""))
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(SPB, cchain, pos, ox, sil, mod)...)
    elseif SulfoSPB in T
        isnothing(rad) || throw(ArgumentError("Invalid fattyacyl chain, \"$rad$cchain\""))
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(SulfoSPB, cchain, pos, ox, sil, mod)...)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(T[i], cchain, pos, ox, sil, mod)...)
    elseif force
        CarbonChain{Tuple{T...}}(parse_carbonchainbody(Radyl, cchain, pos, ox, sil, mod)...)
    elseif isnothing(rad) 
        CarbonChain{Tuple{ntuple(i -> Acyl, length(T))...}}(parse_carbonchainbody(Radyl, cchain, pos, ox, sil, mod)...)
    else
        n, rad = match(r"([d,t,e]?)([P,O,E]?-?)", rad)
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
            "P-" => Alkenyl{ZConfiguration}
            "E-" => Alkenyl{EConfiguration}
            _    => throw(ArgumentError("Invalid fattyacyl chain, \"$n$rad$cchain\""))
        end
        echain[begin:begin + n - 1] .= rad
        Chain = Tuple{echain...}
        CarbonChain{Chain}(parse_carbonchainbody(Radyl, cchain, pos, ox, sil, mod)...)
    end
end

function make_carbonchain(::Type{T}, rad, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: Radyl}
    if force 
        # check n rad ch
        return CarbonChain{T}(parse_carbonchainbody(T, cchain, pos, ox, sil, mod)...)
    end
    if isnothing(rad)
        rad = "" 
    else
        n, rad = match(r"([d,t,e]?)([P,O,E]?-?)", rad)
        isempty(n) || throw(ArgumentError("This should be a single fattyacyl chain, \"$n$rad$cchain\""))
    end  
    Chain = @match rad begin
        ""   => Acyl
        "O-" => Alkyl
        "P-" => Alkenyl{ZConfiguration}
        "E-" => Alkenyl{EConfiguration}
        _    => throw(ArgumentError("Invalid fattyacyl chain, \"$rad$cchain\""))
    end
    Chain <: T || throw(ArgumentError("Fattyacyl chain does not match to class"))
    CarbonChain{Chain}(parse_carbonchainbody(Chain, cchain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, rad, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: AbstractSPB}
    isnothing(rad) || throw(ArgumentError("Invalid sphingoid base, \"$rad$cchain\""))
    CarbonChain{T}(parse_carbonchainbody(T, cchain, pos, ox, sil, mod)...)
end

function make_carbonchain(::Type{T}, rad, cchain::AbstractString, pos, ox, sil, mod; force = false) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchain` not implemented for `$T`"))
end

make_carbonchain(T, carbon::C, doublebond, substituent, chirality, isotopiclabel = nothing) where {C <: Number} = CarbonChain{T}(UInt8(carbon), uint8ize(doublebond), uint8ize(substituent), uint8ize(chirality), isotopiclabel)

"""
    parse_carbonchainbody(T, cchain, pos, ox, sil, mod; force = false) -> (carbon, doublebond, substituent, isotopiclabel)

Parse carbon chain info strings into correct type and format for `CarbonChain` construction
"""
function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod = false) where {T <: Radyl}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    # mod
    sub, chi = parse_chainmodification(T, mod, ox)
    # if !isnothing(ch) 
    #     ch = ch == "R" ? RChirality() : ch == "S" ? SChirality() : ch == "U" ? RSChirality() : throw(ArgumentError("Invalid chirality for carbon chain."))
    #     if isnothing(chi) 
    #         chi = Pair{UInt8, RSSystem}[0xff => ch]
    #     else
    #         push!(chi, 0xff => ch)
    #     end
    # end
    if !isnothing(pos) 
        db = parse_geometricconfig(db, pos)
        if first(last(db)) == cb 
            throw(ArgumentError("The last position cannot contain double bond."))
        end
    end
    cb, db, sub, chi, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod = false) where {T <: AbstractSPB}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    sub, chi = parse_chainmodification(T, mod, ox)
    if !isnothing(pos) 
        db = parse_geometricconfig(db, pos)
        if first(last(db)) == cb 
            throw(ArgumentError("The last position cannot contain double bond."))
        end
    end
    cb, db, sub, chi, nothing
end

function parse_carbonchainbody(::Type{T}, cchain, pos, ox, sil, mod) where {T <: AbstractSTRing}
    throw(ArgumentError("`parse_carbonchainbody` not implemented for `$T`"))
end

function parse_geometricconfig(db, pos)
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    dbm = Dict{UInt8, GeometricConfiguration}()
    for p in ps
        x, e = p
        x = parse(UInt8, x)
        haskey(dbm, x) && throw(ArgumentError("Multiple doublebond at position `$(Int(x))`."))
        e = isnothing(e) ? EZConfiguration() : e == "Z" ? ZConfiguration() : e == "E" ? EConfiguration() : throw(ArgumentError("Invalid double bond configuration, \"$(string(ps))\""))
        v = get!(dbm, x, e)
        # if v == e 
        #     continue 
        # elseif e == EZConfiguration()
        #     continue 
        # elseif v == EZConfiguration()
        #     dbm[x] = e 
        # else 
        #     throw(ArgumentError("Conflict double bond configuration at position `$(Int(x))`"))
        # end
    end
    dbp = Pair{UInt8, GeometricConfiguration}[UInt8(k) => v for (k, v) in dbm]
    n = db - length(dbp)
    if n > 0 
        @warn "Positions are not assigned for some double bonds."
        dbp = vcat(repeat([0x00 => NoEZConfiguration()], n), dbp)
    elseif n < 0 
        @warn "Actual number of double bonds are larger than that indicated in `cb:db` notation"
    end
    sort!(dbp; by = first)
end
