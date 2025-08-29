"""
    makechemical(m; sil)
    makechemical(::Type{T}, m::T; linkage = nothing, config = nothing, sil = nothing)
    makechemical(::Type{T}, ms...; linkage = nothing, config = nothing, sil = nothing)

Construct a chemical of type `T`. If multiple chemicals are given, `T` must be chained chemicals and the first argument contains chemicals and the second argument is linkage information. 

* `m` or `ms`: any chemical(s).
* `linkage`: vector of linkages. The structure depends on chemical type `T`.
* `sil`: stabel isotopic labeling. 
"""
function makechemical(m; sil)
    isnothing(sil) ? m : IsotopiclabeledChemical(m, sil)
end

function makechemical(::Type{T}, m::T; linkage = nothing, config = nothing, sil = nothing) where T
    isnothing(sil) ? m : IsotopiclabeledChemical(m, sil)
end

function makechemical(::Type{T}, ms...; linkage = nothing, config = nothing, sil = nothing) where T
    chemicals = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = getchainlinkage(ms[i]), config = getchainconfig(ms[i]), sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    length(last(chemicals)) == 1 || throw(ArgumentError("Last chemical should be a single chemical"))
    ids = zeros(Int, length(chemicals) - 1)
    i = lastindex(chemicals) - 1
    j = lastindex(chemicals)
    while i >= firstindex(chemicals)
        ids[i] = j
        if length(chemicals[i]) == 1
            j = i
            i -= 1
        else
            i -= 1
        end
    end
    requirelinkage(T) || return chainedchemical(T, chemicals; config)
    lso = collect(Any, map(x -> collect(transformlinkage(T, x)), ms))
    ls = Any[makelinkage(T, last(getchaincomponent(a)), chemicals[b]) for (a, b) in zip(chemicals[begin:end - 1], ids)]
    for i in eachindex(ls)
        if !isnothing(linkage) && !isnothing(linkage[i])
            ls[i] = linkage[i]
        else 
            l = length(lso[i]) > 0 ? last(lso[i]) : nothing
            if isnothing(l) || isnulllinkage(T, l)
                continue
            elseif isnulllinkage(T, first(l))
                ls[i] = first(ls[i]) => last(l)
            elseif isnulllinkage(T, last(l))
                ls[i] = first(l) => last(ls[i])
            else
                ls[i] = l
            end
        end
    end
    chainedchemical(T, chemicals; linkage = ls, config)
end

"""
    concatchemical(m; sil)
    concatchemical(::Type{T}, m::T; linkage = nothing, config = nothing, sil = nothing)
    concatchemical(::Type{T}, ms...; linkage = nothing, config = nothing, sil = nothing)

Construct a chemical of type `T`. Each chemical are considered as linearly chained chemicals and all chemicals are concatenated into a longer linearly chained chemicals `T`. 

* `m` or `ms`: any chemical(s).
* `linkage`: vector of linkages. The structure depends on chemical type `T`.
* `sil`: stabel isotopic labeling. 
"""
function concatchemical(m::T; linkage = nothing, config = nothing, sil) where T
    concatchemical(T, m; linkage, config, sil)
end

function concatchemical(::Type{T}, m::T; linkage = nothing, config = nothing, sil = nothing) where T
    isnothing(sil) && isnothing(linkage) && return m
    ischainedchemical(m) || return IsotopiclabeledChemical(m, s)
    ms = getchaincomponent(m)
    if requirelinkage(T)
        ls = getchainlinkage(m)
        if length(ls) < length(ms) - 1
            throw(ArgumentError("`linkage` for $(m) must have $(length(ms)) or $(length(ms) - 1) elements"))
        end
    end
    chemicals = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = getchainlinkage(ms[i]), config = getchainconfig(ms[i]), sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    requirelinkage(T) || return chainedchemical(T, chemicals; config)
    chainedchemical(T, chemicals; linkage = isnothing(linkage) ? ls : [first(ls, length(chemicals) - 1)..., first(linkage)], config)
end

function concatchemical(::Type{T}, ms...; linkage = nothing, config = nothing, sil = nothing) where T
    # T is chain ?
    ischainedchemical(T) || throw(ArgumentError("`$T` must be a linearly chained chemical type."))
    mss = collect(Any, map(collect ∘ getchaincomponent, ms))
    if requirelinkage(T)
        lss = collect(Any, map(x -> collect(transformlinkage(T, x)), ms))
        lastlk = false
        if !isnothing(linkage)
            if length(linkage) == length(mss) - 1
                lastlk = false
            elseif length(linkage) == length(mss) && isnothing(last(linkage))
                lastlk = false
            elseif length(linkage) == length(mss)
                lastlk = true
            else
                throw(ArgumentError("`linkage` must have $(length(mss)) or $(length(mss) - 1) elements"))
            end
        end
        # linkage ?
        for i in eachindex(mss)
            if length(lss[i]) < length(mss[i]) - 1
                throw(ArgumentError("`linkage` for $(ms[i]) must have $(length(mss[i])) or $(length(mss[i]) - 1) elements"))
            elseif i == lastindex(mss)
                if lastlk
                    lss[i] = [first(lss[i], length(mss[i]) - 1)..., last(linkage)]
                end
            elseif !isnothing(linkage) && !isnothing(linkage[i])
                lss[i] = [first(lss[i], length(mss[i]) - 1)..., linkage[i]]
            elseif length(lss[i]) >= length(mss[i])
                ls = lss[i][length(mss[i])]
                ni = findfirst(x -> length(x) == 1, mss[i + 1])
                nm = mss[i + 1][ni]
                if isnulllinkage(T, ls)
                    lss[i] = [first(lss[i], length(mss[i]) - 1)..., makelinkage(T, last(mss[i]), nm)]
                elseif isnulllinkage(T, first(ls))
                    lss[i] = [first(lss[i], length(mss[i]) - 1)..., first(makelinkage(T, last(mss[i]), nm)) => last(ls)]
                elseif isnulllinkage(T, last(ls))
                    lss[i] = [first(lss[i], length(mss[i]) - 1)..., first(ls) => last(makelinkage(T, last(mss[i]), nm))]
                else
                    lss[i] = first(lss[i], length(mss[i]))
                end
            else
                ni = findfirst(x -> length(x) == 1, mss[i + 1])
                nm = mss[i + 1][ni]
                lss[i] = [lss[i]..., makelinkage(T, last(mss[i]), first(mss[i + 1]))]
            end
        end
    end
    ms = vcat(mss...)
    chemicals = if isnothing(sil)
        (ms..., )
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = getchainlinkage(ms[i]), config = getchainconfig(ms[i]), sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    requirelinkage(T) ? chainedchemical(T, chemicals; linkage = vcat(lss...), config) : chainedchemical(T, chemicals; config)
end

"""
"""
chainedchemical(::Type{<: T}, chemicals; kwargs...) where {T <: AbstractChemical} = throw(MethodError(chainedchemical, (T, chemicals, kwargs...)))
chainedchemical(::Type{<: ChainedChemical}, chemicals; linkage = nothing, config = nothing, kwargs...) = ChainedChemical(chemicals, linkage, config)
chainedchemical(::Type{<: DehydratedChemical}, chemicals; linkage = nothing, config = nothing, kwargs...) = DehydratedChemical(chemicals, linkage, config)

"""
"""
chiralchemical(::Type{C}, ::Type{I}, x...; kwargs...) where {C, I} = C{I}(x...; kwargs...)

"""
"""
function deletechemicalat(chemicals::T, i; linkage = nothing, config = nothing, kwargs...) where T
    i = sort!(unique(i))
    cc = getchaincomponent(chemicals)
    ls = getchainlinkage(chemicals)
    cs = getchainconfig(chemicals)
    if isnothing(linkage) && !isnothing(ls)
        vl = collect(ls)
        del = Int[]
        for j in i
            if j == firstindex(vl)
                push!(del, j)
            elseif j == lastindex(vl) == lastindex(cc)
                # keep last linkage 
                push!(del, j)
                vl[j - 1] = first(vl[j - 1]) => lk(0x00)
            elseif j == lastindex(cc)
                push!(del, j - 1)
            else
                push!(del, j - 1)
                vl[j] = first(vl[j - 1]) => last(vl[j])
            end
        end
        deleteat!(vl, del)
        linkage = (vl..., )
    elseif isnothing(linkage)
        linkage = deepcopy(ls)
    end
    if isnothing(config) && !isnothing(cs)
        config = deepcopy(cs)
        del = Int[] 
        for (k, v) in enumerate(config) 
            first(v) in i && push!(del, k)
        end
        deleteat!(config, del)
    elseif isnothing(config)
        config = deepcopy(cs)
    end
    cc = cc[setdiff(eachindex(cc), i)]
    makechemical(T, cc...; linkage, config)
end

deletechemicalat(chemicals::T, i::Number; linkage = nothing, config = nothing, kwargs...) where T = 
    deletechemicalat(chemicals, [i]; linkage, config, kwargs...)