"""
    makechemical(m; sil)
    makechemical(::Type{T}, m::T; linkage = nothing, sil = nothing)
    makechemical(::Type{T}, ms...; linkage = nothing, sil = nothing)

Construct a chemical of type `T`. If multiple chemicals are given, `T` must be chained chemicals and the first argument contains chemicals and the second argument is linkage information. 

* `m` or `ms`: any chemical(s).
* `linkage`: vector of linkages. The structure depends on chemical type `T`.
* `sil`: stabel isotopic labeling. 
"""
function makechemical(m; sil)
    isnothing(sil) ? m : IsotopiclabeledChemical(m, s)
end

function makechemical(::Type{T}, m::T; linkage = nothing, sil = nothing) where T
    isnothing(sil) ? m : IsotopiclabeledChemical(m, s)
end

function makechemical(::Type{T}, ms...; linkage = nothing, sil = nothing) where T
    chemicals = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
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
    T(chemicals, ls)
end

"""
    concatchemical(m; sil)
    concatchemical(::Type{T}, m::T; linkage = nothing, sil = nothing)
    concatchemical(::Type{T}, ms...; linkage = nothing, sil = nothing)

Construct a chemical of type `T`. Each chemical are considered as linearly chained chemicals and all chemicals are concatenated into a longer linearly chained chemicals `T`. 

* `m` or `ms`: any chemical(s).
* `linkage`: vector of linkages. The structure depends on chemical type `T`.
* `sil`: stabel isotopic labeling. 
"""
function concatchemical(m::T; linkage = nothing, sil) where T
    concatchemical(T, m; linkage, sil)
end

function concatchemical(::Type{T}, m::T; linkage = nothing, sil = nothing) where T
    isnothing(sil) && isnothing(linkage) && return m
    ischainedchemical(m) || return IsotopiclabeledChemical(m, s)
    ms = getchaincomponent(m)
    ls = getchainlinkage(m)
    if length(ls) < length(ms) - 1
        throw(ArgumentError("`linkage` for $(m) must have $(length(ms)) or $(length(ms) - 1) elements"))
    end
    chemicals = if isnothing(sil)
        ms
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    T(chemicals, isnothing(linkage) ? ls : [first(ls, length(chemicals) - 1)..., first(linkage)])
end

function concatchemical(::Type{T}, ms...; linkage = nothing, sil = nothing) where T
    # T is chain ?
    mss = collect(Any, map(collect ∘ getchaincomponent, ms))
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
    ms = vcat(mss...)
    ls = vcat(lss...)
    chemicals = if isnothing(sil)
        (ms..., )
    else
        ntuple(length(ms)) do i
            j = findfirst(x -> i == first(x), sil)
            isnothing(j) ? ms[i] : 
            ischainedchemical(ms[i]) ? makechemical(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sil[j])) : 
            IsotopiclabeledChemical(ms[i], last(sil[j]))
        end
    end
    T(chemicals, ls)
end
