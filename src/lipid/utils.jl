uint8ize(x) = x
uint8ize(x::UInt8) = x
uint8ize(x::Number) = UInt8(x)
uint8ize(x::Pair{M, UInt8}) where M = x 
uint8ize(x::Pair{UInt8, M}) where M = x 
uint8ize(x::Pair{UInt8, UInt8}) = x 
uint8ize(x::Pair{M, <: Number}) where M = first(x) => uint8ize(last(x))
uint8ize(x::Pair{<: Number, M}) where M = uint8ize(first(x)) => last(x)
uint8ize(x::Pair{<: Number, <: Number}) = uint8ize(first(x)) => uint8ize(last(x))
uint8ize(::Nothing) = nothing
uint8ize(x::Vector{UInt8}) = x
uint8ize(x::Vector{<: Number}) = map(uint8, x)
# uint8ize(x::Vector{<: Pair{M, UInt8}}) where M = x
uint8ize(x::Vector{<: Pair{M, UInt8} where M}) = x
uint8ize(x::Vector{<: Pair{UInt8, M} where M}) = x
uint8ize(x::Vector{<: Pair{UInt8, UInt8}}) = x
uint8ize(x::Vector{<: Pair{M, <: Number} where M}) = map(x) do (a, b)
    a => uint8ize(b)
end
uint8ize(x::Vector{<: Pair{<: Number, M} where M}) = map(x) do (a, b)
    uint8ize(a) => b
end
uint8ize(x::Vector{<: Pair{<: Number, <: Number}}) = map(x) do (a, b)
    uint8ize(a) => uint8ize(b)
end
uint8ize(x::Vector) = map(uint8ize, x)

ishydrocarbon(::AbstractFunctionalGroup) = false
ishydrocarbon(::Methyl) = true
ishydrocarbon(::Ethyl) = true
ncarbon(::Methyl) = 1 
ncarbon(::Ethyl) = 2 

ispropertyempty(::Nothing) = true 
ispropertyempty(x::Vector) = isempty(x) 
ispropertyempty(x::Number) = iszero(x)
ispropertyempty(x) = false
ispropertynumber(::Vector{<: Pair{T, UInt8}}) where T = true 
ispropertynumber(::Number) = true 
ispropertynumber(x) = false 
ispropertypartialposition(x::Vector{<: Pair{UInt8}}) = any(x -> first(x) == 0, x)
ispropertypartialposition(x) = false 
ispropertyposition(::Vector{<: Pair{UInt8}}) = true 
ispropertyposition(x) = false 
function ispropertychirality(sub, ch)
    ispropertyposition(sub) || return false 
    isempty(sub) && return true 
    !isempty(filter(x -> last(x) != RSChirality(), ch)) # having sub with [R] or [S] or no [R] [S] nedded
end
function ispropertypartialchirality(sub, ch)
    ispropertypartialposition(sub) && return true
    chi = first.(filter(x -> last(x) != RSChirality(), ch)) # sub with [R] or [S] or no [R] [S] nedded
    for p in unique(first.(sub))
        !in(p, chi) && return true
    end
    false 
end