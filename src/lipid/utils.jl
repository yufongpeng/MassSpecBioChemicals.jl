uint8ize(x::Number) = UInt8(x)
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