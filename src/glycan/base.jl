length(::Monosaccharide) = 1
length(sugar::AbstractGlycan) = sum(length, getchaincomponent(sugar))
length(sugar::GlyComp) = sum(last, getchaincomponent(sugar))

==(::T, ::T) where {T<:GlycanType} = true
==(::T, ::S) where {T<:GlycanType, S<:GlycanType} = false
==(::T, ::T) where {T<:MonosaccharideForm} = true
==(::T, ::S) where {T<:MonosaccharideForm, S<:MonosaccharideForm} = false
hash(::T, h::UInt) where {T<:GlycanType} = hash(T, h)
function hash(x::SM1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
function hash(x::GM1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
function hash(x::GD1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
function hash(x::GT1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
function hash(x::GQ1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
function hash(x::GP1, h::UInt) 
    for i in x.isomer
        h = hash(i, h)
    end
    h
end
hash(::T, h::UInt) where {T<:MonosaccharideForm} = hash(T, h)
copy(::T) where {T<:GlycanType} = T()
copy(x::SM1) = SM1(map(copy, x.isomer))
copy(x::GM1) = GM1(map(copy, x.isomer))
copy(x::GD1) = GD1(map(copy, x.isomer))
copy(x::GT1) = GT1(map(copy, x.isomer))
copy(x::GQ1) = GQ1(map(copy, x.isomer))
copy(x::GP1) = GP1(map(copy, x.isomer))
copy(::T) where {T<:MonosaccharideForm} = T()
