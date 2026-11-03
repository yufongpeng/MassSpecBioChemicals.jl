length(cc::ChainedChemical) = length(getchaincomponent(cc))
length(cc::DehydratedChemical) = length(getchaincomponent(cc))

# Types other than `AttributedChemical`: singletons (no fields) and wrappers/composites of other instances.
const PlainTypes = Union{LeavingGroup, AbstractConfiguration, AbstractLinkageposition, AbstractFunctionalGroup, UnknownChemical, IsotopicLabel}

function ==(x::T, y::T) where {T <: PlainTypes}
    for f in fieldnames(T)
        a = getfield(x, f)
        b = getfield(y, f)
        eq = f === :elements ? chemicalelements_equal(a, b) : a == b
        eq || return false
    end
    true
end

function hash(x::T, h::UInt) where {T <: PlainTypes}
    h = hash(T, h)
    for f in fieldnames(T)
        v = getfield(x, f)
        h = f === :elements ? hash(sum(hash, sort_unique_elements(v); init = zero(h)), h) : hash(v, h)
    end
    h
end

copy(x::T) where {T <: PlainTypes} = fieldcount(T) == 0 ? x : T((safecopy(getfield(x, f)) for f in fieldnames(T))...)
