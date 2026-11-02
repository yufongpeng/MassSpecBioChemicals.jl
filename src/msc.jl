function pushchemicalattr!(m::T, attr::Symbol, value::IsotopicLabel) where {T <: AbstractChemical}
    if hasfield(m, attr)
        setfield!(getfield(m, attr), :isotopiclabel, value)
    else
        attribute = getchemicalattr(m)
        isnothing(attribute) && return pushchemicalattr_new!(m, attr, value)
        if isnothing(attribute.attribute)
            attribute.attribute = Pair{Symbol, Any}[attr => value]
        else
            id = findfirst(==(attr) ∘ first, attribute.attribute)
            if isnothing(id)
                push!(attribute.attribute, attr => value)
            else
                attribute.attribute[id] = attr => value
            end
        end
    end
end
