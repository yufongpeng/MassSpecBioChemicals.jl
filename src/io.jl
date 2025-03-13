repr_linkage(ls::Pair{<: AbstractLinkageposition, <: AbstractLinkageposition}) = string(repr_linkage(first(ls)), "-", repr_linkage(last(ls)))
repr_linkage(l::Linkageposition) = isnothing(l.position) ? "" : l.position > 0 ? string(Int(l.position)) : ""
repr_linkage(::Nothing) = ""