repr_linkage(ls::Pair{<: AbstractLinkageposition, <: AbstractLinkageposition}) = string(repr_linkage(first(ls)), "-", repr_linkage(last(ls)))
repr_linkage(l::Linkageposition) = isnothing(l.position) ? "" : l.position > 0 ? string(Int(l.position)) : ""
repr_linkage(::Nothing) = ""

repr_sil(::Type, x::IsotopicLabel) = repr_sil(x)
repr_sil(::IsotopicLabel) = "sil"

repr_chirality(::Type{T}) where {T<:AbstractChirality} = repr_chirality(T())
repr_chirality(::AbstractChirality) = "" 
repr_chirality(::DForm) = "D-"
repr_chirality(::LForm) = "L-"
repr_chirality(::RChirality) = "(R)-"
repr_chirality(::SChirality) = "(S)-"
repr_chirality(::EConfiguration) = "(E)-"
repr_chirality(::ZConfiguration) = "(Z)-"
repr_chirality(::ClockwiseRotation) = "(+)-"
repr_chirality(::CounterclockwiseRotation) = "(-)-"