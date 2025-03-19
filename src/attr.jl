parentchemical(fg::Substituent) = fg.chemical
parentchemical(sm::SubstitutedChemical) = sm.chemical
parentchemical(::FunctionalGroup{M}) where M = M()
parentchemical(x::XLinkedFunctionalGroup) = parentchemical(x.functionalgroup)
parentchemical(::M) where {M <: UnknownGroup} = UnknownChemical{M}()
leavinggroup(::FunctionalGroup{M, S}) where {M, S} = S()

conjugation(m::AbstractChemical) = Substituent(Ndehydrogen, m)

includeSIL(T) = Union{<: T, <: IsotopiclabeledChemical{<: T}}

leavinggroupelements(::Dehydrogen) = ["H" => -1]
leavinggroupelements(::Odehydrogen) = ["H" => -1]
leavinggroupelements(::Ndehydrogen) = ["H" => -1]
leavinggroupelements(::Didehydrogen) = ["H" => -2]
leavinggroupelements(::Dehydroxy) = ["O" => -1, "H" => -1]
leavinggroupelements(::Deamine) = ["N" => -1, "H" => -2]
leavinggroupelements(::Demethyl) = ["C" => -1, "H" => -3]

dehydrogenposition(a::AbstractChemical) = 0x00
dehydrogenposition(a::DehydratedChemical) = dehydrogenposition(first(getchaincomponent(a)))
dehydroxyposition(a::AbstractChemical) = 0x00
dehydroxyposition(a::DehydratedChemical) = dehydroxyposition(last(getchaincomponent(a)))

isnulllinkage(::Type{DehydratedChemical}, l) = l == lk(0x00)
isnulllinkage(::Type{DehydratedChemical}, l::Pair) = all(==(lk(0x00)), l)
isnulllinkage(::Type{ChainedChemical}, l) = first(l) == lk(0x00)
isnulllinkage(::Type{ChainedChemical}, l::Pair) = all(x -> ==(first(x), lk(0x00)), l)

makelinkage(::Type{ChainedChemical}, a, b) = (lk(dehydroxyposition(a)), Dehydroxy()) => (lk(dehydrogenposition(b)), Dehydrogen())
makelinkage(::Type{DehydratedChemical}, a, b) = lk(dehydroxyposition(a)) => lk(dehydrogenposition(b))

transformlinkage(::Type{<: AbstractChemical}, m::AbstractChemical) = getchainlinkage(m)
function transformlinkage(::Type{DehydratedChemical}, m::ChainedChemical)
    Tuple(first(first(ls)) => first(last(ls)) for ls in getchainlinkage(m))
end
function transformlinkage(::Type{ChainedChemical}, m::DehydratedChemical)
    Tuple((first(ls), Dehydroxy()) => (last(ls), Dehydrogen()) for ls in getchainlinkage(m))
end

getchaincomponent(m::AbstractChemical) = (m, )
getchaincomponent(m::ChainedChemical) = m.chain
getchaincomponent(m::DehydratedChemical) = m.chain
getchainlinkage(m::AbstractChemical) = ()
getchainlinkage(m::ChainedChemical) = m.linkage
getchainlinkage(m::DehydratedChemical) = m.linkage
ischainedchemical(m::AbstractChemical) = false
ischainedchemical(m::ChainedChemical) = true
ischainedchemical(m::DehydratedChemical) = true