struct SaccharideParser{Bool} <: AbstractChemicalParser end
SaccharideParser(; explicit_config = false) = SaccharideParser{explicit_config}()
parse_chemical(::SaccharideParser{true}, x; kwargs...) = parse_glycan(x; explicit_config = true)
parse_chemical(::SaccharideParser{false}, x; kwargs...) = parse_glycan(x)

# specify num -> num
# default 6 -> 1 === dehydrogen -> dehydroxy 
# other than 6, => reverse smiles of sub
function chemicalsmiles(sugar::Hexose; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}C{5}(O1)C{4}(O)C{3}(O)C{2}(O)C{1}1O"
end
function chemicalsmiles(sugar::Glucose; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end
function chemicalsmiles(sugar::Galactose; kwargs...)
    isnothing(sugar.substituent) && return "OC{6}[C{5}@@H](O1)[C{4}@H](O)[C{3}@H](O)[C{2}@@H](O)[C{1}H]1O"
end

function _substituentelements(chain, x)
    ec = [losswaterelements(m, i, i % length(chain)) for (i, m) in enumerate(chain)]
    es = last(ec)
    for (e, n) in leavinggroupelements(leavinggroup(x))
        if n > 0
            i = findfirst(y -> first(y) == e, es)
            es[i] = e => (last(es[i]) + n)
        else
            while n < 0
                i = findfirst(y -> first(y) == e, es)
                es[i] = e => (last(es[i]) - 1)
                n += 1
                filter!(y -> last(y) > 0, es)
            end
        end
    end
    vcat(ec...)
end
chemicalelements(x::Substituent{<: Glycan}; kwargs...) = _substituentelements(getchaincomponent(parentchemical(x)), x)
function chemicalelements(x::Substituent{<: GlyComp}; kwargs...)
    chain = mapreduce(vcat, getchaincomponent(parentchemical(x))) do (m, n)
        fill(m, n)
    end
    _substituentelements(chain, x)
end
