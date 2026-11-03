"""
    AbstractSaccharideParser{E}

Parsers of saccharides. The type parameter `E` is `explicit_config`: if `true`, a monosaccharide without `D-`/`L-` prefix is parsed as `DLForm`.
"""
abstract type AbstractSaccharideParser{E} <: AbstractChemicalParser end
"Parser for monosaccharide, e.g. `Glc`."
struct MonosaccharideParser{E} <: AbstractSaccharideParser{E} end
"Parser for glycan (monosaccharides and glycans connected by linkages), e.g. `Galβ-4Glc`."
struct SaccharideParser{E} <: AbstractSaccharideParser{E} end
"Parser for named series glycan, e.g. `GQ1bα(3[3]NeuGc)`."
struct SeriesGlycanParser{E} <: AbstractSaccharideParser{E} end
"Parser for glycan composition, e.g. `HexNAcHex3`."
struct GlycompParser{E} <: AbstractSaccharideParser{E} end
MonosaccharideParser(; explicit_config = false) = MonosaccharideParser{explicit_config}()
SaccharideParser(; explicit_config = false) = SaccharideParser{explicit_config}()
SeriesGlycanParser(; explicit_config = false) = SeriesGlycanParser{explicit_config}()
GlycompParser(; explicit_config = false) = GlycompParser{explicit_config}()
parse_chemical(p::MonosaccharideParser, x; kwargs...) = parse_monosaccharide(p, x)
parse_chemical(p::SaccharideParser, x; kwargs...) = parse_glycan(p, x)
parse_chemical(p::SeriesGlycanParser, x; kwargs...) = parse_series_glycan(p, x)
parse_chemical(p::GlycompParser, x; kwargs...) = parse_glycomp(p, x)

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
