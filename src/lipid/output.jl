"""
    decode_sn(lipid)

Decode snpositioncode into vector of snposition
"""
function decode_sn(lipid::T) where {T <: Lipid}
    n = length(getlipidchain(lipid))
    position = zeros(Int, n)
    sn = lipid.sn
    ep = n - 1
    bs = nchainposition(T) + 1
    for i in eachindex(position)
        p, sn = divrem(sn, bs ^ ep)
        position[i] = p
        ep -= 1
    end
    position
end

"""
    decode_position(lipid)

Decode headgroup position code into vector of headgroup position
"""
function decode_position(lipid::T) where {T <: SphingoBone}
    p = Int(lipid.position)
    if hascycloheadgroup(lipid)
        p1, p2 = divrem(p, 32)
        iszero(p1) ? (p2, ) : (p1, p2)
    else
        (p, )
    end
end

decode_position(lipid::T) where {T <: SphingoBone{Nothing}} = ()
"""
    repr_db(chain)
    repr_db(dbs)

Decode double bond code into double bond representation
"""
function repr_db(c::CarbonChain) 
    p = repr_db_position(c)
    isempty(p) ? string(ndoublebond(c)) : string(ndoublebond(c), "(", p, ")")
end
"""
    repr_db_position(chain)
    repr_db_position(dbs)

Decode double bond code into double bond position representation
"""
repr_db_position(chain::CarbonChain) = repr_db_position(chain.doublebond)
function repr_db_position(dbs)
    join(decode_db_position(dbs), ",")
end
"""
    decode_db_position(db)

Decode double bond code into double bond position
"""
decode_db_position(chain::CarbonChain) = decode_db_position(chain.doublebond)
function decode_db_position(dbs::Vector)
    filter!(!isnothing, [_decode_db_position(x) for x in dbs])
end
decode_db_position(dbs) = String[]
function _decode_db_position(db::UInt8)
    a, b = divrem(db, 3)
    a == 0x00 && return nothing
    string(a, b == 0 ? "" : b == 1 ? "Z" : b == 2 ? "E" : throw(ArgumentError("Invalid fattyacyl chain")))
end

"""
    repr_sub(sub)

Convert chain substituent(s) into readble representation
"""
function repr_sub(sub::UInt8)
    sub == 0 ? "" : string(";O", sub > 1 ? Int(sub) : "")
end 
repr_sub(sub::Vector{<: Pair{OxygenAtom, UInt8}}) = repr_sub(only(sub))
function repr_sub(sub::Pair{OxygenAtom, UInt8})
    last(sub) == 0 ? "" : string(";O", last(sub) > 1 ? Int(last(sub)) : "")
end

function repr_sub(sub::Vector{<: Pair{F, UInt8} where F})
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ first)
    subs = String[""]
    ns = Int[0]
    for x in sub
        next, p = x
        next = sub_abbr(next)
        if next == last(subs)
            subs[end] += Int(p)
        else
            push!(subs, next)
            push!(ns, Int(p))
        end
    end
    popfirst!(subs)
    popfirst!(ns)
    s = ""
    for (sub, n) in zip(subs, ns)
        s *= n > 1 ? (length(sub) > 1 ? string(";(", sub, ")", n) : string(";", sub, n)) : endswith(sub, r"\d") ? string(";(", sub, ")") : string(";", sub)
    end
    s
end 

function repr_sub(sub::Vector{<: Pair{UInt8, F} where F})
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ last)
    s = ""
    prev = ""
    for x in sub
        p, next = x
        next = sub_abbr(next)
        if next == prev
            s *= string(",", Int(p), next)
        else
            s *= string(";", Int(p), next)
            prev = next
        end
    end
    s
end 
repr_sub(::Nothing) = ""

"""
    repr_singlechain(c::CarbonChain)

Representation of single carbon chain
"""
repr_singlechain(c::CarbonChain{<: SPB}) = string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Acyl}) = string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Alkyl}) = string("O-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Alkenyl}) = string("P-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
function repr_singlechain(c::CarbonChain{<: T}) where {T <: Tuple}
    if SPB in T.parameters
        return string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
    elseif Alkenyl in T.parameters # assume only P and A
        n = count(==(Alkenyl), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkenyl chain"))
        pre = string(pre, "P-")
    elseif Alkyl in T.parameters # assume only O and A
        n = count(==(Alkyl), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkyl chain"))
        pre = string(pre, "O-")
    else
        pre = ""
    end
    string(pre, ncarbon(c), ":", repr_db(c), repr_sub(c.substituent)) 
end

# to full structure level
# delete headgroup
# delete charged chain components
initial_smiles(::CarbonChain{Acyl}) = (["C(=O)"], [0])
initial_smiles(::CarbonChain{Alkyl}) = (["C"], [0])
initial_smiles(::CarbonChain{Alkenyl}) = (["\\C", "=C", "/C"], [0, 0, 0])
initial_smiles(::CarbonChain{SPB}) = (["C", "C(N)"], [0, 0])

function push_smiles_db!(quec, posc, slash)
    if length(quec) > 1
        quec[begin] = string("=", quec[begin])
        quec[begin + 1] = string(slash, quec[begin + 1])
    else
        if isempty(quec)
            push!(quec, "=X")
            push!(posc, 0)
        else
            quec[begin] = string("=", quec[begin])
        end
        push!(quec, string(slash, "X"))
        push!(posc, 0)
    end
end

"""
    chemicalsmiles_carbonchain(chain)

SMILES representation of carbon chains
"""
function chemicalsmiles_carbonchain(chain::CarbonChain{<: CarbonChainType})
    i = j = k = 1
    pos = iszero(chain.doublebond) ? [0x00] : sort!([divrem(db, 3) for db in chain.doublebond]; by = first)
    subs = isnothing(chain.substituent) ? Pair[] : sort(chain.substituent; by = first)
    s = ""
    cn = 0
    quec, posc = initial_smiles(chain)
    while i <= ncarbon(chain)
        if isempty(quec)
            prec = "C"
            numc = 0
        else
            prec = replace(popfirst!(quec), r"[XC]" => "C")
            numc = popfirst!(posc)
        end
        if numc > 0
            prec *= string(numc) 
        end
        if !startswith(prec, r"[\[\\/=#]*C") # -O-, -OO-
            s *= prec
            continue
        end
        if j <= lastindex(pos) && i == first(pos[j])
            if last(pos[j]) == 0
                if isempty(quec)
                    push!(quec, "=X")
                    push!(posc, 0)
                else
                    quec[begin] = string("=", quec[begin])
                end
            elseif last(pos[j]) == 1
                if startswith(prec, '\\')
                    push_smiles_db!(quec, posc, "/")
                elseif startswith(prec, '/')
                    push_smiles_db!(quec, posc, "\\")
                else
                    prec = replace(prec, r"[XC]" => "\\C") # replace ? 
                    push_smiles_db!(quec, posc, "/")
                end
            elseif last(pos[j]) == 2
                if startswith(prec, '\\')
                    push_smiles_db!(quec, posc, "\\")
                elseif startswith(prec, '/')
                    push_smiles_db!(quec, posc, "/")
                else
                    prec = replace(prec, r"[XC]" => "\\C")
                    push_smiles_db!(quec, posc, "\\")
                end
            end
            j += 1
        end
        while k <= lastindex(subs) && first(subs[k]) == i
            smi = chemicalsmiles(last(subs[k]))
            if !endswith(smi, ")")
                bsmi, msmi = match(r"^((?:\([^\)]*\))*)(.*)", smi)
                if isempty(bsmi) || !occursin(r"\d+", bsmi)
                    add = [e for (e, ) in eachmatch(r"(=*#*\[*[A-Z][a-z]*@*H*\]*)\d*", msmi)] # NO OTHER BRANCH
                    for i in eachindex(add)
                        if i > lastindex(posc)
                            push!(quec, add[i])
                            push!(posc, 0)
                        else
                            quec[i] = replace(quec[i], "X" => add[i])
                        end
                    end
                else
                    add = [string(E) => isempty(n) ? 0 : parse(Int, string(n)) for (E, n) in eachmatch(r"(=*#*\[*[A-Z][a-z]*@*H*\]*)(\d*)", msmi)] # NO OTHER BRANCH
                    new = last.(add)
                    org = new .- cn
                    for i in eachindex(new)
                        if i > lastindex(posc)
                            push!(quec, first(add[i]))
                            push!(posc, new[i])
                        elseif posc[i] > 0
                            if new[i] - posc[i] > 0
                                new[i:end] .-= 1
                                new[i] = posc[i]
                            end
                            quec[i] = replace(quec[i], "X" => first(add[i]))
                        else
                            posc[i] = new[i]
                            quec[i] = replace(quec[i], "X" => first(add[i]))
                        end
                    end
                    id = findall(>(0), org)
                    bsmi = replace(bsmi, (string.(org[id]) .=> string.(new[id]))...)
                    cn = maximum(posc)
                end
                prec = string(prec, bsmi)
            else
                prec = string(prec, smi)
            end
            k += 1
        end
        i += 1
        s *= prec
    end
    s
end

function Base.show(io::IO, level::AnnotationLevel)
    print(io, @match level begin
        &specieslevel                   => "specieslevel"
        &molecularspecieslevel          => "molecularspecieslevel"
        &phosphatepositionlevel         => "phosphatepositionlevel"
        &snpositionlevel                => "snpositionlevel"
        &dbpositionpartiallevel         => "dbpositionpartiallevel"
        &dbpositionlevel                => "dbpositionlevel"
        &dbconfigpartiallevel           => "dbconfigpartiallevel"
        &dbconfiglevel                  => "dbconfiglevel"
        &structuredefinedpartiallevel   => "structuredefinedpartiallevel"
        &structuredefinedlevel          => "structuredefinedlevel"
        &structurepositionpartiallevel  => "structurepositionpartiallevel"
        &structurepositionlevel         => "structurepositionlevel"
        &structureconfigpartiallevel    => "structureconfigpartiallevel"
        &structureconfiglevel           => "structureconfiglevel"
        &fullstructurelevel             => "fullstructurelevel"
        &completestructurelevel         => "completestructurelevel"
    end)
end

function Base.show(io::IO, level::PassLevel)
    print(io, @match level begin
        &passphosphatepositionlevel     => "passphosphatepositionlevel"
        &passsnpositionlevel            => "passsnpositionlevel"
    end)
end