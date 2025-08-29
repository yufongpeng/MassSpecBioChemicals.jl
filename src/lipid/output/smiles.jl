
# MassSpecChemicals, SMILES
function getchemicalattr(lipid::FattyAcyl, ::Val{:SMILES}; onlycarbonchain = false)
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    smi = chemicalsmiles_carbonchain(lipid.chain)
    # reverse temporal solution
    lg, smib = match(r"^\(*([A-Z][a-z]*)(.*?)\)*$", chemicalsmiles(dehydrogengroup(lipid.backbone)))
    lg = lg == "H" ? "" : lg
    isempty(smib) ? string(lg, smi) : string(lg, "(", smib, ")", smi)
end

function getchemicalattr(lipid::FattyAcylEstolid, ::Val{:SMILES}; onlycarbonchain = false)
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    smi = chemicalsmiles_carbonchain(lipid.chain)
    smib = chemicalsmiles(lipid.backbone)
    cc = 0 
    i = 0
    while cc < lipid.position 
        i = nextind(smib, i)
        if smib[i] == 'C'
            cc += 1 
        end
    end
    string(smib[begin:i + 2], smi, smib[i + 3:end])
end

function getchemicalattr(lipid::WaxEster, ::Val{:SMILES}; onlycarbonchain = false)
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    smi = chemicalsmiles_carbonchain(lipid.chain)
    smib = chemicalsmiles(lipid.backbone)
    lipid.position <= 1 && return string("O(", smi, ")", smib[begin + 1:end])
    cc = 0 
    i = 0
    while cc < lipid.position 
        i = nextind(smib, i)
        if smib[i] == 'C'
            cc += 1 
        end
    end
    string(smib[begin:i + 2], smi, smib[i + 3:end])
end

function getchemicalattr(lipid::Monoradylglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    p = only(decode_sn(lipid))
    c = chemicalsmiles_carbonchain(lipid.chain)
    smi = chemicalsmiles(lipid.backbone; chain = (:C1, :C3))
    r = collect(eachmatch(r"\(O\)", smi))
    next = r[lastindex(r) - p + 1].match.offset
    string(smi[begin:next + 2], c, smi[next + 3:end])
end
function getchemicalattr(lipid::Diradylglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    position = decode_sn(lipid)
    p = sortperm(position; rev = true)
    position = position[p]
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain[p]]
    smi = chemicalsmiles(lipid.backbone; chain = (:C1, :C3))
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[[length(r) - t + 1 for t in position]], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

function getchemicalattr(lipid::Triradylglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    position = decode_sn(lipid)
    p = sortperm(position; rev = true)
    position = position[p]
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain[p]]
    smi = chemicalsmiles(lipid.backbone; chain = (:C1, :C3))
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[[length(r) - t + 1 for t in position]], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

function getchemicalattr(lipid::Radyldiglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    position = decode_sn(lipid)
    p = sortperm(position; rev = true, by = x -> (min(x, 3), x % 2)) # 3 4 2 1 
    position = map(position[p]) do x 
        findfirst(==(x), [3, 4, 2, 1])
    end 
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain[p]]
    smi = string(chemicalsmiles(first(getchaincomponent(lipid.backbone)); chain = (:C3, :O1))[begin:end - 1], chemicalsmiles(last(getchaincomponent(lipid.backbone)); chain = (:O1, :C3)))
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[position], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

function getchemicalattr(lipid::Radyltriglycerol, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    position = decode_sn(lipid)
    p = sortperm(position; rev = true, by = x -> (max(x, 2), x % 2)) # 4 3 1 2
    position = map(position[p]) do x 
        findfirst(==(x), [4, 3, 0, 1, 2])
    end 
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain[p]]
    smi = string(chemicalsmiles(first(getchaincomponent(lipid.backbone)); chain = (:C3, :O1))[begin:end - 1], chemicalsmiles(getchaincomponent(lipid.backbone)[begin + 1]; chain = (:O3, :O1))[begin:end - 1], chemicalsmiles(last(getchaincomponent(lipid.backbone)); chain = (:O1, :C3)))
    r = collect(eachmatch(r"\(O\)", smi))
    s = ""
    prev = firstindex(smi)
    for (m, c) in zip(r[position], cs)
        next = m.match.offset
        s *= smi[prev:next + 2]
        prev = next + 3
        s *= c
    end
    s *= smi[prev:end]
end

function getchemicalattr(lipid::SphingoidBase, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    chemicalsmiles_carbonchain(lipid.chain)
end

function getchemicalattr(lipid::Lysosulfonolipid, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    chemicalsmiles_carbonchain(lipid.chain)
end

function getchemicalattr(lipid::Ceramide, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain]
    cc = 0 
    i = 0
    bb = 0
    while cc < 2
        i = nextind(first(cs), i)
        if bb == 0 && first(cs)[i] == 'C'
            cc += 1 
        elseif first(cs)[i] == '('
            bb += 1
        elseif first(cs)[i] == ')'
            bb -= 1
        end
    end
    ids = findall("(N)", first(cs))
    id = ids[findfirst(>(i:i), ids)]
    string(first(cs)[begin:id[2]], cs[begin + 1], first(cs)[last(id):end])
end

function getchemicalattr(lipid::Sulfonolipid, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain]
    cc = 0 
    i = 0
    bb = 0
    while cc < 2
        i = nextind(first(cs), i)
        if bb == 0 && first(cs)[i] == 'C'
            cc += 1 
        elseif first(cs)[i] == '('
            bb += 1
        elseif first(cs)[i] == ')'
            bb -= 1
        end
    end
    ids = findall("(N)", first(cs))
    id = ids[findfirst(>(i:i), ids)]
    string(first(cs)[begin:id[2]], cs[begin + 1], first(cs)[last(id):end])
end

function getchemicalattr(lipid::Acylceramide, ::Val{:SMILES}; onlycarbonchain = false)
    # check
    lipid = onlycarbonchain ? dissociate_carbonchain(lipid) : lipid
    any(>=(fullstructurelevel), annotationlevel(lipid)) || return ""
    cs = [chemicalsmiles_carbonchain(c) for c in lipid.chain]
    cc = 0 
    i = 0
    bb = 0
    while cc < 2
        i = nextind(first(cs), i)
        if bb == 0 && first(cs)[i] == 'C'
            cc += 1 
        elseif first(cs)[i] == '('
            bb += 1
        elseif first(cs)[i] == ')'
            bb -= 1
        end
    end
    ids = findall("(N)", first(cs))
    id = ids[findfirst(>(i:i), ids)]
    string(first(cs)[begin:id[2]], cs[begin + 1], first(cs)[last(id):end])
end

getchemicalattr(lipid::Omodifiedradylglycerol, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_carbonchain(lipid); onlycarbonchain = true) : ""
getchemicalattr(lipid::Glycerophospholipid, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_carbonchain(lipid); onlycarbonchain = true) : ""
getchemicalattr(lipid::MixSphingoBone, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_carbonchain(lipid); onlycarbonchain = true) : ""
getchemicalattr(lipid::SphingoBone, ::Val{:SMILES}; onlycarbonchain = false) = onlycarbonchain ? chemicalsmiles(dissociate_carbonchain(lipid); onlycarbonchain = true) : ""

# to full structure level
# delete headgroup
# delete charged chain components
initial_smiles(::CarbonChain{Acyl}) = (["C(=O)"], [0])
initial_smiles(::CarbonChain{Alkyl}) = (["C"], [0])
initial_smiles(::CarbonChain{Alkenyl{ZConfiguration}}) = (["\\C", "=C", "/C"], [0, 0, 0])
initial_smiles(::CarbonChain{Alkenyl{EConfiguration}}) = (["\\C", "=C", "\\C"], [0, 0, 0])
initial_smiles(::CarbonChain{SPB}) = (["C", "C(N)"], [0, 0])
initial_smiles(::CarbonChain{SulfoSPB}) = (["S(=O)(=O)(O)C", "C(N)"], [0, 0])

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

getchemicalattr(chain::CarbonChain{<: CarbonChainType}, ::Val{:SMILES}; kwargs...) = chemicalsmiles_carbonchain(chain)

"""
    chemicalsmiles_carbonchain(chain)

SMILES representation of carbon chains
"""
function chemicalsmiles_carbonchain(chain::CarbonChain{<: CarbonChainType})
    i = j = k = 1
    pos = (isempty(chain.doublebond) || chain.doublebond == 0) ? [0x00] : sort!([db for db in chain.doublebond]; by = first)
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
            if last(pos[j]) == EZConfiguration() || last(pos[j]) == NoEZConfiguration()
                if isempty(quec)
                    push!(quec, "=X")
                    push!(posc, 0)
                else
                    quec[begin] = string("=", quec[begin])
                end
            elseif last(pos[j]) == ZConfiguration()
                if startswith(prec, '\\')
                    push_smiles_db!(quec, posc, "/")
                elseif startswith(prec, '/')
                    push_smiles_db!(quec, posc, "\\")
                else
                    prec = replace(prec, r"[XC]" => "\\C") # replace ? 
                    push_smiles_db!(quec, posc, "/")
                end
            elseif last(pos[j]) == EConfiguration()
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
            if smi == "(H)"
                nothing
            elseif !endswith(smi, ")")
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
