include("abbr.jl")
getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: LipidChemical} = chemicalname(lipid.chemical; kwargs...)
# Lipid
# MassSpecChemicals
function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Glycerolipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(chiral_class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(chiral_class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(chiral_class_abbr(lipid), " ", join(rp, "/"))
end
getchemicalattr(lipid::GlycerophosphoNacylethanolamine, ::Val{:name}; kwargs...) = string(chiral_class_abbr(lipid), " ", repr_singlechain(lipid.chain))
function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Glycerophospholipid}
    position = decode_sn(lipid)
    if any(==(0), position) && any(>(0), position)
        pos = chainposition(lipid)
        return string(chiral_class_abbr(lipid), " ", 
                        join([string(repr_singlechain(c), p > 0 ? string("(", pos[p], ")") : "") for (p, c) in zip(position, getlipidchain(lipid))], "_"))
    elseif any(==(0), position)
        return string(chiral_class_abbr(lipid), " ", join(repr_singlechain.(getlipidchain(lipid)), "_"))
    end
    maxsn = nchainposition(T)
    rp = String[]
    cs = repr_singlechain.(getlipidchain(lipid))
    for i in 1:maxsn
        j = findfirst(==(i), position)
        push!(rp, isnothing(j) ? "0:0" : cs[j])
    end
    string(chiral_class_abbr(lipid), " ", join(rp, "/"))
end

# function getchemicalattr(lipid::T) where {T <: FattyAcyl{B, <: Tuple}} where B
#     string(chiral_class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
# end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: FattyAcyl}
    string(chiral_class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Union{NacylAlkylAmine, FattyAcylEster}}
    if ncarbon(lipid.backbone.chain) == 0
        string(chiral_class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    else
        string(replace(chemicalname(lipid.backbone), chiral_class_abbr(lipid.backbone) => chiral_class_abbr(lipid)), "/", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    end
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: FattyAcylEstolid}
    if ncarbon(lipid.backbone.chain) == 0
        string(chiral_class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""))
    elseif iszero(lipid.position)
        string(chiral_class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", replace(chemicalname(lipid.backbone), string(chiral_class_abbr(lipid.backbone), " ") => ""))
    else
        i = Int(lipid.position)
        l = lipid.backbone
        id = findfirst(x -> (Int(first(x)) == i) && (last(x) == Hydroxy()), l.chain.substituent)
        sub = l.chain.substituent[setdiff(eachindex(l.chain.substituent), [id])]
        id = findfirst(x -> Int(first(x)) == i, l.chain.chirality)
        rs = last(l.chain.chirality[id])
        l = MonoFattyAcyl(l.backbone, CarbonChain{Acyl}(l.chain.carbon, l.chain.doublebond, sub, l.chain.chirality, l.chain.isotopiclabel))
        string(chiral_class_abbr(lipid), " ", replace(repr_singlechain(lipid.chain), r"[OP]-" => ""), "/", i, "O(", l, ")", repr_chirality(rs))
    end
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Sphingolipid{H, <: Tuple}} where H
    isnothing(lipid.headgroup) && return string(chiral_class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), "/"))
    position = repr_position(lipid)
    string(chiral_class_abbr(lipid), position, " ", join(repr_singlechain.(lipid.chain), "/"))
end

function getchemicalattr(lipid::T, ::Val{:name}; kwargs...) where {T <: Sphingolipid{H, <: CarbonChain}} where H
    isnothing(lipid.headgroup) && return string(chiral_class_abbr(lipid), " ", repr_singlechain(lipid.chain))
    position = repr_position(lipid)
    string(chiral_class_abbr(lipid), position, " ", repr_singlechain(lipid.chain))
end

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
    repr_position(lipid)

Decode headgroup position code into headgroup position representation
"""
function repr_position(lipid::T) where {T <: SphingoBone}
    p = Int(lipid.position)
    if length(tuplize(lipid.chirality)) == 2
        p1, p2 = divrem(p, 32)
        if iszero(p1) && iszero(p2)
            ""
        elseif iszero(p1) 
            rc = repr_chirality(lipid.chirality; bracket = false)
            string("(", p2, rc, ")") 
        else
            rc = repr_chirality.(lipid.chirality; bracket = false)
            string("(", p1, rc[1], ",", p2, rc[2], ")")
        end
    elseif iszero(p)
        ""
    else
        rc = repr_chirality(lipid.chirality; bracket = false)
        string("(", p, rc, ")") 
    end
end

repr_position(lipid::T) where {T <: SphingoBone{Nothing}} = ""

function repr_position(lipid::T) where {T <: MixSphingoBone}
    ps = Int.(lipid.position)
    s = Vector{String}(undef, length(ps))
    for (i, (p, c)) in enumerate(zip(ps, lipid.chirality))
        if length(tuplize(c)) == 2
            p1, p2 = divrem(p, 32)
            s[i] = if iszero(p1) && iszero(p2)
                ""
            elseif iszero(p1) 
                rc = repr_chirality(c; bracket = false)
                string(p2, rc) 
            else
                rc = repr_chirality.(c; bracket = false)
                string(p1, rc[1], ",", p2, rc[2])
            end
        elseif iszero(p)
            s[i] = ""
        else
            rc = repr_chirality(c; bracket = false)
            s[i] = string(p, rc) 
        end
    end
    if all(isempty, s)
        "" 
    else
        string("(", join(s, ","), ")")
    end
end

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
function _decode_db_position(db::Pair)
    a, b = db
    a == 0x00 && return nothing
    string(a, b == NoEZConfiguration() ? "" : b == EZConfiguration() ? "" : b == ZConfiguration() ? "Z" : b == EConfiguration() ? "E" : throw(ArgumentError("Invalid E/Z configuration.")))
end

"""
    repr_sub(sub, chirality)

Convert chain substituent(s) into readble representation
"""
function repr_sub(sub::UInt8, chirality)
    sub == 0 ? "" : string(";O", sub > 1 ? Int(sub) : "")
end 
repr_sub(sub::Vector{<: Pair{OxygenAtom, UInt8}}, chirality) = repr_sub(only(sub), chirality)
function repr_sub(sub::Pair{OxygenAtom, UInt8}, chirality)
    last(sub) == 0 ? "" : string(";O", last(sub) > 1 ? Int(last(sub)) : "")
end

function repr_sub(sub::Vector{<: Pair{F, UInt8} where F}, chirality)
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

function repr_sub(sub::Vector{<: Pair{UInt8, F} where F}, chirality)
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ last)
        del = Int[]
    for (i, v) in enumerate(sub)
        if last(v) isa Hydrogen 
            if (first(v) => RSChirality()) in chirality || (first(v) => AChirality()) in chirality
                push!(del, i)
            end
        end
    end
    deleteat!(sub, del)
    ch = Dict(chirality)
    s = ""
    prev = ""
    for x in sub
        p, next = x
        next = sub_abbr(next)
        c = get(ch, p, AChirality())
        rc = repr_chirality(c)
        if next == prev
            s *= string(",", Int(p), next, rc)
        else
            s *= string(";", Int(p), next, rc)
            prev = next
        end
    end
    s
end 
repr_sub(::Nothing, chirality) = ""

repr_chirality(::AChirality; bracket = true) = ""
repr_chirality(::RSChirality; bracket = true) = ""
repr_chirality(::RChirality; bracket = true) = bracket ? "[R]" : "R"
repr_chirality(::SChirality; bracket = true) = bracket ? "[S]" : "S"
repr_chirality(::NoDLForm; bracket = true) = ""
repr_chirality(::DLForm; bracket = true) = "" # Other than Glycan
repr_chirality(::DForm; bracket = true) = "D"
repr_chirality(::LForm; bracket = true) = "L"

# function repr_chainchirality(chirality) 
#     i = findfirst(x -> ==(first(x), 0xff), chirality)
#     isnothing(i) ? "" : repr_chirality(last(chirality[i]))
# end
# repr_chainchirality(chirality::Nothing) = ""

"""
    repr_singlechain(c::CarbonChain)

Representation of single carbon chain
"""
repr_singlechain(c::CarbonChain{<: AbstractSPB}) = string(ncarbon(c), ":", repr_db(c), repr_sub(skip_olinkage(c.substituent), c.chirality)) # skip XLinkedFunctionalGroup{OLinkage}
repr_singlechain(c::CarbonChain{<: Acyl}) = string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
repr_singlechain(c::CarbonChain{<: Alkyl}) = string("O-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
repr_singlechain(c::CarbonChain{<: Alkenyl{ZConfiguration}}) = string("P-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
repr_singlechain(c::CarbonChain{<: Alkenyl{EConfiguration}}) = string("E-", ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
function repr_singlechain(c::CarbonChain{<: T}) where {T <: Tuple}
    if SPB in T.parameters
        return string(ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
    elseif Alkenyl{ZConfiguration} in T.parameters # assume only P and A
        n = count(==(Alkenyl{ZConfiguration}), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkenyl chain"))
        pre = string(pre, "P-")
    elseif Alkenyl{EConfiguration} in T.parameters # assume only P and A
        n = count(==(Alkenyl{EConfiguration}), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkenyl chain"))
        pre = string(pre, "E-")
    elseif Alkyl in T.parameters # assume only O and A
        n = count(==(Alkyl), T.parameters)
        pre = n == 1 ? "" : n == 2 ? "d" : n == 3 ? "t" : n == 4 ? "e" : throw(ArgumentError("Too many alkyl chain"))
        pre = string(pre, "O-")
    else
        pre = ""
    end
    string(pre, ncarbon(c), ":", repr_db(c), repr_sub(c.substituent, c.chirality)) 
end

function skip_olinkage(sub::Vector{<: Pair{<: AbstractFunctionalGroup, UInt8}})
    del = Int[]
    for (i, v) in enumerate(sub)
        first(v) isa XLinkedFunctionalGroup{OLinkage} && push!(del, i)
    end
    deleteat!(deepcopy(sub), del)
end

function skip_olinkage(sub::Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}})
    del = Int[]
    for (i, v) in enumerate(sub)
        last(v) isa XLinkedFunctionalGroup{OLinkage} && push!(del, i)
    end
    deleteat!(deepcopy(sub), del)
end
skip_olinkage(sub) = sub

repr_class_rs(lipid::Lipid) = ""
function repr_class_rs(lipid::Union{Glycerolipid, Glycerophospholipid})
    ch = map(x -> repr_chirality(x; bracket = false), tuplize(lipid.chirality))
    dot = ["'" ^ (i - 1) for i in eachindex(ch)]
    s = "("
    for (c, d) in zip(ch, dot)
        isempty(c) && continue
        s = string(s, c, d, ",")
    end
    s = s[begin:end - 1]
    isempty(s) ? "" : string(s, ")-")
end
for x in NAAA
    T = PROTEIN_3LETTER_AA[string(x)]
    @eval function repr_class_rs(c::NacylAmine{<: includeSIL($T)}) 
        rs = repr_chirality(first(typeof(first(getlipidbody(c))).parameters)())
        isempty(rs) ? "" : string(rs, "-")
    end
end

function repr_class_rs(c::FattyAcylCarnitine) 
    rs = repr_chirality(first(typeof(first(getlipidbody(c))).parameters)())
    isempty(rs) ? "" : string(rs, "-")
end

include("smiles.jl")

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