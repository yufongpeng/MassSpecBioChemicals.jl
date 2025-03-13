chainposition(::Type{<: Hydrocarbon}) = ["hydrocarbon"]
chainposition(::Type{<: FattyAcid}) = ["acid"]
chainposition(::Type{<: FattyAlcohol}) = ["alcohol"]
chainposition(::Type{<: FattyAldehyde}) = ["acyl"]
chainposition(::Type{<: FattyAmide}) = ["nacyl"]
chainposition(::Type{<: FattyAmine}) = ["amine"]
chainposition(::Type{<: FattyAcylCarnitine}) = ["oacyl"]
chainposition(::Type{<: FattyAcylCoA}) = ["sacyl"]
chainposition(::Type{<: NacylAmine}) = ["amine", "nacyl"]
chainposition(::Type{<: FattyAcylEster}) = ["alcohol", "oacyl"]
chainposition(::Type{<: WaxEster}) = ["alcohol", "oacyl"]
chainposition(::Type{<: Glycerophospholipid}) = ["sn-1", "sn-2"]
chainposition(::Type{<: Sphingolipid}) = ["lcb", "nacyl"]
chainposition(::Type{<: Sphingolipid{H, <: CarbonChain{SPB}}}) where H = ["lcb"]
chainposition(::Type{<: Glycerolipid}) = ["sn-1", "sn-2", "sn-3"]
chainposition(::Type{<: Omodifiedradylglycerol}) = ["sn-1", "sn-2"]
chainposition(::Type{<: Bisradylglycerophosphoglycerol}) = ["sn-1", "sn-2", "sn-1'", "sn-2'"]
chainposition(::Type{<: Bisradylglycerophosphate}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
chainposition(lipid::T) where {T <: Lipid} = chainposition(T)
nchainposition(::Type{T}) where {T <: Lipid} = length(chainposition(T))
nchainposition(lipid::T) where {T <: Lipid} = length(chainposition(T))
ncarbonchain(lipid::Lipid) = ncarbonchain(lipid.chain)
ncarbonchain(chains::Tuple) = sum(ncarbonchain, chains)
ncarbonchain(::CarbonChain{<: CarbonChainType}) = 1
ncarbonchain(::CarbonChain{<: T}) where {T <: Tuple} = length(T.parameters)
ncarbon(chains::Tuple) = sum(ncarbon, chains)
ncarbon(chain::CarbonChain) = chain.carbon
ndoublebond(chain::CarbonChain{S, UInt8}) where S = Int(chain.doublebond)
ndoublebond(chain::CarbonChain{S}) where S = length(chain.doublebond)

hascycloheadgroup(::Lipid) = false
hascycloheadgroup(::Sphingolipid) = false
hascycloheadgroup(::CeramidePhosphate) = true
hascycloheadgroup(::SphingoidBasePhosphate) = true

getlipidchain(lipid::Lipid) = tuplize(lipid.chain)
getlipidbody(lipid::Lipid) = tuplize(lipid.backbone)
getlipidbody(lipid::Sphingolipid) = tuplize(lipid.headgroup)

dehydroxyposition(::FattyAcid) = nothing
dehydrogenposition(::FattyAcid) = missing
dehydroxyposition(::FattyAlcohol) = nothing
dehydrogenposition(::FattyAlcohol) = nothing
dehydroxyposition(::FattyAmine) = nothing
dehydrogenposition(::FattyAmine) = nothing
# 5 isspecieslevel
# 3/4 isphosphatepositionlevel
# 4 ismolecularspecieslevel
# 3 issnpositionlevel
# 3 isdbpositionlevel
# 3 isstructuredefinedlevel
# 2 isfullstructurelevel
# 1 iscompletestructurelevel

# levels are positive listed, not listed => not in the level

# @enum Level begin
#     completestructurelevel
#     structureconfiglevel
#     fullstructurelevel
#     structureconfigpartiallevel
#     structurepositionlevel
#     structurepositionpartiallevel
#     structuredefinedlevel
#     structuredefinedpartiallevel
#     dbconfiglevel
#     dbconfigpartiallevel
#     dbpositionlevel
#     dbpositionpartiallevel
#     snpositionlevel
#     phosphatepositionlevel
#     molecularspecieslevel
#     specieslevel
# end

abstract type LipidAnnotationLevel end
struct AnnotationLevel <: LipidAnnotationLevel
    level::UInt8
end
struct PassLevel <: LipidAnnotationLevel
    level::UInt8
end
const completestructurelevel = AnnotationLevel(0x00)
const structureconfiglevel = AnnotationLevel(0x01)
const fullstructurelevel = AnnotationLevel(0x02)
const structureconfigpartiallevel = AnnotationLevel(0x03)
const structurepositionlevel = AnnotationLevel(0x04)
const structurepositionpartiallevel = AnnotationLevel(0x05)
const structuredefinedlevel = AnnotationLevel(0x06)
const structuredefinedpartiallevel = AnnotationLevel(0x07)
const dbconfiglevel = AnnotationLevel(0x08)
const dbconfigpartiallevel = AnnotationLevel(0x09)
const dbpositionlevel = AnnotationLevel(0x0a)
const dbpositionpartiallevel = AnnotationLevel(0x0b)
const snpositionlevel = AnnotationLevel(0x0c)
const phosphatepositionlevel = AnnotationLevel(0x0d)
const molecularspecieslevel = AnnotationLevel(0x0e)
const specieslevel = AnnotationLevel(0x0f)
const passsnpositionlevel = PassLevel(0x0c)
const passphosphatepositionlevel = PassLevel(0x0d)

function levelpriority(level::LipidAnnotationLevel)
    @match level begin
        &specieslevel                   => (0x00, 0x00, 0x00, 0x00, 0x00, 0x00)
        &molecularspecieslevel          => (0x01, 0x01, 0x01, 0x01, 0x01, 0x00)
        &phosphatepositionlevel         => (0x00, 0x00, 0x00, 0x00, 0x00, 0x01)
        &passphosphatepositionlevel     => (0x00, 0x00, 0x00, 0x00, 0x00, 0x02)
        &snpositionlevel                => (0x01, 0x01, 0x01, 0x01, 0x02, 0x00)
        &passsnpositionlevel            => (0x01, 0x01, 0x01, 0x01, 0x03, 0x00)
        &dbpositionpartiallevel         => (0x01, 0x01, 0x02, 0x01, 0x01, 0x00)
        &dbpositionlevel                => (0x01, 0x01, 0x02, 0x02, 0x01, 0x00)
        &dbconfigpartiallevel           => (0x01, 0x01, 0x03, 0x01, 0x01, 0x00)
        &dbconfiglevel                  => (0x01, 0x01, 0x03, 0x03, 0x01, 0x00)
        &structuredefinedpartiallevel   => (0x02, 0x01, 0x01, 0x01, 0x01, 0x00) # headgroup comp/chain mod comp
        &structuredefinedlevel          => (0x02, 0x02, 0x01, 0x01, 0x01, 0x00) # all comp
        &structurepositionpartiallevel  => (0x03, 0x01, 0x01, 0x01, 0x01, 0x02) # headgroup position/chain mod position
        &structurepositionlevel         => (0x03, 0x03, 0x01, 0x01, 0x01, 0x02) # all position
        &structureconfigpartiallevel    => (0x04, 0x01, 0x01, 0x01, 0x01, 0x02)
        &structureconfiglevel           => (0x04, 0x04, 0x01, 0x01, 0x01, 0x02)
        &fullstructurelevel             => (0x03, 0x03, 0x03, 0x03, 0x03, 0x02)
        &completestructurelevel         => (0x04, 0x04, 0x03, 0x03, 0x03, 0x02)
    end
end

function isless(x::LipidAnnotationLevel, y::LipidAnnotationLevel)
    result = false
    for (i, j) in zip(levelpriority(x), levelpriority(y))
        isless(j, i) && return false
        result = result || isless(i, j) 
    end
    result
end

function ispartial(level::LipidAnnotationLevel)
    level in (dbpositionpartiallevel, dbconfigpartiallevel, structuredefinedpartiallevel, structurepositionpartiallevel, structureconfigpartiallevel)
end
function isaddtional(level::LipidAnnotationLevel)
    level in (dbconfigpartiallevel, dbconfiglevel, structurepositionpartiallevel, structurepositionlevel, structureconfigpartiallevel, structureconfiglevel)
end

function transform_additional(level::LipidAnnotationLevel)
    @match level begin
        &specieslevel                   => specieslevel
        &molecularspecieslevel          => molecularspecieslevel
        &phosphatepositionlevel         => phosphatepositionlevel
        &passphosphatepositionlevel     => passphosphatepositionlevel
        &snpositionlevel                => snpositionlevel
        &passsnpositionlevel            => passsnpositionlevel
        &dbpositionpartiallevel         => dbpositionpartiallevel
        &dbpositionlevel                => dbpositionlevel
        &dbconfigpartiallevel           => dbpositionpartiallevel
        &dbconfiglevel                  => dbpositionlevel
        &structuredefinedpartiallevel   => structuredefinedpartiallevel # headgroup comp/chain mod comp
        &structuredefinedlevel          => structuredefinedlevel # all comp
        &structurepositionpartiallevel  => structuredefinedpartiallevel # headgroup position/chain mod position
        &structurepositionlevel         => structuredefinedlevel # all 
        &structureconfigpartiallevel    => structuredefinedpartiallevel
        &structureconfiglevel           => structuredefinedlevel
        &fullstructurelevel             => fullstructurelevel
        &completestructurelevel         => completestructurelevel
    end
end
function transform_partial(level::LipidAnnotationLevel)
    @match level begin
        &specieslevel                   => specieslevel
        &molecularspecieslevel          => molecularspecieslevel
        &phosphatepositionlevel         => phosphatepositionlevel
        &passphosphatepositionlevel     => passphosphatepositionlevel
        &snpositionlevel                => snpositionlevel
        &passsnpositionlevel            => passsnpositionlevel
        &dbpositionpartiallevel         => molecularspecieslevel
        &dbpositionlevel                => dbpositionlevel
        &dbconfigpartiallevel           => molecularspecieslevel
        &dbconfiglevel                  => dbpositionlevel
        &structuredefinedpartiallevel   => molecularspecieslevel # headgroup comp/chain mod comp
        &structuredefinedlevel          => structuredefinedlevel # all comp
        &structurepositionpartiallevel  => molecularspecieslevel # headgroup position/chain mod position
        &structurepositionlevel         => structurepositionlevel # all position
        &structureconfigpartiallevel    => molecularspecieslevel
        &structureconfiglevel           => structureconfiglevel
        &fullstructurelevel             => fullstructurelevel
        &completestructurelevel         => completestructurelevel
    end
end

const ANNOTATIONLEVELGAPCHECK = [
    fullstructurelevel,
    structureconfiglevel,
    structureconfigpartiallevel,
    structurepositionlevel,
    structurepositionpartiallevel,
    structuredefinedlevel,
    structuredefinedpartiallevel,
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel,
    snpositionlevel,
    phosphatepositionlevel,
    molecularspecieslevel,
    specieslevel
]

const ANNOTATIONLEVELMASTER = Dict(
    specieslevel                    => [molecularspecieslevel, phosphatepositionlevel],
    molecularspecieslevel           => [structuredefinedpartiallevel, dbpositionpartiallevel, snpositionlevel],
    phosphatepositionlevel          => [passphosphatepositionlevel],
    passphosphatepositionlevel      => [structurepositionpartiallevel], # Not check
    snpositionlevel                 => [passsnpositionlevel],
    passsnpositionlevel             => [fullstructurelevel], # Not check
    dbpositionpartiallevel          => [dbpositionlevel, dbconfigpartiallevel],
    dbpositionlevel                 => [dbconfiglevel],
    dbconfigpartiallevel            => [dbconfiglevel],
    dbconfiglevel                   => [fullstructurelevel],
    structuredefinedpartiallevel    => [structuredefinedlevel, structurepositionpartiallevel],
    structuredefinedlevel           => [structurepositionlevel],
    structurepositionpartiallevel   => [structurepositionlevel, structureconfigpartiallevel],
    structurepositionlevel          => [fullstructurelevel, structureconfiglevel], 
    structureconfigpartiallevel     => [structureconfiglevel],
    structureconfiglevel            => [completestructurelevel],
    fullstructurelevel              => [completestructurelevel],
    completestructurelevel          => [completestructurelevel]
)

function trim_level_gap!(chain::Vector{<: LipidAnnotationLevel})
    for l in ANNOTATIONLEVELGAPCHECK
        in(l, chain) || delete_level_above!(chain, l)
    end
    chain
end

function branch_intersect!(chain::Vector{<: LipidAnnotationLevel}, branchin, branch = branchin)
    newbranch = [x for x in branch if x in chain]
    setdiff!(chain, newbranch)
    union!(chain, intersect(newbranch, branchin))
end

function delete_level_above!(chain::Vector{<: LipidAnnotationLevel}, level::LipidAnnotationLevel)
    i = findfirst(==(level), chain)
    if !isnothing(i)
        deleteat!(chain, i)
    end
    level == completestructurelevel && return chain
    for l in ANNOTATIONLEVELMASTER[level]
        delete_level_above!(chain, l)
    end
    chain
end

# chain, like climb stair
function annotationchain(lipid::Lipid)
    chain = LipidAnnotationLevel[specieslevel]
    pilevel = testphosphatepositionlevel(lipid)
    if any(x -> ncarbonchain(x) > 1, getlipidchain(lipid))
        return union!(chain, first(pilevel))
    end
    push!(chain, molecularspecieslevel)
    snlevel = testsnpositionlevel(lipid)
    hglevel = testheadgrouppositionlevel(lipid)
    # intersect lb_n_lc
    levels = intersect((annotationchain(x) for x in getlipidbody(lipid))..., (annotationchain(x) for x in getlipidchain(lipid))...)
    # union chain_u_(lb_n_lc), 
    # omit phosphate, snposition in lb_n_lc, trim_level_equiv!
    union!(chain, levels)
    branch_intersect!(chain, pilevel...)
    branch_intersect!(chain, snlevel...)
    branch_intersect!(chain, hglevel...)
    # trim levels above gap, from high level to low level
    trim_level_gap!(chain)
end

function annotationchain(chain::CarbonChain)
    dbs = chain.doublebond
    checkdb = true
    if iszero(dbs) 
        dbs_certain = [dbpositionpartiallevel, dbconfigpartiallevel]
        dbs_result = [dbpositionlevel, dbconfiglevel]
    elseif isa(dbs, UInt8) || any(<(0x03), dbs)
        dbs_certain = LipidAnnotationLevel[]
        dbs_result = LipidAnnotationLevel[]
        checkdb = false
    elseif any(x -> iszero(x % 3), dbs)
        dbs_certain = [dbpositionpartiallevel]
        dbs_result = [dbpositionlevel]
    else
        dbs_certain = [dbpositionpartiallevel, dbconfigpartiallevel]
        dbs_result = [dbpositionlevel, dbconfiglevel]
    end
    sub = chain.substituent
    if sub isa Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}} 
        sub_result = [structuredefinedlevel, structurepositionlevel]
        checksub = true
        for (p, m) in sub
            if checksub || checkdb
                il = annotationchain(originalmolecule(m))
            else
                break 
            end
            checksub && intersect!(sub_result, il)
            checkdb && intersect!(dbs_result, il)
            checksub = checksub && length(sub_result) > 0
        end
        push!(sub_result, structuredefinedpartiallevel)
        push!(sub_result, structurepositionpartiallevel)
    elseif isnothing(sub) || isempty(sub)
        sub_result = [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel]
    else
        if any(x -> first(x) isa UnknownGroup, sub)
            sub_certain = LipidAnnotationLevel[]
            sub_result = LipidAnnotationLevel[]
            sub = filter(x -> !isa(first(x), UnknownGroup), sub)
            checksub = false
        else
            sub_certain = [structuredefinedpartiallevel]
            sub_result = [structuredefinedlevel]
            checksub = true
        end
        if checksub || checkdb
            for (m, n) in sub
                il = annotationchain(originalmolecule(m))
                checkdb && intersect!(dbs_result, il)
                checksub && intersect!(sub_result, il)
                checksub = checksub && length(sub_result) > 0
            end
        end
        union!(sub_result, sub_certain)
    end
    union!(dbs_result, dbs_certain)
    if structurepositionlevel in sub_result && dbconfiglevel in dbs_result
        union(sub_result, dbs_result, [specieslevel, molecularspecieslevel, passphosphatepositionlevel, phosphatepositionlevel, passsnpositionlevel, snpositionlevel, fullstructurelevel])
    else
        union(sub_result, dbs_result, [specieslevel, molecularspecieslevel, passphosphatepositionlevel, phosphatepositionlevel, passsnpositionlevel, snpositionlevel])
    end
end

function annotationchain(dc::Union{DehydratedChemical, Glycan})
    chain = intersect(annotationchain.(getchaincomponent(dc))...)
    ip = findall(>=(structurepositionpartiallevel), chain)
    isempty(ip) && return chain
    for ((a, b), l) in zip(IterTools.partition(getchaincomponent(dc), 2), getchainlinkage(dc))
        pa = dehydroxyposition(a)
        pb = dehydrogenposition(b)
        # glycan to check a/b
        if !ismissing(pa) && !isnothing(pa) && first(l).position == 0
            deleteat!(chain, ip)
            break
        elseif !ismissing(pb) && !isnothing(pb) && last(l).position == 0
            deleteat!(chain, ip)
            break
        end
    end
    trim_level_gap!(chain)
end

function annotationchain(c::T) where {T <: AbstractGlycan}
    if hasfield(T, :isomer)
        [
            structuredefinedlevel,
            structuredefinedpartiallevel,
            dbconfiglevel,
            dbconfigpartiallevel,
            dbpositionlevel,
            dbpositionpartiallevel,
            passsnpositionlevel,
            snpositionlevel,
            passphosphatepositionlevel,
            phosphatepositionlevel,
            molecularspecieslevel,
            specieslevel
        ]
    else
        [
            completestructurelevel,
            fullstructurelevel,
            structureconfiglevel,
            structureconfigpartiallevel,
            structurepositionlevel,
            structurepositionpartiallevel,
            structuredefinedlevel,
            structuredefinedpartiallevel,
            dbconfiglevel,
            dbconfigpartiallevel,
            dbpositionlevel,
            dbpositionpartiallevel,
            passsnpositionlevel,
            snpositionlevel,
            passphosphatepositionlevel,
            phosphatepositionlevel,
            molecularspecieslevel,
            specieslevel
        ]
    end
end
# annotationchain(c::Glycan) = intersect!(annotationchain.(getchaincomponent(c))...)
# check linkage
annotationchain(c::GlyComp) = [
    structuredefinedlevel, 
    structuredefinedpartiallevel, 
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel, 
    passsnpositionlevel, 
    snpositionlevel, 
    passphosphatepositionlevel, 
    phosphatepositionlevel,
    molecularspecieslevel, 
    specieslevel
]
annotationchain(c::T) where {S <: Nothing, T <: Monosaccharide{S}} = [
    completestructurelevel,
    fullstructurelevel,
    structureconfiglevel,
    structureconfigpartiallevel,
    structurepositionlevel,
    structurepositionpartiallevel,
    structuredefinedlevel,
    structuredefinedpartiallevel,
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel,
    passsnpositionlevel,
    snpositionlevel,
    passphosphatepositionlevel,
    phosphatepositionlevel,
    molecularspecieslevel,
    specieslevel
] 
annotationchain(c::T) where {S <: Vector{<: Pair{<: FunctionalGroup, UInt8}}, T <: Monosaccharide{<: S}} = [
    structuredefinedlevel,
    structuredefinedpartiallevel,
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel,
    passsnpositionlevel,
    snpositionlevel,
    passphosphatepositionlevel,
    phosphatepositionlevel,
    molecularspecieslevel,
    specieslevel
]
function annotationchain(c::T) where {S <: Vector{<: Pair{UInt8, <: FunctionalGroup}}, T <: Monosaccharide{<: S}}
    if any(x -> first(x) == 0, c.substituent)
        [
            structuredefinedlevel,
            structuredefinedpartiallevel,
            dbconfiglevel,
            dbconfigpartiallevel,
            dbpositionlevel,
            dbpositionpartiallevel,
            passsnpositionlevel,
            snpositionlevel,
            passphosphatepositionlevel,
            phosphatepositionlevel,
            molecularspecieslevel,
            specieslevel
        ]
    else
        [
            completestructurelevel,
            fullstructurelevel,
            structureconfiglevel,
            structureconfigpartiallevel,
            structurepositionlevel,
            structurepositionpartiallevel,
            structuredefinedlevel,
            structuredefinedpartiallevel,
            dbconfiglevel,
            dbconfigpartiallevel,
            dbpositionlevel,
            dbpositionpartiallevel,
            passsnpositionlevel,
            snpositionlevel,
            passphosphatepositionlevel,
            phosphatepositionlevel,
            molecularspecieslevel,
            specieslevel
        ] 
    end
end

annotationchain(c::AbstractPeptide) = [
    completestructurelevel,
    fullstructurelevel,
    structureconfiglevel,
    structureconfigpartiallevel,
    structurepositionlevel,
    structurepositionpartiallevel,
    structuredefinedlevel,
    structuredefinedpartiallevel,
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel,
    passsnpositionlevel,
    snpositionlevel,
    passphosphatepositionlevel,
    phosphatepositionlevel,
    molecularspecieslevel,
    specieslevel
]
annotationchain(::Union{Nothing, <: Metabolite, <: BasicCompound}) = [
    completestructurelevel,
    fullstructurelevel,
    structureconfiglevel,
    structureconfigpartiallevel,
    structurepositionlevel,
    structurepositionpartiallevel,
    structuredefinedlevel,
    structuredefinedpartiallevel,
    dbconfiglevel,
    dbconfigpartiallevel,
    dbpositionlevel,
    dbpositionpartiallevel,
    passsnpositionlevel,
    snpositionlevel,
    passphosphatepositionlevel,
    phosphatepositionlevel,
    molecularspecieslevel,
    specieslevel
]
annotationchain(::UnknownChemical) = [molecularspecieslevel, specieslevel, passphosphatepositionlevel, phosphatepositionlevel, passsnpositionlevel, snpositionlevel]

# [] -> [phosphatepositionlevel] -> [phosphatepositionlevel, passphosphatepositionlevel]
function testphosphatepositionlevel(c::Union{Phosphatidylinositol, Lysophosphatidylinositol})
    pi = first(getchaincomponent(c.backbone))
    if isnothing(pi.substituent)
        ([phosphatepositionlevel], 
        [phosphatepositionlevel, passphosphatepositionlevel])
    else
        n = 0
        p = String[]
        for i in first(getchaincomponent(c.backbone)).substituent
            if i isa Phosphate
                n += 1
            elseif i isa Pair && last(i) isa Phosphate
                n += 1
                push!(p, string(Int(first(i)), "'"))
            end
        end
        if n < 1
            ([phosphatepositionlevel, passphosphatepositionlevel], 
            [phosphatepositionlevel, passphosphatepositionlevel])
        else
            (length(p) == n ? [phosphatepositionlevel] : LipidAnnotationLevel[], 
            [phosphatepositionlevel, passphosphatepositionlevel])
        end
    end
end
testphosphatepositionlevel(::Lipid) = (
    [phosphatepositionlevel, passphosphatepositionlevel], 
    [phosphatepositionlevel, passphosphatepositionlevel]
)
# [] -> [snpositionlevel] -> [snpositionlevel, passsnpositionlevel]
function testsnpositionlevel(lipid::Union{<: Glycerolipid, <: Glycerophospholipid})
    (any(iszero, decode_sn(lipid)) ? LipidAnnotationLevel[] : [snpositionlevel],
    [snpositionlevel, passsnpositionlevel])
end
testsnpositionlevel(::Lipid) = ([snpositionlevel, passsnpositionlevel], [snpositionlevel, passsnpositionlevel])
testsnpositionlevel(::GlycerophosphoNacylethanolamine) = ([snpositionlevel, passsnpositionlevel], [snpositionlevel, passsnpositionlevel])
# [structuredefinedpartiallevel, structuredefinedlevel] -> [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel]
testheadgrouppositionlevel(::Lipid) = (
    [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel],
    [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel]
)
function testheadgrouppositionlevel(lipid::SphingoBone)
    (isnothing(lipid.headgroup) ? [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel] : 
        iszero(lipid.position) ? [structuredefinedpartiallevel, structuredefinedlevel] : 
        [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel],
        [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel]
    )
end
function testheadgrouppositionlevel(lipid::MixSphingoBone)
    (isnothing(lipid.headgroup) ? [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel] : 
        any(iszero, lipid.position) ? [structuredefinedpartiallevel, structuredefinedlevel] : 
        [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel],
        [structuredefinedpartiallevel, structuredefinedlevel, structurepositionpartiallevel, structurepositionlevel]
    )
end

function annotationlevel(c; partial = false, additional = false, pass = false)
    chain = annotationchain(c)
    if !partial
        chain = transform_partial.(chain)
    end
    if !additional
        chain = transform_additional.(chain)
    end
    if !pass
        i = findall(x -> isa(x, PassLevel), chain)
        j = [findfirst(y -> ==(x.level, y.level), chain) for x in @view chain[i]]
        deleteat!(chain, sort!(unique!(vcat(i, filter!(!isnothing, j)))))
    end
    maximal_annotationlevel(chain)
end
function maximal_annotationlevel(level::Vector{T}) where {T <: LipidAnnotationLevel}
    newlevel = T[]
    for l in level 
        if l in newlevel || any(>(l), newlevel)
            continue
        else
            del = findall(<(l), newlevel)
            deleteat!(newlevel, del)
            push!(newlevel, l)
        end
    end
    newlevel
end

# depre
#==
function partialize(level::AnnotationLevel)
    @match level begin
        specieslevel                    => specieslevel
        molecularspecieslevel           => specieslevel
        phosphatepositionlevel          => specieslevel
        passphosphatepositionlevel      => specieslevel
        snpositionlevel                 => molecularspecieslevel
        passsnpositionlevel             => molecularspecieslevel
        dbpositionpartiallevel          => dbpositionpartiallevel
        dbpositionlevel                 => dbpositionpartiallevel
        dbconfigpartiallevel            => dbconfigpartiallevel
        dbconfiglevel                   => dbconfigpartiallevel
        structuredefinedpartiallevel    => structuredefinedpartiallevel # headgroup comp/chain mod comp
        structuredefinedlevel           => structuredefinedpartiallevel # all comp
        structurepositionpartiallevel   => structurepositionpartiallevel # headgroup position/chain mod position
        structurepositionlevel          => structurepositionpartiallevel # all position
        structureconfigpartiallevel     => structureconfigpartiallevel
        structureconfiglevel            => structureconfigpartiallevel
        _                               => throw(ArgumentError("$level cannot be partialized"))
    end
end

function maximal_annotationlevel_sorted(level::Vector{AnnotationLevel})
    newlevel = AnnotationLevel[]
    for l in level
        while !isempty(newlevel) && l > last(newlevel)
            pop!(newlevel)
        end
        push!(newlevel, l)
    end
    newlevel
end
function maximal_annotationlevel(level::Vector{AnnotationLevel})
    newlevel = AnnotationLevel[]
    for l in level 
        if l in newlevel || any(>(l), newlevel)
            continue
        else
            del = findall(<(l), newlevel)
            deleteat!(newlevel, del)
            push!(newlevel, l)
        end
    end
    newlevel
end
function minimal_annotationlevel(level::Vector{<: LipidAnnotationLevel})
    newlevel = AnnotationLevel[]
    pending = AnnotationLevel[]
    for l in level 
        if l in newlevel || any(<(l), newlevel)
            continue
        elseif l isa LeastLevel
            i = findall(>(l.level), newlevel)
            if isnothing(i)
                push!(newlevel, l.level)
            else
                push!(pending, l.level)
            end
        else
            del = findall(>(l), newlevel)
            deleteat!(newlevel, del)
            push!(newlevel, l)
        end
    end
    for l in pending
        if all(!>(l), newlevel)
            push!(newlevel, l)
        end
    end
    newlevel
end
# function contain(l1::AnnotationLevel, l2::AnnotationLevel)
#     if l1 == completestructurelevel
#         l2 == fullstructurelevel || l2 == structureconfiglevel || contain(fullstructurelevel, l2)
#     elseif l1 == fullstructurelevel
#         l2 == snpositionlevel || l2 == passsnpositionlevel || l2 == structurepositionlevel || l2 == dbconfiglevel
#     else
#         false
#     end
# end

# # push not exclude highest
# function push_annotationlevel(level::Vector{AnnotationLevel}, l::AnnotationLevel; highest = l)
#     if !(highest > l)
#     elseif any(x -> l > x || contain(x, l), level) 
#         return level
#     else
#         minimal_annotationlevel(vcat(level, l))
#     end
# end
# exact level
function annotationlevel(lipid::Lipid; partial = false, sublevel = false)
    # testphosphatepositionlevel
    pilevel = testphosphatepositionlevel(lipid)
    any(x -> ncarbonchain(x) > 1, getlipidchain(lipid)) && return sublevel ? [pilevel.level] : [transform_sublevel(pilevel.level)]
    lb = map(x -> annotationlevel(x; partial = false, sublevel = true), getlipidbody(lipid))
    lb = unique!(vcat(lb...))
    # check snpositionlevel
    snlevel = testsnpositionlevel(lipid)
    any(==(molecularspecieslevel), lb) && return sublevel ? [snlevel.level] : [transform_sublevel(snlevel.level)]
    # mix findmin
    lc = map(x -> annotationlevel(x; partial = true, sublevel = true), getlipidchain(lipid))
    if !partial
        lc = map(lc) do l
            maximal_annotationlevel(transform_partial.(l))
        end
    end
    lc = unique!(vcat(lc...))
    any(==(molecularspecieslevel), lc) && return sublevel ? [snlevel.level] : [transform_sublevel(snlevel.level)]
    l = vcat(lc, lb, snlevel, testheadgrouppositionlevel(lipid))
    result = minimal_annotationlevel(l)
    # result = push_annotationlevel(result, snlevel; highest = snpositionlevel)
    # result = push_annotationlevel(result, testheadgrouppositionlevel(lipid); highest = structurepositionlevel)
    sublevel ? result : maximal_annotationlevel(transform_sublevel.(result))
end

# assume snpositionlevel/passsnpositionlevel, phosphatepositionlevel/passphosphatepositionlevel 
function annotationlevel(chain::CarbonChain; partial = false, sublevel = false)
    # molecularspecieslevel dbpositionlevel dbconfiglevel structuredefinedlevel structurepositionlevel fullstructurelevel
    # implement structureconfigpartiallevel structureconfiglevel completestructurelevel after R/S available
    dbs = chain.doublebond
    if iszero(dbs) 
        dbs_result = [dbpositionlevel, dbconfiglevel]
        checkdbs = [true, true]
    elseif isa(dbs, UInt8) || any(<(0x03), dbs)
        dbs_result = [molecularspecieslevel]
        checkdbs = [false]
    elseif any(x -> iszero(x % 3), dbs)
        dbs_result = [dbpositionlevel]
        checkdbs = [true]
    else
        dbs_result = [dbpositionlevel, dbconfiglevel]
        checkdbs = [true, true]
    end

    checkdbs = dbs_result .> molecularspecieslevel
    sub = chain.substituent
    if chain.substituent isa Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}} 
        sub_result = [structuredefinedlevel, structurepositionlevel]
        checksub = [true, true]
        for (p, m) in chain.substituent
            if any(checksub) || any(checkdbs)
                il = annotationlevel(originalmolecule(m); partial = false, sublevel = true)
            else
                break 
            end
            for i in eachindex(sub_result)
                if checksub[i] && any(<(sub_result[i]), il)
                    sub_result[i:end] .= partialize.(sub_result[i:end])
                    checksub[i:end] .= false
                end
            end
            for i in eachindex(dbs_result)
                if checkdbs[i] && any(<(dbs_result[i]), il)
                    dbs_result[i:end] .= partialize.(dbs_result[i:end])
                    checkdbs[i:end] .= false
                end
            end
        end
    elseif isnothing(sub) || isempty(sub)
        sub_result = [structurepositionlevel]
    else
        if any(x -> first(x) isa UnknownGroup, sub)
            sub_result = [molecularspecieslevel]
            checksub = false
        else
            sub_result = [structuredefinedlevel]
            checksub = true
        end
        if checksub || any(checkdbs)
            for (m, n) in chain.substituent
                il = annotationlevel(originalmolecule(m); partial = false, sublevel = true)
                if checksub && any(<(structuredefinedlevel), il)
                    sub_result = [structuredefinedpartiallevel]
                    checksub = false
                end
                for i in eachindex(dbs_result)
                    if checkdbs[i] && any(<(dbs_result[i]), il)
                        dbs_result[i:end] .= partialize.(dbs_result[i:end])
                        checkdbs[i:end] .= false
                    end
                end
            end
        end
    end
    if partial 
        sub_result = maximal_annotationlevel_sorted(sub_result)
        dbs_result = maximal_annotationlevel_sorted(dbs_result)
    else
        i = findfirst(!ispartial, reverse(sub_result))
        sub_result = isnothing(i) ? [molecularspecieslevel] : [sub_result[end - i + 1]]
        i = findfirst(!ispartial, reverse(dbs_result))
        dbs_result = isnothing(i) ? [molecularspecieslevel] : [dbs_result[end - i + 1]]
    end
    (sub_result == [structurepositionlevel] && dbs_result == [dbconfiglevel]) && return [fullstructurelevel] 
    if !sublevel
        sub_result = map(transform_sublevel, sub_result)
        dbs_result = map(transform_sublevel, dbs_result)
    end
    maximal_annotationlevel(vcat(sub_result, dbs_result))
end

function annotationlevel(dc::DehydratedChemical; partial = false, sublevel = true)
    lv = unique!(vcat(annotationlevel.(getchaincomponent(dc); partial, sublevel)...))
    for ((a, b), l) in zip(IterTools.partition(getchaincomponent(dc), 2), getchainlinkage(dc))
        pa = dehydroxyposition(a)
        pb = dehydrogenposition(b)
        if !ismissing(pa) && !isnothing(pa)
            if first(l).position == 0
                push!(lv, structuredefinedlevel)
            elseif a isa Saccharide && first(l) isa Linkageposition
                push!(lv, fullstructurelevel)
            end
        end
        if !ismissing(pb) && !isnothing(pb) && last(l).position == 0
            push!(lv, structuredefinedlevel)
        end
    end
    unique!(lv)
    if !sublevel
        lv = transform_sublevel.(lv)
    end
    minimal_annotationlevel(lv)
end
function annotationlevel(dc::Glycan; partial = false, sublevel = true)
    lv = unique!(vcat(annotationlevel.(getchaincomponent(dc))...))
    for ((a, b), l) in zip(IterTools.partion(getchaincomponent(dc)), getchainlinkage(dc))
        pa = dehydroxyposition(a)
        pb = dehydrogenposition(b)
        if !ismissing(pa) && !isnothing(pa)
            if first(l).position == 0
                push!(lv, structuredefinedlevel)
            elseif first(l) isa Linkageposition
                push!(lv, fullstructurelevel)
            end
        end
        if !ismissing(pb) && !isnothing(pb) && last(l).position == 0
            push!(lv, structuredefinedlevel)
        end
    end
    if !partial
        filter!(!ispartial, lv)
    end
    unique!(lv)
    if !sublevel
        lv = transform_sublevel.(lv)
    end
    minimal_annotationlevel(lv)
end
annotationlevel(c::AbstractGlycan; partial = false, sublevel = true) = [completestructurelevel]
annotationlevel(c::GlyComp; partial = false, sublevel = true) = [structuredefinedlevel]
annotationlevel(c::Monosaccharide; partial = false, sublevel = true) = [completestructurelevel] # check sub position

annotationlevel(c::AbstractPeptide; partial = false, sublevel = true) = [completestructurelevel]
annotationlevel(::Metabolite; partial = false, sublevel = true) = [completestructurelevel]
annotationlevel(::BasicCompound; partial = false, sublevel = true) = [completestructurelevel]
annotationlevel(::UnknownChemical; partial = false, sublevel = true) = [molecularspecieslevel]

function isstructuredefinedlevel(lipid::Lipid)
    overmolecularspecieslevel(lipid) || return false
    all(isstructuredefinedlevel, getlipidbody(lipid)) ? !any(!isstructuredefinedlevel, getlipidchain(lipid)) : false
    # one not full, all over sd
end

function isstructuredefinedlevel(chain::CarbonChain)
    sub = chain.substituent
    isfullstructurelevel(chain) && return false
    for (f, n) in sub
        f isa UnknownGroup && return false # specieslevel
    end
    true
end

function isstructuredefinedlevel(c::AbstractChemical)
    # check linkage
    all(isstructuredefinedlevel, getchaincomponent(c))
end

isstructuredefinedlevel(c::Saccharide) = true # check linkage
isstructuredefinedlevel(c::AbstractPeptide) = true # check linkage
isstructuredefinedlevel(::Metabolite) = true
isstructuredefinedlevel(::BasicCompound) = true

function issubstructuredefinedlevel(chain::CarbonChain)
    sub = chain.substituent
    issublevelstructurelevel(chain) && return false
    for (f, n) in sub
        f isa UnknownGroup && return false # specieslevel
    end
    true
end

function issublevelstructurelevel(chain::CarbonChain)
    sub = chain.substituent
    isnothing(sub) || isempty(sub) || chain.substituent isa Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}} 
end

function isfullstructurelevel(c::AbstractChemical)
    chain = getchaincomponent(c)
    if all(isfullstructurelevel, chain)
        length(chain) == 1 && return true
        default_lk = [makelinkage(DehydratedChemical, a, b) for (a, b) in IterTools.partition(chain, 2, 1)]
        real_lk = getchainlinkage(c)
        isnothing(real_lk) && return false
        for (d, r) in zip(default_lk, real_lk)
            dl, dr = d 
            rl, rr = r 
            dl == lk(nothing) || iszero(rl.position) && return false
            dr == lk(nothing) || iszero(rr.position) && return false
        end
        return true
    end
    false
end

function isfullstructurelevel(lipid::Lipid)
    overmolecularspecieslevel(lipid) || return false
    all(isfullstructurelevel, getlipidbody(lipid)) ? !any(!isfullstructurelevel, getlipidchain(lipid)) : false
end

function isfullstructurelevel(chain::CarbonChain)
    sub = chain.substituent
    if isnothing(sub) || isempty(sub) || chain.substituent isa Vector{<: Pair{UInt8, <: AbstractFunctionalGroup}} 
        dbs = chain.doublebond
        iszero(dbs) && return true
        isa(dbs, UInt8) && return false
        return !any(x -> iszero(x % 3), dbs)
    end
    false
end

isfullstructurelevel(c::Monosaccharide) = true
isfullstructurelevel(c::αAminoAcid) = true
isfullstructurelevel(::Metabolite) = true
isfullstructurelevel(::BasicCompound) = true

function isspecieslevel(chain::CarbonChain)
    sub = chain.substituent
    isfullstructurelevel(chain) && return false
    for (f, n) in sub
        f isa UnknownGroup && return true # specieslevel
    end
    false
end

function issnpositionlevel(lipid::Union{<: Glycerolipid, <: Glycerophospholipid})
    overmolecularspecieslevel(lipid) || return false
    position = decode_sn(lipid)
    !any(==(0), position)
    # one not full, all over sn
end

function ismolecularspecieslevel(lipid::Lipid)
    !any(x -> ncarbonchain(x) > 1, getlipidchain(lipid))
end

function isdbpositionlevel(c::AbstractChemical)
    # check linkage
    all(isdbpositionlevel, getchaincomponent(c))
end

function isdbpositionlevel(lipid::Lipid)
    overmolecularspecieslevel(lipid) || return false
    all(isdbpositionlevel, getlipidbody(lipid)) ? !any(!isdbpositionlevel, getlipidchain(lipid)) : false
    # one not full, all over db
end

function isdbpositionlevel(chain::CarbonChain)
    dbs = chain.doublebond
    sub = chain.substituent
    isfullstructurelevel(chain) && return false
    for (f, n) in sub
        f isa UnknownGroup && return false # specieslevel
    end
    iszero(dbs) || !isa(dbs, UInt8)
end

isdbpositionlevel(c::Saccharide) = true
isdbpositionlevel(c::AbstractPeptide) = true
isdbpositionlevel(::Metabolite) = true
isdbpositionlevel(::BasicCompound) = true


function overmolecularspecieslevel(lipid::Lipid)
    !any(x -> ncarbonchain(x) > 1, getlipidchain(lipid))
end
=#
