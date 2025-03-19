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
                il = annotationchain(parentchemical(m))
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
                il = annotationchain(parentchemical(m))
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
function annotationchain(c::T) where {S <: Vector{<: Pair{<: FunctionalGroup, UInt8}}, T <: Monosaccharide{<: S}}
    if any(x -> first(x) == Phosphate(), c.substituent)
        [
        structuredefinedlevel,
        structuredefinedpartiallevel,
        dbconfiglevel,
        dbconfigpartiallevel,
        dbpositionlevel,
        dbpositionpartiallevel,
        passsnpositionlevel,
        snpositionlevel,
        molecularspecieslevel,
        specieslevel
    ]
    else
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
    end 
end

function annotationchain(c::T) where {S <: Vector{<: Pair{UInt8, <: FunctionalGroup}}, T <: Monosaccharide{<: S}}
    i = findall(x -> first(x) == 0, c.substituent)
    if isempty(i)
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
    elseif any(j -> last(c.substituent[j]) == Phsphate(), i)
        [
            structuredefinedlevel,
            structuredefinedpartiallevel,
            dbconfiglevel,
            dbconfigpartiallevel,
            dbpositionlevel,
            dbpositionpartiallevel,
            passsnpositionlevel,
            snpositionlevel,
            molecularspecieslevel,
            specieslevel
        ]
    else
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
