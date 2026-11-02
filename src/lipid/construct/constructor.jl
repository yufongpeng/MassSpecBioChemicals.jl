function defaultchemicalkwargs(::Type{T}, core...) where {T <: Lipid}
    name = lipiddefaultname(T, core...)
    elements = lipiddefaultelements(T, core...)
    (; name, elements)
end

function make_lipid(::Type{T}, core::Tuple; new_attr = true, kwargs...) where {T <: Lipid}
    S = Core.apply_type(T, typeof(core[begin]), typeof(core[begin+1]))
    kw = (; pairs(defaultchemicalkwargs(S, core...))..., kwargs...)
    attr = Base.structdiff(kw, NamedTuple{(:name, :elements)})
    S(core..., String(kw.name), collect(Pair{String, Int}, kw.elements),
      buildattribute(new_attr, pairs(attr)))
end

# Rewrap the unparameterized type since sorting may change the type of the chains
function make_radylglycerol(::Type{T}, backbone, chain, sn; kwargs...) where {T <: Lipid}
    chain, code = arrange_lipid_chains(T, chain, sn)
    make_lipid(T, check_structure!(T, backbone, chain, code); kwargs...)
end

CycloChain(position, atom, doublebond, substituent, chirality;
        name = string("CycloChain(", atom, ")"), elements = carbonchaindefaultelements()) =
    CycloChain(position, atom, doublebond, substituent, chirality, String(name),
        collect(Pair{String, Int}, elements))

IsoprenoidChain{N}(isotopiclabel = nothing; name = string("IsoprenoidChain{", N, "}"),
        elements = carbonchaindefaultelements()) where {N} =
    IsoprenoidChain{N, typeof(isotopiclabel)}(isotopiclabel, String(name),
        collect(Pair{String, Int}, elements))

MonoFattyAcyl(backbone, chain; kwargs...) =
    make_lipid(MonoFattyAcyl, check_structure!(MonoFattyAcyl, backbone, chain); kwargs...)

NacylAmine(backbone, chain; kwargs...) =
    make_lipid(NacylAmine, check_structure!(NacylAmine, backbone, chain); kwargs...)

FattyAcylEster(backbone, chain; position = 0x00, kwargs...) =
    make_lipid(FattyAcylEster, check_structure!(FattyAcylEster, backbone, chain, encode_position(position)); kwargs...)
FattyAcylEster(backbone, chain, position; kwargs...) =
    make_lipid(FattyAcylEster, check_structure!(FattyAcylEster, backbone, chain, encode_position(position)); kwargs...)

# Glycerolipids and glycerophospholipids: `sn` can be the merged code, a symbol, a string like "(sn-1)",
# or a collection of them for each chain; see `arrange_lipid_chains`.
Monoradylglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Monoradylglycerol, backbone, chain, sn; kwargs...)
Monoradylglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Monoradylglycerol, backbone, chain, sn; kwargs...)

Diradylglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Diradylglycerol, backbone, chain, sn; kwargs...)
Diradylglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Diradylglycerol, backbone, chain, sn; kwargs...)

Triradylglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Triradylglycerol, backbone, chain, sn; kwargs...)
Triradylglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Triradylglycerol, backbone, chain, sn; kwargs...)

Omodifiedradylglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Omodifiedradylglycerol, backbone, chain, sn; kwargs...)
Omodifiedradylglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Omodifiedradylglycerol, backbone, chain, sn; kwargs...)

Radyldiglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Radyldiglycerol, backbone, chain, sn; kwargs...)
Radyldiglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Radyldiglycerol, backbone, chain, sn; kwargs...)

Radyltriglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Radyltriglycerol, backbone, chain, sn; kwargs...)
Radyltriglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Radyltriglycerol, backbone, chain, sn; kwargs...)

Radylglycerophosphate(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Radylglycerophosphate, backbone, chain, sn; kwargs...)
Radylglycerophosphate(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Radylglycerophosphate, backbone, chain, sn; kwargs...)

Radylglycerophosphoaminoacid(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Radylglycerophosphoaminoacid, backbone, chain, sn; kwargs...)
Radylglycerophosphoaminoacid(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Radylglycerophosphoaminoacid, backbone, chain, sn; kwargs...)

Radylglycerophosphoglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Radylglycerophosphoglycerol, backbone, chain, sn; kwargs...)
Radylglycerophosphoglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Radylglycerophosphoglycerol, backbone, chain, sn; kwargs...)

Bisradylglycerophosphate(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Bisradylglycerophosphate, backbone, chain, sn; kwargs...)
Bisradylglycerophosphate(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Bisradylglycerophosphate, backbone, chain, sn; kwargs...)

Bisradylglycerophosphoglycerol(backbone, chain; sn = :noorder, kwargs...) =
    make_radylglycerol(Bisradylglycerophosphoglycerol, backbone, chain, sn; kwargs...)
Bisradylglycerophosphoglycerol(backbone, chain, sn; kwargs...) =
    make_radylglycerol(Bisradylglycerophosphoglycerol, backbone, chain, sn; kwargs...)

GlycerophosphoNacylethanolamine(backbone, chain; kwargs...) =
    make_lipid(GlycerophosphoNacylethanolamine, check_structure!(GlycerophosphoNacylethanolamine, backbone, chain, nothing); kwargs...)

SphingoBone(headgroup, chain, position, chirality; kwargs...) =
    make_lipid(SphingoBone, check_structure!(SphingoBone, headgroup, chain, encode_position_chirality(SphingoBone, headgroup, position, chirality)); kwargs...)
SphingoBone(headgroup, chain; position = nothing, chirality = nothing, kwargs...) = 
    make_lipid(SphingoBone, check_structure!(SphingoBone, headgroup, chain, encode_position_chirality(SphingoBone, headgroup, position, chirality)); kwargs...)
SphingoBone(headgroup, chain, position; chirality = nothing, kwargs...) =
    make_lipid(SphingoBone, check_structure!(SphingoBone, headgroup, chain, encode_position_chirality(SphingoBone, headgroup, position, chirality)); kwargs...)

MixSphingoBone(headgroup, chain, position, chirality; kwargs...) =
    make_lipid(MixSphingoBone, check_structure!(MixSphingoBone, headgroup, chain, encode_position_chirality(MixSphingoBone, headgroup, position, chirality)); kwargs...)
MixSphingoBone(headgroup, chain; position = nothing, chirality = nothing, kwargs...) = 
    make_lipid(MixSphingoBone, check_structure!(MixSphingoBone, headgroup, chain, encode_position_chirality(MixSphingoBone, headgroup, position, chirality)); kwargs...)
MixSphingoBone(headgroup, chain, position; chirality = nothing, kwargs...) =
    make_lipid(MixSphingoBone, check_structure!(MixSphingoBone, headgroup, chain, encode_position_chirality(MixSphingoBone, headgroup, position, chirality)); kwargs...)

SterolBone(chain; kwargs...) = make_lipid(SterolBone, (chain,); kwargs...)
SubstitutedSterol(backbone; kwargs...) =
    make_lipid(SubstitutedSterol,(backbone,); kwargs...)
Retinylester(chain; kwargs...) = make_lipid(Retinylester, (chain,); kwargs...)
CoenzymeQ(chain; kwargs...) = make_lipid(CoenzymeQ, (chain,); kwargs...)
