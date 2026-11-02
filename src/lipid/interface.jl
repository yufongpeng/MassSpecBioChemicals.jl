# MassSpecBioChemicals
parentchemical(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
parentchemical(chain::CarbonChain{Acyl}) = FattyAcid(chain)
parentchemical(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00, AChirality())
leavinggroup(::CarbonChain{Acyl}) = Dehydroxy()
leavinggroup(::CarbonChain{SPB}) = Dehydrogen()
parentchemical(chain::CarbonChain{<: AbstractSTRing}) = SterolBone(chain)
leavinggroup(::CarbonChain{<: AbstractSTRing}) = Dehydrogen()

dehydroxyposition(::FattyAcid) = nothing
dehydrogenposition(::FattyAcid) = missing
dehydroxyposition(::FattyAlcohol) = nothing
dehydrogenposition(::FattyAlcohol) = nothing
dehydroxyposition(::FattyAmine) = nothing
dehydrogenposition(::FattyAmine) = nothing

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
chainposition(::Type{<: Radyldiglycerol}) = ["sn-1", "sn-2", "sn-1'", "sn-2'"]
chainposition(::Type{<: Bisradylglycerophosphate}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
chainposition(::Type{<: Radyltriglycerol}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
chainposition(lipid::T) where {T <: Lipid} = chainposition(T)
nchainposition(::Type{T}) where {T <: Lipid} = length(chainposition(T))
nchainposition(lipid::T) where {T <: Lipid} = length(chainposition(T))
ncarbonchain(lipid::Lipid) = ncarbonchain(lipid.chain)

hascycloheadgroup(::Lipid) = false
hascycloheadgroup(::Sphingolipid) = false
hascycloheadgroup(::CeramidePhosphate) = true
hascycloheadgroup(::SphingoidBasePhosphate) = true

getlipidchain(lipid::Lipid) = tuplize(lipid.chain)
getlipidbody(lipid::Lipid) = tuplize(lipid.backbone)
getlipidbody(lipid::Sphingolipid) = tuplize(lipid.headgroup)

chiralitize(::Type{T}, args, kwargs) where {T <: Lipid} = T, args, kwargs

ncarbonchain(chains::Tuple) = sum(ncarbonchain, chains)
ncarbonchain(::CarbonChain{<: CarbonChainType}) = 1
ncarbonchain(::CarbonChain{<: T}) where {T <: Tuple} = length(T.parameters)
ncarbon(chains::Tuple) = sum(ncarbon, chains)
ncarbon(chain::CarbonChain) = chain.carbon
ndoublebond(chain::CarbonChain{S, UInt8}) where S = Int(chain.doublebond)
ndoublebond(chain::CarbonChain{S}) where S = length(chain.doublebond)

deisomerize(lipid::MonoFattyAcyl) = MonoFattyAcyl(deisomerize(lipid.backbone), deisomerize(lipid.chain))
deisomerize(lipid::NacylAmine) = NacylAmine(deisomerize(lipid.backbone), deisomerize(lipid.chain))
deisomerize(lipid::FattyAcylEster) = FattyAcylEster(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00)
deisomerize(lipid::Monoradylglycerol) = Monoradylglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, RSChirality())
deisomerize(lipid::Diradylglycerol) = Diradylglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, RSChirality())
deisomerize(lipid::Triradylglycerol) = Triradylglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, RSChirality())
deisomerize(lipid::Radyldiglycerol) = Radyldiglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x0000, (RSChirality(), RSChirality()))
deisomerize(lipid::Radyltriglycerol) = Radyltriglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x0000, (RSChirality(), RSChirality(), RSChirality()))
deisomerize(lipid::Radylglycerophosphate) = Radylglycerophosphate(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, RSChirality())
deisomerize(lipid::Radylglycerophosphoaminoacid) = Radylglycerophosphoaminoacid(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, (RSChirality(), RSChirality()))
deisomerize(lipid::Radylglycerophosphoglycerol) = Radylglycerophosphoglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x00, (RSChirality(), RSChirality()))
deisomerize(lipid::Bisradylglycerophosphate) = Bisradylglycerophosphate(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x0000, (RSChirality(), RSChirality()))
deisomerize(lipid::Bisradylglycerophosphoglycerol) = Bisradylglycerophosphoglycerol(deisomerize(lipid.backbone), deisomerize(lipid.chain), 0x0000, (RSChirality(), RSChirality(), RSChirality()))
deisomerize(lipid::SphingoBone) = SphingoBone(deisomerize(lipid.headgroup), deisomerize(lipid.chain), 0x00, lipid.chirality isa RSSystem ? RSChirality() : map(x -> RSChirality(), lipid.chirality))
deisomerize(lipid::MixSphingoBone) = MixSphingoBone(deisomerize(lipid.headgroup), deisomerize(lipid.chain), map(x -> 0x00, lipid.position), lipid.chirality isa Vector{<: RSSystem} ? map(x -> RSChirality(), lipid.chirality) : map(x -> map(y -> RSChirality(), x), lipid.chirality))
