dissociate_headgroup(lipid::Glycerolipid) = lipid
function dissociate_headgroup(lipid::Union{<: Glycerophospholipid, <: Omodifiedradylglycerol})
    sn = decode_sn(lipid)
    n = ncarbonchain(lipid)
    cls = @match n begin
        1   => Monoradylglycerol
        2   => Diradylglycerol
        3   => Triradylglycerol
    end
    sn = sn .* ((nchainposition(cls) + 1) .^ ((n - 1):-1:0))
    cls(last(getchaincomponent(lipid.backbone)), lipid.chain, UInt8(sum(sn)))
end

dissociate_headgroup(lipid::CeramideBone) = SphingoBone(nothing, lipid.chain, 0x00)
dissociate_headgroup(lipid::SphingoidBaseBone) = SphingoBone(nothing, lipid.chain, 0x00)
dissociate_headgroup(lipid::Ceramide) = lipid
dissociate_headgroup(lipid::SphingoidBase) = lipid

# to INCHI, Search lipid map?
# fatty acyl to INCHI, Search lipid map?