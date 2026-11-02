"""
    make_lipid(Con, bone, pos, chain, sn)

Construct lipids from parsed information
"""
make_lipid(Con::Type{<: FattyAcyl}, bone, pos, chain, sn) = 
    make_lipid(Con, check_structure!(Con, bone, first(chain)))

function make_lipid(Con::Type{<: NacylAmine}, bone, pos, chain, sn)
    if bone isa Tuple
        if length(chain) == 1
            hbone = parse_lipid("FN 0:0")
            make_lipid(Con, check_structure!(Con, hbone, first(chain)))
        elseif length(chain) == nchainposition(Con)        
            bhead, bbone, bchain = bone
            hbone = make_lipid(bhead, bbone, nothing, [first(chain)], 0x00)
            make_lipid(Con, check_structure!(Con, hbone, last(chain)))
        else
            throw(ArgumentError("Invalid number of chains"))
        end
    else
        make_lipid(Con, check_structure!(Con, bone, first(chain)))
    end
end

function make_lipid(Con::Type{<: FattyAcylEster}, bone, pos, chain, sn)
    bhead, bbone, bechain = bone
    if length(chain) == 1
        if bechain == (Acyl, )
            hbone = parse_lipid("FA 0:0")
        else
            hbone = parse_lipid("FOH 0:0")
        end
        make_lipid(Con, check_structure!(Con, hbone, first(chain), sn))
    elseif length(chain) == nchainposition(Con)     
        hbone = make_lipid(bhead, bbone, nothing, [first(chain)], 0x00)
        make_lipid(Con, check_structure!(Con, hbone, last(chain), sn))
    else
        throw(ArgumentError("Invalid number of chains"))
    end
end

# Radylglycerophosphoglycerol Single -> Tuple
make_lipid(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, bone, pos, chain, sn) = 
    make_lipid(Con, check_structure!(Con, bone, chain, sn))

make_lipid(Con::Type{<: Sphingolipid}, bone, pos, chain, sn) = 
    make_lipid(Con, check_structure!(Con, bone, chain, pos...))

function make_lipid(Con::Type{<: Sterol}, bone, pos, chain, sn)
    throw(ArgumentError("`make_lipid` not implemented for `Sterol`"))
end

function make_lipid(Con::Type{<: Prenol}, bone, pos, chain, sn)
    throw(ArgumentError("`make_lipid` not implemented for `Prenol`"))
end