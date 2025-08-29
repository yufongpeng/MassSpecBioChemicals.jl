function class_struct_HC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `HC` is not suppoted"))
    isnothing(sil) || @warn "Ignore stable isotope labeling on class"
    (MonoFattyAcyl, makechemical(Dihydrogen(); sil), (Alkyl, ))
end
function class_struct_FA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FA` is not suppoted"))
    (MonoFattyAcyl, makechemical(HydrogenOxide(); sil), (Acyl, ))
end
function class_struct_FAL(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FAL` is not suppoted"))
    isnothing(sil) || @warn "Ignore stable isotope labeling on class"
    (MonoFattyAcyl, makechemical(Dihydrogen(); sil), (Acyl, ))
end
function class_struct_FOH(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FOH` is not suppoted"))
    (MonoFattyAcyl, makechemical(HydrogenOxide(); sil), (Alkyl, ))
end
function class_struct_FAM(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FAM` is not suppoted"))
    (MonoFattyAcyl, makechemical(Ammonia(); sil), (Acyl, ))
end
function class_struct_FN(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FN` is not suppoted"))
    (MonoFattyAcyl, makechemical(Ammonia(); sil), (Alkyl, ))
end
function class_struct_CAR(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `CAR` is not suppoted"))
    (MonoFattyAcyl, makechemical(chiralchemical(Carnitine, parse_singlechirality(pre)); sil), (Acyl, ))
end
function class_struct_CoA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `CoA` is not suppoted"))
    (MonoFattyAcyl, makechemical(CoA(); sil), (Acyl, ))
end
function class_struct_NAE(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `NAE` is not suppoted"))
    (NacylAmine, makechemical(DehydratedChemical, Ethanolamine(); sil), (Acyl, ))
end
function class_struct_NAT(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `NAT` is not suppoted"))
    (NacylAmine, makechemical(Taurine(); sil), (Acyl, ))
end
function class_struct_NAAA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `NA` is not suppoted"))
    # pre
    (NacylAmine, makechemical(parse_aa(string(pre, replace(cls, r"^NA" => ""))); sil), (Acyl, ))
end
function class_struct_NA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `NA` is not suppoted"))
    (NacylAmine, class_struct_FN(head, pre, cls, post, pos, sil), (Alkyl, Acyl))
end
function class_struct_WE(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `WE` is not suppoted"))
    (FattyAcylEster, class_struct_FOH(head, pre, cls, post, pos, sil), (Alkyl, Acyl))
end
function class_struct_FAHFA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `FAHFA` is not suppoted"))
    isnothing(sil) || @warn "Ignore stable isotope labeling on class"
    (FattyAcylEster, class_struct_FA(head, pre, cls, post, pos, sil), (Acyl, Acyl)) # second one ;O
end
function class_struct_MG(head, pre, cls, post, pos, sil)
    gly = parse_glycerol_rs(pre)
    (isnothing(head) ? Monoradylglycerol : Omodifiedradylglycerol, 
        isnothing(head) ? makechemical(gly; sil) : 
        concatchemical(DehydratedChemical, parse_headgroup(head), gly; sil), # makechemical => concatchemical? 
        (Radyl, )
    )
end
function class_struct_DG(head, pre, cls, post, pos, sil)
    gly = parse_glycerol_rs(pre)
    (isnothing(head) ? Diradylglycerol : Omodifiedradylglycerol, 
        isnothing(head) ? makechemical(gly; sil) : 
        concatchemical(DehydratedChemical, parse_headgroup(head), gly; sil), # makechemical => concatchemical?
        (Radyl, Radyl)
    )
end
function class_struct_TG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `TG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Triradylglycerol, makechemical(gly; sil), (Radyl, Radyl, Radyl))
end
function class_struct_SQMG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `SQMG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Quinovose([0x06 => Sulfo()]; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x03)]), 
        (Radyl, )
    )
end
function class_struct_SQDG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `SQDG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Quinovose([0x06 => Sulfo()]; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x03)]), 
        (Radyl, Radyl)
    )
end
function class_struct_MGMG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `MGMG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Galactose(; DL = DForm), gly; sil, linkage = [β(0x01) => lk(0x03)]), 
        (Radyl, )
    )
end
function class_struct_MGDG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `MGDG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Galactose(; DL = DForm), gly; sil, linkage = [β(0x01) => lk(0x03)]), 
        (Radyl, Radyl)
    )
end
function class_struct_DGMG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `DGMG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Galactose(; DL = DForm), Galactose(; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x06), β(0x01) => lk(0x03)]), 
        (Radyl, )
    )
end
function class_struct_DGDG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `DGDG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, Galactose(; DL = DForm), Galactose(; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x06), β(0x01) => lk(0x03)]), 
        (Radyl, Radyl)
    )
end
function class_struct_GlcAMG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `GlcAMG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, GlucuronicAcid(; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x03)]), 
        (Radyl, )
    )
end
function class_struct_GlcADG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `GlcADG` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (Omodifiedradylglycerol, 
        makechemical(DehydratedChemical, GlucuronicAcid(; DL = DForm), gly; sil, linkage = [α(0x01) => lk(0x03)]), (Radyl, Radyl))
end

function class_struct_PA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PA` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPA` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_PC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PC` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Choline(), PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPC` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Choline(), PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_PE(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PE` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPE(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPE` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_PE_NMe(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PE-NMe` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, Nmethylethanolamine(), PhosphoricAcid(), gly; sil),
        (Radyl, Radyl)
    )
end
function class_struct_LPE_NMe(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPE-NMe` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, Nmethylethanolamine(), PhosphoricAcid(), gly; sil), 
        (Radyl, )
    )
end
function class_struct_PE_NMe2(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PE-NMe2` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, NNdimethylethanolamine(), PhosphoricAcid(), gly; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_LPE_NMe2(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPE-NMe2` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, NNdimethylethanolamine(), PhosphoricAcid(), gly; sil), 
        (Radyl, )
    )
end
function class_struct_PE_N(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PE-N` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    post = replace(post, r"^\(" => "", r"\)$" => "")
    if post == "FA"
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FA 0:0"), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl, Acyl)
        )
    elseif post == "Alk"
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FOH 0:0"), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl, Alkyl)
        )
    else
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup(post), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl)
        )
    end
end
function class_struct_LPE_N(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPE-N` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    post = replace(post, r"^\(" => "", r"\)$" => "")
    if post == "FA"
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FA 0:0"), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, Acyl)
        )
    elseif post == "Alk"
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FOH 0:0"), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, Alkyl)
        )
    else
        (Radylglycerophosphate, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup(post), 
                makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil)), 
            (Radyl, )
        )
    end
end
function class_struct_PS(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PS` is not suppoted, please use `GP` instead"))
    ser, gly = parse_aa_glycerol_rs(pre)
    (Radylglycerophosphoaminoacid, makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPS(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPS` is not suppoted, please use `LGP` instead"))
    ser, gly = parse_aa_glycerol_rs(pre)
    (Radylglycerophosphoaminoacid, makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_PS_N(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PS-N` is not suppoted, please use `GP` instead"))
    ser, gly = parse_aa_glycerol_rs(pre)
    post = replace(post, r"^\(" => "", r"\)$" => "")
    if post == "FA"
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FA 0:0"), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl, Acyl)
        )
    elseif post == "Alk"
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FOH 0:0"), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl, Alkyl)
        )
    else
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup(post), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, Radyl)
        )
    end
end
function class_struct_LPS_N(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPS-N` is not suppoted, please use `LGP` instead"))
    ser, gly = parse_aa_glycerol_rs(pre)
    post = replace(post, r"^\(" => "", r"\)$" => "")
    if post == "FA"
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FA 0:0"), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, Acyl)
        )
    elseif post == "Alk"
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup("FOH 0:0"), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, Alkyl)
        )
    else
        (Radylglycerophosphoaminoacid, 
            concatchemical(DehydratedChemical, 
                parse_tailgroup(post), 
                makechemical(DehydratedChemical, ser, PhosphoricAcid(), gly; sil)), 
            (Radyl, )
        )
    end
end
function class_struct_PI(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PI` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, Inositol(), PhosphoricAcid(), gly; sil, linkage = [α(0x01) => lk(nothing), lk(nothing) => lk(0x01)]), 
        (Radyl, Radyl)
    )
end
function class_struct_LPI(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPI` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, Inositol(), PhosphoricAcid(), gly; sil, linkage = [α(0x01) => lk(nothing), lk(nothing) => lk(0x01)]), 
        (Radyl, )
    )
end
function class_struct_PMeOH(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PMeOH` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Methanol(), PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPMeOH(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPMeOH` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Methanol(), PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_PEtOH(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PEtOH` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Ethanol(), PhosphoricAcid(), gly; sil), (Radyl, Radyl))
end
function class_struct_LPEtOH(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPEtOH` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, makechemical(DehydratedChemical, Ethanol(), PhosphoricAcid(), gly; sil), (Radyl, ))
end
function class_struct_GP(head, pre, cls, post, pos, sil)
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
    concatchemical(DehydratedChemical, parse_headgroup(head), PhosphoricAcid(), gly; sil), # makechemical => concatchemical?
        (Radyl, Radyl)
    )
end
function class_struct_LGP(head, pre, cls, post, pos, sil)
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
    concatchemical(DehydratedChemical, parse_headgroup(head), PhosphoricAcid(), gly; sil), # makechemical => concatchemical?
        (Radyl, Radyl)
    )
end
function class_struct_PIP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PIP` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(1)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_LPIP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPIP` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(1)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_PIP2(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PIP2` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(2)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin + 1]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_LPIP2(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPIP2` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(2)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin + 1]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, )
    )
end
function class_struct_PIP3(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PIP3` is not suppoted, please use `GP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(3)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin + 2]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_LPIP3(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPIP3` is not suppoted, please use `LGP` instead"))
    gly = parse_glycerol_rs(pre)
    (Radylglycerophosphate, 
        makechemical(DehydratedChemical, 
            isnothing(pos) ? Inositol([Phosphoryl() => UInt8(3)]) : 
                Inositol([parse(UInt8, x.match) => Phosphoryl() for x in eachmatch(r"\d+", pos)][begin:begin + 2]), 
            PhosphoricAcid(), 
            gly; sil), 
        (Radyl, )
    )
end
function class_struct_PG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PG` is not suppoted, please use `GP` instead"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Radylglycerophosphoglycerol, makechemical(DehydratedChemical, gly2, PhosphoricAcid(), gly1; sil), (Radyl, Radyl))
end
function class_struct_LPG(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPG` is not suppoted, please use `LGP` instead"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Radylglycerophosphoglycerol, makechemical(DehydratedChemical, gly2, PhosphoricAcid(), gly1; sil), (Radyl, ))
end
function class_struct_PGP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `PGP` is not suppoted, please use `GP` instead"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Radylglycerophosphoglycerol, 
        makechemical(DehydratedChemical, PhosphoricAcid(), gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_LPGP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LPGP` is not suppoted, please use `LGP` instead"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Radylglycerophosphoglycerol, 
        makechemical(DehydratedChemical, PhosphoricAcid(), gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, )
    )
end
function class_struct_BPA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `BPA` is not suppoted"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Bisradylglycerophosphate, 
        makechemical(DehydratedChemical, gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl, Radyl, Radyl)
    )
end
function class_struct_SLBPA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `SLBPA` is not suppoted"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Bisradylglycerophosphate, 
        makechemical(DehydratedChemical, gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl, Radyl)
    )
end
function class_struct_LBPA(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LBPA` is not suppoted"))
    gly2, gly1 = parse_glycerol2_rs(pre)
    (Bisradylglycerophosphate, 
        makechemical(DehydratedChemical, gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_CL(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `CL` is not suppoted"))
    gly3, gly2, gly1 = parse_glycerol3_rs(pre)
    (Bisradylglycerophosphoglycerol, 
        makechemical(DehydratedChemical, gly3, PhosphoricAcid(), gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl, Radyl, Radyl)
    )
end
function class_struct_MLCL(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `MLCL` is not suppoted"))
    gly3, gly2, gly1 = parse_glycerol3_rs(pre)
    (Bisradylglycerophosphoglycerol, 
        makechemical(DehydratedChemical, gly3, PhosphoricAcid(), gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl, Radyl)
    )
end
function class_struct_DLCL(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `DLCL` is not suppoted"))
    gly3, gly2, gly1 = parse_glycerol3_rs(pre)
    (Bisradylglycerophosphoglycerol, 
        makechemical(DehydratedChemical, gly3, PhosphoricAcid(), gly2, PhosphoricAcid(), gly1; sil), 
        (Radyl, Radyl)
    )
end
function class_struct_GP_NAE(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `GP-NAE` is not suppoted"))
    gly = parse_glycerol_rs(pre)
    (GlycerophosphoNacylethanolamine, 
        makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), gly; sil), 
        (Acyl, )
    )
end

function class_struct_Cer(head, pre, cls, post, pos, sil)
    if isnothing(head)
        isnothing(sil) || @warn "Ignore stable isotope labeling on class"
        return (SphingoBone, nothing, (SPB, Acyl))
    end
    head = try
        parse_headgroup(head)
    catch
        parse_glycomp(head)
    end
    if head isa Tuple
        (MixedSphingoBone, 
            getchaincomponent(makechemical(DehydratedChemical, head...; sil)), # distribute sil directly?
            (SPB, Acyl)
        )
    elseif head isa GlyComp
        isnothing(pos) || @warn "Ignore headgroup position for unstructured glycan"
        (SphingoBone, head, (SPB, Acyl))
    else
        (SphingoBone, concatchemical(typeof(head), head; sil), (SPB, Acyl))
    end
end
function class_struct_SPB(head, pre, cls, post, pos, sil)
    if isnothing(head)
        isnothing(sil) || @warn "Ignore stable isotope labeling on class"
        return (SphingoBone, nothing, (SPB, ))
    end
    head = try
        parse_headgroup(head)
    catch
        parse_glycomp(head)
    end
    if head isa Tuple
        (MixedSphingoBone, 
            getchaincomponent(makechemical(DehydratedChemical, head...; sil)), 
            (SPB, )
        )
    elseif head isa GlyComp
        isnothing(pos) || @warn "Ignore headgroup position for unstructured glycan"
        (SphingoBone, head, (SPB, Acyl))
    else
        (SphingoBone, concatchemical(typeof(head), head; sil), (SPB, ))
    end
end
function class_struct_CerP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `CerP` is not suppoted, please use `Cer` instead"))
    (SphingoBone, makechemical(PhosphoricAcid(); sil), (SPB, Acyl))
end
function class_struct_SPBP(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `SPBP` is not suppoted, please use `SPB` instead"))
    (SphingoBone, makechemical(PhosphoricAcid(); sil), (SPB, ))
end
function class_struct_EPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `EPC` is not suppoted, please use `Cer` instead"))
    (SphingoBone, makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), (SPB, Acyl))
end
function class_struct_LEPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LEPC` is not suppoted, please use `SPB` instead"))
    (SphingoBone, makechemical(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), (SPB, ))
end
function class_struct_IPC(head, pre, cls, post, pos, sil)
    isnothing(head) && return (SphingoBone, 
        makechemical(DehydratedChemical, Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(nothing)]), 
        (SPB, Acyl)
    ) 
    head = parse_headgroup(head)
    head isa Tuple && throw(ArgumentError("Multiple headgroups on SPB are not suppoted for subclass `IPC`, please use `Cer` instead"))
    (SphingoBone, 
        makechemical(DehydratedChemical, head, Inositol(), PhosphoricAcid(); sil, linkage = [last(getchainlinkage(head)), α(0x01) => lk(nothing)]), 
        (SPB, Acyl)
    )
end
function class_struct_LIPC(head, pre, cls, post, pos, sil)
    isnothing(head) && return (SphingoBone, 
        makechemical(DehydratedChemical, Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(nothing)]), 
        (SPB, )
    ) 
    head = parse_headgroup(head)
    head isa Tuple && throw(ArgumentError("Multiple headgroups on SPB are not suppoted for subclass `LIPC`, please use `SPB` instead"))
    (SphingoBone, 
        makechemical(DehydratedChemical, head, Inositol(), PhosphoricAcid(); sil, linkage = [last(getchainlinkage(head)), α(0x01) => lk(nothing)]), 
        (SPB, )
    )
end
function class_struct_MIPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `MIPC` is not suppoted, please use `IPC` instead"))
    (SphingoBone, 
        makechemical(DehydratedChemical, Mannose(), Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x02), α(0x01) => lk(nothing)]), 
        (SPB, Acyl)
    )
end
function class_struct_LMIPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LMIPC` is not suppoted, please use `LIPC` instead"))
    (SphingoBone, 
        makechemical(DehydratedChemical, Mannose(), Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x02), α(0x01) => lk(nothing)]), 
        (SPB, )
    )
end
function class_struct_MIPIPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `M(IP)2C` is not suppoted, please use `IPC` instead"))
    (SphingoBone, 
            makechemical(DehydratedChemical, 
            Inositol(), 
            PhosphoricAcid(), 
            Mannose(), 
            Inositol(), 
            PhosphoricAcid(); sil, linkage = [α(0x01) => lk(nothing), lk(nothing) => lk(0x06), α(0x01) => lk(0x02), α(0x01) => lk(nothing)]), 
        (SPB, Acyl)
    )
end
function class_struct_LMIPIPC(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LM(IP)2C` is not suppoted, please use `LIPC` instead"))
    (SphingoBone, 
            makechemical(DehydratedChemical, 
            Inositol(), 
            PhosphoricAcid(), 
            Mannose(), 
            Inositol(), 
            PhosphoricAcid(); sil, linkage = [α(0x01) => lk(nothing), lk(nothing) => lk(0x06), α(0x01) => lk(0x02), α(0x01) => lk(nothing)]), 
        (SPB, )
    )
end
function class_struct_SM(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `SM` is not suppoted, please use `Cer` instead"))
    (SphingoBone, makechemical(DehydratedChemical, Choline(), PhosphoricAcid(); sil), (SPB, Acyl))
end
function class_struct_LSM(head, pre, cls, post, pos, sil)
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `LSM` is not suppoted, please use `SPB` instead"))
    (SphingoBone, makechemical(DehydratedChemical, Choline(), PhosphoricAcid(); sil), (SPB, ))
end
function class_struct_GSL(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    (SphingoBone, makechemical(parse_gsl_isomer(cls, post); sil), (SPB, Acyl))
end
function class_struct_LGSL(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `SPB` instead"))
    (SphingoBone, makechemical(parse_gsl_isomer(cls, post); sil), (SPB, ))
end
function class_struct_GlcCer(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    class_struct_Cer("D-Glcβ-", pre, "Cer", post, pos, sil)
end
function class_struct_GlcSPB(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `SPB` instead"))
    class_struct_Cer("D-Glcβ-", pre, "SPB", post, pos, sil)
end
function class_struct_GalCer(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    class_struct_Cer("D-Galβ-", pre, "Cer", post, pos, sil)
end
function class_struct_GalSPB(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `SPB` instead"))
    class_struct_Cer("D-Galβ-", pre, "SPB", post, pos, sil)
end
function class_struct_LacCer(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    class_struct_Cer("D-Galβ-4D-Glcβ-", pre, "Cer", post, pos, sil)
end
function class_struct_LacSPB(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `SPB` instead"))
    class_struct_Cer("D-Galβ-4Glcβ-", pre, "SPB", post, pos, sil)
end
function class_struct_SL(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    (SphingoBone, nothing, (SulfoSPB, Acyl))
end
function class_struct_LSL(head, pre, cls, post, pos, sil) 
    isnothing(head) || throw(ArgumentError("Headgroup extension on subclass `$cls` is not suppoted, please use `Cer` instead"))
    (SphingoBone, nothing, (SulfoSPB, ))
end

function class_struct_ACer(head, pre, cls, post, pos, sil)
    isnothing(sil) || @warn "Ignore stable isotope labeling on class"
    (SphingoBone, parse_lipid(pre), (SPB, Acyl, Acyl))
end
class_struct_ASM(head, pre, cls, post, pos, sil) = (MixSphingoBone, (parse_lipid(pre), makechemical(DehydratedChemical, Choline(), PhosphoricAcid(); sil)), (SPB, Acyl, Acyl))

class_struct_ST(head, pre, cls, post, pos, sil) = (SterolBone, nothing, (STRing, ))
class_struct_SE(head, pre, cls, post, pos, sil) = (SterolBone, nothing, (STRing, Acyl))
class_struct_SG(head, pre, cls, post, pos, sil) = (SubstitutedSterol, nothing, (STRing, ))
class_struct_ASG(head, pre, cls, post, pos, sil) = (SubstitutedSterol, nothing, (STRing, ))
class_struct_FC(head, pre, cls, post, pos, sil) = (SterolBone, nothing, (CRing, ))
class_struct_SST(head, pre, cls, post, pos, sil) = (SterolBone, nothing, (eval(string(cls, "Ring")), ))
class_struct_SSE(head, pre, cls, post, pos, sil) = (SterolBone, nothing, (eval(string(replace(cls, r"E$" => ""), "Ring"), Acyl)))
class_struct_SSG(head, pre, cls, post, pos, sil) = (SubstitutedSterol, nothing, (eval(string(replace(cls, r"G$" => ""), "Ring")), ))
class_struct_SASG(head, pre, cls, post, pos, sil) = (SubstitutedSterol, nothing, (eval(string(replace(cls, r"G$" => "", r"^A" => ""), "Ring")), ))
class_struct_BA(head, pre, cls, post, pos, sil) = (SubstitutedSterol, nothing, (BARing, ))