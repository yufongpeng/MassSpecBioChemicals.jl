using MassSpecChemicals
using MassSpecBioChemicals
using MassSpecBioChemicals.Lipids
const MBL = MassSpecBioChemicals.Lipids
using Pkg
Pkg.activate("test")
using JSON3
Pkg.activate(".")

lpc = "LPC 36:1;O2"
lpc = "LPC 18:1(2);OH"
s = "DG 16:0/18:1/0:0"
s = "PC 18:1(2E);3OH/18:1(1)"
s = "(S)-PC 18:1(2E);3OH[R],4OH[S]/18:1(16);17Me"
s = "(S)-PC 18:1(2E);3OH[R],4OH[S]/18:1(16);17Me;10Ep[R]"
s = "(S)-PC 18:1(2E);3OH[R],4OH[S]/18:1(16);17Me;10Ep[R];11H[R]"
s = "MIPC(1) 18:1(4E);3OH/18:0"
s = "GM3(1) 18:1(4E);3OH/18:0"
s = "FAHFA 18:1(4E)/5O(FA 18:0)"
s = "WE 36:1"

s = "FA 18:0;18COOH"
s = "FOH 18:1(16);17Me"
s = "FA 19:0;19COOH;10OH"
s = "HC 19:2(1,18);10OH[R]"
s = "FOH 19:2(3Z,16Z);10OH[R],19OH"
s = "FN 19:2(3Z,16Z);19NH2;10OH"
s = "FAM 19:2(3Z,16Z);19CONH2;9OH,10OH,11OH" # Full
s = "FAM 19:2(3Z,16Z);19CONH2;9OH[R],10OH,11OH[R]"
s = "FAM 19:2(3Z,16Z);19CONH2;9OH[R],10OH[R],11OH[R]"
s = "FAL 18:0;9oxo,10oxo,18oxo"
s = "L-NASer 16:0"
s = "L-NASer 19:2(3Z,16Z);19oxo;19COSer;9OH[R],10OH[R],11OH[R]" # Full
s = "L-NASer 19:2(3Z,16Z);19oxo;19CO(L-Ser);9OH[R],10OH[R],11OH[R]"
s = "L-NASer 19:2(3Z,16Z);19CO(L-Ser-L-Ala);9OH[R],10OH[R],11OH[R]"

s = "FAM 19:2(3Z,15Z);18oxo;18NH2;10OH" # X

s = "FA 24:1-ASM 18:1(4E);3OH/24:1(15Z)"
s = "L-NASer 19:2(3Z,16Z);19CO(L-Ser);9OH[R],10OH[R],11OH[R]"
s = "Glcβ-SPB(1) 18:1(4E);3OH[S];2H[R]"
s = "DL-Glcβ-SPB(1) 18:1(4E);3OH[S];2H[R]"
s = "HexCer 18:1;O2/24:1"
s = "Hex-Hex-Cer 18:1;O2/24:1"
s = "D-Glcβ-DG 16:0/20:4(5Z,8Z,11Z,14Z)"
s = "SPBP(1,3R) 18:1(4E);2H[R]"
s = "SPBP(1,3R) 18:1(4E)"
s = "D-Gal2D-GlcD-GlcNAcCer 18:1(4E);3OH/24:1(15Z)"
s = "TG 18:1(9Z)/18:1(16);17Me/20:0"
s = "TG 18:1(16);17Me_18:1(9Z)_20:0(sn-1)"
s = "TG 18:2(5,16);17Me_18:1(16Z)_20:0(sn-1)"
s = "TG 18:2(5,16);17Me_18:2(5,16);5Me_20:0(sn-1)"
s = "TG 18:2(5,16);10Me_18:2(5,16);10Me[R]_20:0(sn-1)"
s = "TG 18:2(5,16);9OH;10Me_18:2(5,16);10Me[R],11Me_20:0(sn-1)"
s = "TG 18:2(5,16);9OH;10Me_18:2(5,16);Me2_20:0(sn-1)"
s = "BPA 16:0(sn-2)_20:0(sn-3')_40:4"

lipid = parse_chemical(LipidParser(), s)
lipid2 = parse_chemical(LipidParser(), s)
lipid = parse_chemical(LipidParser(; explicit_config = true), s)
parse_chemical(LipidParser(), s) == parse_chemical(LipidParser(), s)

MBL.annotationlevel(lipid)
chemicalsmiles(lipid; onlycarbonchain = true)

test_lipid_js = JSON3.read(joinpath("test", "data", "test_lipid.json"))
test_lipid = Dict{UnionAll, Dict}()
for (c, s) in test_lipid_js
    dict = Dict{String, Any}()
    for (l, a) in s
        l = string(l)
        print(l, " -- ")
        lipid = parse_chemical(LipidParser(), l)
        show(lipid)
        first(a) == chemicalname(lipid) || println("\n---------", first(a))
        println()
        push!(dict, l => (object = lipid, name = first(a), annotationlevel = [eval(Meta.parse(aa)) for aa in last(a)]))
    end
    push!(test_lipid, eval(c) => dict)
end

function test_annotationlevel(x)
    tal = MBL.annotationlevel(x.object; partial = true, additional = true, pass = true)
    length(tal) == length(x.annotationlevel) && all(t -> in(t, x.annotationlevel), tal)
end

all(test_annotationlevel(x) for (k, x) in test_lipid[FattyAcyl])
function show_error(L)
    i = findall(!test_annotationlevel(x) for (k, x) in test_lipid[L])
    for j in i
        println(j, " -- ")
        println(" ---- ", MBL.annotationlevel(test_lipid[L][j].object; partial = true, additional = true, pass = true))
        println(" -------- ", test_lipid[L][j].annotationlevel)
    end
end
show_error(FattyAcyl)
all(test_annotationlevel(x) for (k, x) in test_lipid[Glycerolipid])
show_error(Glycerolipid)
all(test_annotationlevel(x) for (k, x) in test_lipid[Glycerophospholipid])
show_error(Glycerophospholipid)
all(test_annotationlevel(x) for (k, x) in test_lipid[Sphingolipid])
show_error(Sphingolipid)

chemicalsmiles(MBL.parse_lipid("DG 18:1(9Z)/18:1(16);17Me/0:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("DG P-18:1(9Z);11oxo;12Ep;13OO/18:1(2Z);2OH;10OMe;18COOH/0:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("DG O-18:1(9Z);2oxy;11oxo;12Ep;13H[R]/18:0;2OH;10OFo;18COCl/0:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("Cer 18:1(4E);1OH,3OH;2H[R]/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("Glc-Cer(1) 18:1(4E);3OH;2H[R]/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("FA 10:0-ACer(1) 18:1(4E);3OH;2H[R]/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("SM(1) 18:1(4E);3OH;2H[R]/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("PC 18:1(9Z)/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("(S)-PG 18:1(9Z)/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("LPC 0:0/18:1(16);17Me"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("FA 18:1(9Z)"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("FA 18:1(9Z)"))
chemicalsmiles(MBL.parse_lipid("FN 18:1(9Z)"))
chemicalsmiles(MBL.parse_lipid("FAM 18:1(9Z)"))
chemicalsmiles(MBL.parse_lipid("FOH 18:1(1Z)"))
chemicalsmiles(MBL.parse_lipid("FAHFA 18:1(4E)/5O(FA 18:0)"))
chemicalsmiles(MBL.parse_lipid("WE 18:1(1Z)/18:0"))
chemicalsmiles(MBL.parse_lipid("NA 18:1(1Z)/18:0"))
chemicalsmiles(MBL.parse_lipid("NASer 18:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("NAE 18:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("CAR 18:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("GP-NAE 18:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("(R,S')-LBPA 18:0/0:0/16:0/0:0"); onlycarbonchain = true)
chemicalsmiles(MBL.parse_lipid("(R,S',R'')-MLCL 16:0/20:4(5Z,8Z,11Z,14Z)/18:0/0:0"); onlycarbonchain = true)

for (k, v) in test_lipid[FattyAcyl]
    print(k, " -- ")
    println(chemicalsmiles(v.object; onlycarbonchain = true))
end

for (k, v) in test_lipid[Glycerolipid]
    print(k, " -- ")
    println(chemicalsmiles(v.object; onlycarbonchain = true))
end

for (k, v) in test_lipid[Glycerophospholipid]
    print(k, " -- ")
    println(chemicalsmiles(v.object; onlycarbonchain = true))
end

for (k, v) in test_lipid[Sphingolipid]
    print(k, " -- ")
    println(chemicalsmiles(v.object; onlycarbonchain = true))
end

for (k, v) in test_lipid[FattyAcyl]
    print(k, " -- ")
    println(MBL.dissociate_carbonchain_group(v.object))
end

for (k, v) in test_lipid[Glycerolipid]
    print(k, " -- ")
    println(MBL.dissociate_carbonchain_group(v.object))
end

for (k, v) in test_lipid[Glycerophospholipid]
    print(k, " -- ")
    println(MBL.dissociate_carbonchain_group(v.object))
end

for (k, v) in test_lipid[Sphingolipid]
    print(k, " -- ")
    println(MBL.dissociate_carbonchain_group(v.object))
end