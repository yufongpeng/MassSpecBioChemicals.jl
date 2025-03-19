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
s = "LPC 18:1(2E);3OH/0:0"
s = "MIPC(1) 18:1(4E);3OH/18:0"
s = "GM3(1) 18:1(4E);3OH/18:0"

lipid = MBL.parse_lipid(s)
chemicalsmiles(lipid; onlycarbonchain = true)

test_lipid_js = JSON3.read(joinpath("test", "data", "test_lipid.json"))
test_lipid = Dict{UnionAll, Dict}()
for (c, s) in test_lipid_js
    dict = Dict{String, Any}()
    for (l, a) in s
        l = string(l)
        lipid = MBL.parse_lipid(l)
        show(lipid)
        print(", ")
        push!(dict, l => (object = lipid, annotationlevel = [eval(Meta.parse(aa)) for aa in a]))
    end
    push!(test_lipid, eval(c) => dict)
end

function test_annotationlevel(x)
    tal = MBL.annotationlevel(x.object; partial = true, additional = true, pass = true)
    length(tal) == length(x.annotationlevel) && all(t -> in(t, x.annotationlevel), tal)
end

all(test_annotationlevel(x) for (k, x) in test_lipid[FattyAcyl])
all(test_annotationlevel(x) for (k, x) in test_lipid[Glycerolipid])
all(test_annotationlevel(x) for (k, x) in test_lipid[Glycerophospholipid])
all(test_annotationlevel(x) for (k, x) in test_lipid[Sphingolipid])
findall(!test_annotationlevel(x) for (k, x) in test_lipid[Sphingolipid])