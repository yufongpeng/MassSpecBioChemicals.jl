# MassSpecBioChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecBioChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=test)](https://github.com/yufongpeng/MassSpecBioChemicals.jl/actions/workflows/CI.yml?query=branch%3Atest)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecBioChemicals.jl/branch/test/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecBioChemicals.jl)

MassSpecBioChemicals.jl is package for representation of biochemicals with four submodules, `Lipids`, `Glycans`, `Proteins`, and `Nucleotides`. It is built on [MassSpecChemicals.jl](https://github.com/yufongpeng/MassSpecChemicals.jl); however, it is still in early-stage, so interfaces are not fully implemented. Currently, the development is focused on submodule `Lipids`, especially the tranformation of standard lipid nomenclatures to internal lipid objects, generation of SMILES, and manipulation of internal lipid objects. The following features to be implemented is the transformation of lipids at less exact annotation level to full structure level to generate SMILES. The transformation will be based on either public lipid libraries filtered by taxonomy Information, or user-defined fatty acid composition. 

## Lipids
### Input
Standard lipid nomenclatures are used in this package; other lipids formats are currently not supported[1]. Nomenclatures for glycans are mainly based on the Symbol Nomenclature For Glycans (SNFG) guidline[2].
```julia
julia> using MassSpecBioChemicals
       using MassSpecBioChemicals.Lipids
       const MBL = MassSpecBioChemicals.Lipids
MassSpecBioChemicals.Lipids

julia> fahfa = parse_chemical(LipidParser(), "FAHFA 18:1(4E)/5O(FA 18:0)")
FAHFA 18:1(4E)/5O(FA 18:0)

julia> naser = parse_chemical(LipidParser(), "L-NASer 19:2(3Z,16Z);19CO(L-Ser-L-Ala);9OH[R],10OH[R],11OH[R]")
L-NASer 19:2(3Z,16Z);19CO(L-Ser-L-Ala);9OH[R],10OH,11OH[R]

julia> fam = parse_chemical(LipidParser(), "FAM 19:2(3Z,16Z);19CONH2;9OH,10OH,11OH")
FAM 19:2(3Z,16Z);19CONH2;9OH,10OH,11OH

julia> pc = parse_chemical(LipidParser(), "PC 18:1(2E);3OH/18:1(15)")
PC 18:1(2E);3OH/18:1(15)

julia> pip2 = parse_chemical(LipidParser(), "PIP2(4',5') 20:4(5Z,8Z,11Z,14Z)/16:0")
PIP2(4',5') 20:4(5Z,8Z,11Z,14Z)/16:0

julia> cl = parse_chemical(LipidParser(), "CL 16:0(sn-1)_20:0(sn-2')_40:4")
CL 40:4_16:0(sn-1)_20:0(sn-2')

julia> pen = parse_chemical(LipidParser(), "PE-N(FA 18:1(12Z)) 20:4(5Z,8Z,11Z,14Z)/16:0")
PE-N(FA 18:1(12Z)) 20:4(5Z,8Z,11Z,14Z)/16:0

julia> gsl = parse_chemical(LipidParser(), "Galβ-3GalNAcβ-4[Neu5Acα-3]D-Galβ-4Glcβ-Cer(1) 18:1(4E);3OH[S];2H[R]/24:1(15Z)")
D-Galβ1-3D-GalNAcβ1-4[Neu5Acα1-3]D-Galβ1-4D-Glcβ1-Cer(1) 18:1(4E);2H[R];3OH[S]/24:1(15Z)

julia> gq1 = parse_chemical(LipidParser(), "GQ1bα(1) 18:1(4E);3OH[S];2H[R]/24:1(15Z)")
GQ1bα(1) 18:1(4E);2H[R];3OH[S]/24:1(15Z)

julia> acer = parse_chemical(LipidParser(), "FA 24:1(15)-ACer(1) 18:1(4E);3OH[R];2H[R]/24:1(15Z)")
FA 24:1(15)-ACer(1) 18:1(4E);2H[R];3OH[R]/24:1(15Z)

```
The function checks the number and position of functional groups and double bond, as well as chirality. 
```julia
julia> parse_chemical(LipidParser(), "LPC 18:1(1)")
ERROR: ArgumentError: Position `1` cannot contain double bond.
…

julia> parse_chemical(LipidParser(), "LPC 18:1(2);18oxo;18COOH")
ERROR: ArgumentError: Too many functional groups at position `18`.
…

julia> parse_chemical(LipidParser(), "LPC 18:2(3,3)")
ERROR: ArgumentError: Multiple doublebond at position `3`.
…

```
### Class
|Category|Class|Pre class modifier|Post class modifier|Description|Example|
|--------|-----|------------------|-------------------|-----------|-------|
|Fatty acyl|FA|||Fatty acid|FA 18:0|
|Fatty acyl|HC|||Hydrocarbon|HC 16:0|
|Fatty acyl|FAL|||Fatty aldehyde|FAL 18:0|
|Fatty acyl|FOH|||Fatty alcohol|FOH 16:0|
|Fatty acyl|FN|||Fatty amine|FN 18:0|
|Fatty acyl|FAM|||Fatty amide|FAM 16:0|
|Fatty acyl|CAR|L, D, (R), (S)||Acylcarnitine|(R)-CAR 18:0|
|Fatty acyl|CoA|||Acyl Coenzyme A|CoA 16:0|
|Fatty acyl|NA[AA]|L, D, (R), (S)||Nacylaminoacid|L-NASer 16:0|
|Fatty acyl|NAE|||Nacylethanolamine|NAE 18:0|
|Fatty acyl|NAT|||Nacyltaurine|NAT 16:0|
|Fatty acyl|WE|||Wax ester|WE 16:0/18:0|
|Fatty acyl|NA|||Nacylalkylamine|NA 16:0/18:0|
|Fatty acyl|FAHFA|||Fatty acyl estolid|FAHFA 20:4/7O(FA 16:0)|
|Glycerolipid|MG|Glycan-Glycerol chirality||Monoradylglycerol and its derivative|MG 16:0/0:0/0:0|
|Glycerolipid|DG|Glycan-Glycerol chirality||Diradylglycerol and its derivative|Glcβ-(R)-DG 16:0_18:0|
|Glycerolipid|TG|Glycerol chirality||Triradylglycerol|TG 16:0_18:0_20:1(sn-1)|
|Glycerolipid|SQMG|Glycerol chirality||Sulfoquinovosylmonoradylglycerol|SQMG 16:0|
|Glycerolipid|SQDG|Glycerol chirality||Sulfoquinovosyldiradylglycerol|SQDG 16:0/18:0|
|Glycerolipid|MGMG|Glycerol chirality||Monogalactosylmonoradylglycerol|MGMG 16:0|
|Glycerolipid|MGDG|Glycerol chirality||Monogalactosyldiradylglycerol|MGDG 16:0/18:0|
|Glycerolipid|DGMG|Glycerol chirality||Digalactosylmonoradylglycerol|DGMG 16:0|
|Glycerolipid|DGDG|Glycerol chirality||Digalactosyldiradylglycerol|DGDG 16:0/18:0|
|Glycerolipid|GlcAMG|Glycerol chirality||Glucuronosylmonoradylglycerol|GlcAMG 16:0|
|Glycerolipid|GlcADG|Glycerol chirality||Glucuronosyldiradylglycerol|GlcADG 16:0/18:0|
|Glycerophospholipid|GP|Headgroup-Glycerol chirality||Glycerophospholipid|D-Glcβ-(R)-GP 18:0/16:0|
|Glycerophospholipid|PA|Glycerol chirality||Phosphatidic acid|PA 16:0/18:0|
|Glycerophospholipid|PC|Glycerol chirality||Phosphatidylcholine|(R)-PC O-16:0/18:0|
|Glycerophospholipid|PE|Glycerol chirality||Phosphatidylethanolamine|PE P-16:0/18:0|
|Glycerophospholipid|PS|Glycerol chirality-Serine chirality||Phosphatidylserine|(S)-L-PS 16:0/18:0|
|Glycerophospholipid|PI|Glycerol chirality||Phosphatidylinositol|PI 16:0_18:0|
|Glycerophospholipid|PG|Glycerol chirality||Phosphatidylglycerol|(R,R')-PG 16:0/18:0|
|Glycerophospholipid|PGP|Glycerol chirality||Phosphatidylglycerol phophate|PGP 16:0/18:0|
|Glycerophospholipid|PMeOH|Glycerol chirality||Phosphatidylmethanol|PMeOH 16:0_18:0|
|Glycerophospholipid|PEtOH|Glycerol chirality||Phosphatidylethanol|PEtOH 16:0_18:0|
|Glycerophospholipid|PIP|Glycerol chirality|Phosphoryl position|Phosphatidylinositol phosphate|PIP(3') 20:4;O_16:0|
|Glycerophospholipid|PIP2|Glycerol chirality|Phosphoryl position|Phosphatidylinositol diphosphate|PIP2(3',4') 20:4;O_16:0|
|Glycerophospholipid|PIP3|Glycerol chirality|Phosphoryl position|Phosphatidylinositol triphosphate|PIP3(3',4', 5') 20:4;O_16:0|
|Glycerophospholipid|PE-N|Glycerol chirality|Fatty acyl|N-modified phosphatidylethanolamine|PE-N(FA 18:0) 16:0/20:4|
|Glycerophospholipid|PE-NMe|Glycerol chirality||N-methyl phosphatidylethanolamine|PE-NMe 16:0/18:0|
|Glycerophospholipid|PE-NMe2|Glycerol chirality||N,N-dimethyl phosphatidylethanolamine|PE-NMe2 16:0/18:0|
|Glycerophospholipid|PS-N|Glycerol chirality-Serine chirality|Fatty acyl|N-modified phosphatidylserine|PS-N(18:0) 16:0/18:0|
|Glycerophospholipid|PS-NMe|Glycerol chirality-Serine chirality||N-methyl phosphatidylserine|PS-N(FA 16:0) 18:0/16:0|
|Glycerophospholipid|PS-NMe2|Glycerol chirality-Serine chirality||N,N-dimethyl phosphatidylserine|PS-N(FA 16:0) 18:0/16:0|PS-N(FA 16:0) 18:0/22:2|
|Glycerophospholipid|L[GP]|Glycerol chirality-Serine chirality||Lysopholipid|(S)-L-LPS-NMe2 16:0/0:0|
|Glycerophospholipid|BPA|Glycerol chirality||Bisphosphatidic acid|(R,R')-BPA 16:0/18:0/18:0/16:0|
|Glycerophospholipid|SLBPA|Glycerol chirality||Semilysobisphosphatidic acid|(R,R')-SLBPA 16:0/18:0/0:0/16:0|
|Glycerophospholipid|LBPA|Glycerol chirality||Lysobisphosphatidic acid|(R,R')-BPA 16:0/0:0/0:0/16:0|
|Glycerophospholipid|CL|Glycerol chirality||Cardiolipin|(R,S',R'')-CL 16:0/18:0/18:0/16:0|
|Glycerophospholipid|MLCL|Glycerol chirality||Monolysocardiolipin|(R,S',R'')-MLCL 16:0/18:0/18:0/0:0|
|Glycerophospholipid|DLCL|Glycerol chirality||Dilysocardiolipin|(R,S',R'')-DLCL 0:0/18:0/18:0/0:0|
|Glycerophospholipid|GP-NAE|Glycerol chirality||Glycerophosphaethanolamine|GP-NAE 16:0|
|Sphingolipid|SPB|Headgroup|Headgroup position|Sphingoid base and its derivative|PI-SPB(1) 18:1(4E);3OH|
|Sphingolipid|Cer|Headgroup|Headgroup position|Ceramide and its derivative|Glcβ-Cer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SPBP|Headgroup|Headgroup position|Sphingoid base- phosphate|SPBP(1,3R) 18:1(4E);2H[R]
|Sphingolipid|CerP|Headgroup|Headgroup position|Ceramide-phosphate|CerP(1) 18:1(4E);3OH/18:0|
|Sphingolipid|EPC||Headgroup position|Ethanolaminephosphorylceramide|EPC(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LEPC||Headgroup position|Ethanolaminephosphoryl-sphingoid base|LEPC(1) 18:1(4E);3OH|
|Sphingolipid|IPC||Headgroup position|Inositolphosphorylceramide|IPC(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LIPC||Headgroup position|Inositolphosphoryl-sphingoid base|LIPC(1) 18:1(4E);3OH|
|Sphingolipid|MIPC||Headgroup position|Mannosyl-inositolphosphorylceramide|MIPC(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LMIPC||Headgroup position|Mannosyl-inositolphosphoryl-sphingoid base|LMIPC(1) 18:1(4E);3OH|
|Sphingolipid|M(IP)2C||Headgroup position|Mannosyl-diinositolphosphorylceramide|M(IP)2C(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LM(IP)2C||Headgroup position|Mannosyl-diinositolphosphoryl-sphingoid base|LM(IP)2C(1) 18:1(4E);3OH|
|Sphingolipid|SM||Headgroup position|Sphingomyelin|SM(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LSM||Headgroup position|Lysosphingomyelin|LSM(1) 18:1(4E);3OH|
|Sphingolipid|SL||Headgroup position|Sulfonlipid|SL(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LSL||Headgroup position|Lysosulfonolipid|LSL(1) 18:1(4E);3OH|
|Sphingolipid|HexCer||Headgroup position|Hexosylceramide|HexCer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|HexSPB||Headgroup position|Hexosyl-sphingoid base|HexSPB(1) 18:1(4E);3OH|
|Sphingolipid|GlcCer||Headgroup position|Glucosylceramide|GlcCer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GlcSPB||Headgroup position|Glucosyl-sphingoid base|GlcSPB(1) 18:1(4E);3OH|
|Sphingolipid|GalCer||Headgroup position|Galactosylceramide|GalCer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GalSPB||Headgroup position|Galactosyl-sphingoid base|GalSPB(1) 18:1(4E);3OH|
|Sphingolipid|Hex2Cer||Headgroup position|Dihexosylceramide|Hex2Cer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Hex2SPB||Headgroup position|Dihexosyl-sphingoid base|Hex2SPB(1) 18:1(4E);3OH|
|Sphingolipid|LacCer||Headgroup position|Lactosylceramide|LacCer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LacSPB||Headgroup position|Lactosyl-sphingoid base|LacSPB(1) 18:1(4E);3OH|
|Sphingolipid|SHexCer||Headgroup position|Sulfatide|SHexCer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SHexSPB||Headgroup position|Lysosulfatide|SHexSPB(1) 18:1(4E);3OH|
|Sphingolipid|ACer|Fatty acyl|Headgroup position|Acylceramide|FA 18:0-ACer(1) 18:1(4E);3OH/18:0|
|Sphingolipid|ASM|Fatty acyl|Headgroup position|Acylsphingomyelin|FA 18:0-ASM(3R,1) 18:1(4E)/18:0|
|Sphingolipid|GM3||Headgroup position|Monosialoganglioside|GM3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GM2||Headgroup position|Monosialoganglioside|GM2(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GM1||Headgroup position|Monosialoganglioside|GM1a(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GD3||Headgroup position|Disialoganglioside|GD3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GD2||Headgroup position|Disialoganglioside|GD2(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GD1||Headgroup position|Disialoganglioside|GD1α(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GT3||Headgroup position|Trisialoganglioside|GT3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GT2||Headgroup position|Trisialoganglioside|GT2(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GT1||Headgroup position|Trisialoganglioside|GT1bα(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GQ1||Headgroup position|Tetrasialoganglioside|GQ1bα(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GP1||Headgroup position|Pentasialoganglioside|GP1c(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GA2||Headgroup position|Asialoganglioside|GA2(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GA1||Headgroup position|Asialoganglioside|GA1(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Gb3||Headgroup position|Globotriaosylceramide|Gb3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Gb4||Headgroup position|Globoside 4|Gb4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Gb5||Headgroup position|Globoside 5|Gb5(1) 18:1(4E);3OH/18:0|
|Sphingolipid|iGb3||Headgroup position|Isoglobotriaosylceramide|iGb3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|iGb4||Headgroup position|Isogloboside 4|iGb4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|iGb5||Headgroup position|Isogloboside 5|iGb5(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Lc3||Headgroup position|Lactotriaosylceramide|Lc3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|LM1||Headgroup position|Lactotetraosylceramide|LM1(1) 18:1(4E);3OH/18:0|
|Sphingolipid|Lc4||Headgroup position|Lactotetraosylceramide|Lc4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|nLc4||Headgroup position|Neolactotetraosylceramide|nLc4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|nLc5||Headgroup position|Neolactopentaosylceramide|nLc5(1) 18:1(4E);3OH/18:0|
|Sphingolipid|GM4||Headgroup position|Monosialoganglioside|GM4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SM4||Headgroup position|Sulfatide|SM4(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SM3||Headgroup position|Sulfatide|SM3(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SM2||Headgroup position|Sulfatide|SM2(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SM1||Headgroup position|Sulfatide|SM1a(1) 18:1(4E);3OH/18:0|
|Sphingolipid|SB1||Headgroup position|Sulfatide|SB1b(1) 18:1(4E);3OH/18:0|
|Sphingolipid|L[SP]||Headgroup position|Lysosphingolipid|LGM3(1) 18:1(4E);3OH|

* Glycerol chirality can only use R/S System.
* Serine chirality can only use D/L System.

### Annotation Level
```julia
julia> MBL.annotationlevel(naser)
1-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 completestructurelevel

julia> MBL.annotationlevel(fam)
1-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 fullstructurelevel

julia> MBL.annotationlevel(pc)
3-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 dbpositionlevel
 snpositionlevel
 structuredefinedlevel

julia> MBL.annotationlevel(cl)
1-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 specieslevel

julia> MBL.annotationlevel(gsl)
1-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 completestructurelevel

julia> MBL.annotationlevel(gq1)
1-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 completestructurelevel

julia> MBL.annotationlevel(acer)
2-element Vector{MassSpecBioChemicals.Lipids.LipidAnnotationLevel}:
 structuredefinedlevel
 dbpositionlevel

```
### SMILES
SMILES can only be generated from lipids at full strucure level or above. Full SMILES are available only for simple `FattyAcyl`; for other lipids, complex (polar) groups are dissociated, i.e. only carbon chain parts are presserved, and the carbon chain only SMILES are generated.  
```julia
julia> chemicalsmiles(fam)
"NC(=O)C\\C=C/CCCCC(O)C(O)C(O)CCCC\\C=C/CC(C(=O)N)"

julia> chemicalsmiles(naser; onlycarbonchain = true)
"OC(=O)C\\C=C/CCCCC(O)C(O)C(O)CCCC\\C=C/CC(C(=O)O)"

julia> chemicalsmiles(pc; onlycarbonchain = true)
""

julia> chemicalsmiles(pip2; onlycarbonchain = true)
"C(O)C(OC(=O)CCCCCCCCCCCCCCC)C(OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC)"

julia> chemicalsmiles(gsl; onlycarbonchain = true)
"C(O)C(NC(=O)CCCCCCCCCCCCC\\C=C/CCCCCCCC)C(O)\\C=C\\CCCCCCCCCCCCC"

julia> chemicalsmiles(acer; onlycarbonchain = true)
""

```
### Database
This package includes LMSD database from LIPID MAPS[3].

## Reference
1. G. Liebisch, E. Fahy, J. Aoki, E.A. Dennis, T. Durand, C.S. Ejsing, M. Fedorova, I. Feussner, W.J. Griffiths, H. Kofeler, A.H. Merrill Jr., R.C. Murphy, V.B. O'Donnell, O. Oskolkova, S. Subramaniam, et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures, J. Lipid Res., 61 (2020), pp. 1539-1555
2. Neelamegham S, Aoki-Kinoshita K, Bolton E, Frank M, Lisacek F, Lütteke T, et al. Updates to the symbol nomenclature for Glycans guidelines. Glycobiology. 2019; 29: 620–624.
3. LIPID MAPS: update to databases and tools for the lipidomics community, Conroy MJ, Andrews RM, Andrews S, Cockayne L, Dennis, EA, Fahy E, Gaud C, Griffiths WJ, Jukes G, Kolchin M, Mendivelso K, Lopez-Clavijo AF, Ready C, Subramaniam S, O'Donnell, VB, Nucleic Acids Research, 2023, DOI: 10.1093/nar/gkad896 , PMID: 37855672