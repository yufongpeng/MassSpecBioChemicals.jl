# MassSpecBioChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecBioChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=test)](https://github.com/yufongpeng/MassSpecBioChemicals.jl/actions/workflows/CI.yml?query=branch%3Atest)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecBioChemicals.jl/branch/test/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecBioChemicals.jl)

MassSpecBioChemicals.jl is package for representation of biochemicals with four submodules, `Lipids`, `Glycans`, `Proteins`, and `Nucleotides`. It is built on [MassSpecChemicals.jl](https://github.com/yufongpeng/MassSpecChemicals.jl); however, it is still in early-stage, so interfaces are not fully implemented. Currently, the development is focused on submodule `Lipids`, especially the tranformation of standard lipid nomenclatures to internal lipid objects, generation of SMILES, and manipulation of internal lipid objects. The following features to be implemented is the transformation of lipids at less exact annotation level to full structure level to generate SMILES. The transformation will be based on either public lipid libraries filtered by taxonomy Information, or user-defined fatty acid composition. 

## Reference
1. G. Liebisch, E. Fahy, J. Aoki, E.A. Dennis, T. Durand, C.S. Ejsing, M. Fedorova, I. Feussner, W.J. Griffiths, H. Kofeler, A.H. Merrill Jr., R.C. Murphy, V.B. O'Donnell, O. Oskolkova, S. Subramaniam, et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures, J. Lipid Res., 61 (2020), pp. 1539-1555
2. Neelamegham S, Aoki-Kinoshita K, Bolton E, Frank M, Lisacek F, Lütteke T, et al. Updates to the symbol nomenclature for Glycans guidelines. Glycobiology. 2019; 29: 620–624.