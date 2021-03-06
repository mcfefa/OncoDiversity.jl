# OncoDiversity

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mcfefa.github.io/OncoDiversity.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mcfefa.github.io/OncoDiversity.jl/dev)
[![Build Status](https://github.com/mcfefa/OncoDiversity.jl/workflows/CI/badge.svg)](https://github.com/mcfefa/OncoDiversity.jl/actions)

Calculating the generalized diversity index for individual patients across a variety of different cancer datasets, including:

- clustered single cell data
- individual CDR3 sequences recovered from bulk sequencing files (BAM files)
- SomaScan analyte intensities from blood samples 


### Installation
``` 
Pkg.develop(PackageSpec(url="https://github.com/mcfefa/OncoDiversity.jl"))

using OncoDiversity

```

