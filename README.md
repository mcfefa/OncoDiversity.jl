# OncoDiversity

## Citation

If you are using OncoDiversity.jl in your research, please cite this paper: 

MC Ferrall-Fairbanks, NH Chakiryan, BI Chobrutskiy, Y Kim, JK Teer, A Berglund, JJ Mule, M Fournier, EM Siegel, J Dhillon, SSA Falasiri, JF Arturo, EN Katende, G Blanck, BJ Manley, PM Altrock. (2022) _Cancer Res._ 2022 Mar 1; 82(5):929-942. doi: [10.1158/0008-5472.CAN-21-1747](https://aacrjournals.org/cancerres/article/82/5/929/681768/Quantification-of-T-and-B-cell-Immune-Receptor). PMID: [35031572](https://pubmed.ncbi.nlm.nih.gov/35031572/). PMCID: [PMC8898280](https://pmc.ncbi.nlm.nih.gov/articles/PMC8898280/).

## Cancer Genomics Cloud Integration

The Seven Bridges [Cancer Genomics Cloud (CGC)](https://www.cancergenomicscloud.org/), powered by Velsera and funded by the NCI, is a flexible cloud platform that enables analysis, storage, and computation of large cancer datasets. The CGC provides a user-friendly portal to access and analyze cancer data where it lives. We developed a cloud-based workflow performs end-to-end processing from RNA-seq data retrieval through CDR3 V(D)J sequence recovery and quantification using OncoDiversity.jl by connecting multiple applications within the Cancer Genomics Cloud (CGC) environment and allows users to directly leverage short read archive (SRA) accession numbers for publicly available cohorts. 

**CDR3 V(D)J Recovery Diversity Public Project** is available at: https://cgc.sbgenomics.com/u/sevenbridges/cdr3-v-d-j-recovery. 

## OncoDiversity.jl Package

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

