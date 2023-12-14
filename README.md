# musicParam

`lute` param class, method, and generic definitions for the DeconRNAseq deconvolution algorithm. See `?deconrnaseqParam-class` for more information.

## Install

Install this package from GitHub using:

```
devtools::install_github("metamaden/deconrnaseqParam")
```

## Dependencies

This param class requires the `DeconRNASeq` software to run (available from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html)).

A YML file has been included to set up a conda environment to run the main dependencies `./inst/yml/deconrnaseq.yml`.

## Citations

Racle, Julien, and David Gfeller. 2020. “EPIC: A Tool to Estimate the Proportions of Different Cell Types from Bulk Gene Expression Data.” In 
Bioinformatics for Cancer Immunotherapy: Methods and Protocols, edited by Sebastian Boegel, 233–48. Methods in Molecular Biology. New York, NY: 
Springer US. https://doi.org/10.1007/978-1-0716-0327-7_17.