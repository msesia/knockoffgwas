# KnockoffGWAS

 [![Build Status](https://travis-ci.org/msesia/knockoffgwas.svg?branch=master)](https://travis-ci.org/msesia/knockoffgwas)

A powerful and versatile statistical method for the analysis of genome-wide association data with population structure.
This method localizes causal variants while controlling the false discovery rate, and is valid even if the samples have diverse ancestries and familial relatedness.

Accompanying paper:
> *Controlling the false discovery rate in GWAS with population structure* <br />
> M. Sesia, S. Bates, E. Candès, J. Marchini, C. Sabatti <br />
> preprint at bioRxiv; https://doi.org/10.1101/2020.08.04.236703

For more information, visit: [https://msesia.github.io/knockoffgwas](https://msesia.github.io/knockoffgwas).

For an earlier version of this method restricted to homogeneous populations, see also [KnockoffZoom](https://github.com/msesia/knockoffzoom).


## Overview

The goal of *KnockoffGWAS* is to identify causal variants for complex traits effectively and precisely through genome-wide fine-mapping, accounting for linkage disequilibrium and controlling the false discovery rate.
The results leverage the genetic models used for phasing and are equally valid for quantitative and binary traits.
The main innovation *KnockoffGWAS* is to support the analysis of diverse populations, with different ancestries and possibly close familial relatedness.
Furthermore, *KnockoffGWAS* includes a highly efficient standalone C++ program for generating genetic knockoffs for large data sets, which facilitates applications compared to *KnockoffZoom*.

The code contained in this repository is designed to allow the application of *KnockoffGWAS* to large datasets, such as the [UK Biobank](https://www.ukbiobank.ac.uk/).
Some of the code is provided in the form of Bash and R scripts, while the core algorithms for Monte Carlo knockoff sampling are implemented in C++.

The *KnockoffGWAS* methodology is divided into different modules, each corresponding to a separate Bash script contained in the directory `knockoffgwas/`.

## Dependencies

Recommended OS: Linux. Mac OS is not supported but should be compatible.

The following software should be available from your user path:

   - [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)

The following [R](https://www.r-project.org/) (version 4.0.2) packages are required:

   - [fastcluster](https://CRAN.R-project.org/package=fastcluster ) 1.1.25
   - [bigsnpr](https://privefl.github.io/bigsnpr/) 1.4.4
   - [bigstatsr](https://privefl.github.io/bigstatsr/) 1.2.3
   - [tidyverse](https://www.tidyverse.org/) 1.3.0
   
The above version numbers correspond to the configuration on which this software was tested. Newer version are likely to be compatible, but have not been tested.

## Installation

Clone this repository on your system and install any missing dependencies. Estimated installation time (dependencies): 5-15 minutes.
Compile the C++ program for knockoff generation by entering the directory `snpknock2` and running `make`.

## Toy dataset

A toy dataset containing 1000 artificial samples typed at 2000 loci (divided between chromosome 21 and 22) is offered as a toy example to test *KnockoffGWAS*. To run the example, simply execute the script `analyze.sh`.

```{bash}
./analyze.sh
```

This script will also verify whether required R packages are available and install them otherwise.

The analysis should take less than 5 minutes on a personal computer. The results can be visualized interactively with the script `visualize.sh`, which will launch a [Shiny](https://shiny.rstudio.com/) app in your browser. Some additional R packages are required by the visualization tool, and will be automatically installed if not found.

```{bash}
./visualize.sh
```

The expected results for the analysis of this toy dataset are provided in the directory `results/` and can be visualized by running the script `visualize.sh` before running `analyze.sh`. Note that the script `analyze.sh` will overwrite the default results. 

## Large-scale applications

*KnockoffGWAS* is computationally efficient and we have successfully applied it to the analysis of the genetic data in the UK Biobank. For more information, visit [https://msesia.github.io/knockoffgwas/ukbiobank.html](https://msesia.github.io/knockoffgwas/ukbiobank.html).
The analysis of large datasets cannot be carried out on a personal computer. The computational resources required for the analysis of the UK Biobank data are summarized in the [accompanying paper](https://doi.org/10.1101/2020.08.04.236703).

The modular nature of our method allows the code contained in each of the 4 main scripts to be easily deployed on a computing cluster for large-scale applications. This task will require some additional user effort compared to the toy example, but the scripts for each module are documented and quite intuitive.


## Authors

   - [Matteo Sesia](https://msesia.github.io/) (University of Southern California).

## Contributors

   - [Stephen Bates](https://stephenbates19.github.io/) (Stanford University).

## License

This software is distributed under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Further references

Read more about:
 - [KnockoffGWAS](https://github.com/msesia/knockoffgwas);
 - the broader framework of [knockoffs](https://web.stanford.edu/group/candes/knockoffs/).
