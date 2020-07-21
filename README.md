# KnockoffZoom v2

A powerful and versatile statistical method for the analysis of genome-wide association data.

Accompanying paper:
> *Controlling the false discovery rate in GWAS with diverse and related samples* <br />
> M. Sesia, S. Bates, E. Candès, J. Marchini, C. Sabatti <br />
> bioRxiv preprint

For more information, visit: [https://msesia.github.io/knockoffzoom](https://msesia.github.io/knockoffzoom).

## Overview

The goal of *KnockoffZoom* is to identify causal variants for complex traits effectively and precisely through genome-wide fine-mapping, accounting for linkage disequilibrium and controlling the false discovery rate.
The results leverage the genetic models used for phasing and are equally valid for quantitative and binary traits.


The code contained in this repository is designed to allow the application of *KnockoffZoom* to large datasets, such as the [UK Biobank](https://www.ukbiobank.ac.uk/).
Some of the code is provided in the form of Bash and R scripts, while the core algorithms for Monte Carlo knockoff sampling are implemented in the R package [SNPknock](https://github.com/msesia/snpknock/), which should be installed separately.

The *KnockoffZoom* methodology is divided into different modules, each corresponding to a separate Bash script contained in the directory `knockoffzoom/`.

## Dependencies

Recommended OS: Linux. Mac OS is not supported but should be compatible.

The following software should be available from your user path:

   - [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
   - [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) alpha
   - [fastPHASE](http://scheet.org/software.html) 1.4.8
   - [GNU datamash](https://www.gnu.org/software/datamash/) 1.3
   - [GNU awk](https://github.com/onetrueawk/awk) 4.0.2
   - [GNU Core Utilities](https://www.gnu.org/software/coreutils/) 8.22

The following [R](https://www.r-project.org/) (version 3.5.1) packages are required:

   - [SNPknock](https://msesia.github.io/snpknock/) 0.8.2
   - [adjclust](https://CRAN.R-project.org/package=adjclust ) 0.5.6
   - [bigsnpr](https://privefl.github.io/bigsnpr/) 0.9.1
   - [bigstatsr](https://privefl.github.io/bigstatsr/) 0.8.4
   - [snpStats](https://doi.org/doi:10.18129/B9.bioc.snpStats) 1.32
   - [Matrix](https://CRAN.R-project.org/package=Matrix ) 1.2.15
   - [data.table](https://CRAN.R-project.org/package=data.table) 1.12.0
   - [tidyverse](https://www.tidyverse.org/) 1.2.1
   - [devtools](https://CRAN.R-project.org/package=devtools) 1.13.6

The above version numbers correspond to the configuration on which this software was tested. Newer version are likely to be compatible, but have not been tested.

## Installation

Clone this repository on your system and install any missing dependencies. Estimated installation time (dependencies): 5-15 minutes.

## Toy dataset

A toy dataset containing 1000 artificial samples typed at 2000 loci (divided between chromosome 21 and 22) is offered as a toy example to test *KnockoffZoom*. To run the example, simply execute the script `analyze.sh`.

```{bash}
./analyze.sh
```

This script will also verify whether required R packages are available and install them otherwise.

The analysis should take approximately 5 minutes on a personal computer. The results can be visualized interactively with the script `visualize.sh`, which will launch a [Shiny](https://shiny.rstudio.com/) app in your browser. Some additional R packages are required by the visualization tool, and will be automatically installed if not found.

```{bash}
./visualize.sh
```

The expected results for the analysis of this toy dataset are provided in the directory `results/` and can be visualized by running the script `visualize.sh` before running `analyze.sh`. Note that the script `analyze.sh` will overwrite the default results. 

## Tutorial

A guided step-by-step analysis of the above toy dataset using *KnockoffZoom* is available at:
[https://msesia.github.io/knockoffzoom/tutorial.html](https://msesia.github.io/knockoffzoom/tutorial.html).

## Large-scale applications

*KnockoffZoom* is computationally efficient and we have successfully applied it to the analysis of the genetic data in the UK Biobank. For more information, visit [https://msesia.github.io/knockoffzoom/ukbiobank.html](https://msesia.github.io/knockoffzoom/ukbiobank.html).
The analysis of large datasets cannot be carried out on a personal computer. The computational resources required for the analysis of the UK Biobank data are summarized in the [accompanying paper](https://doi.org/10.1101/631390).

The modular nature of our method allows the code contained in each of the 5 main scripts to be easily deployed on a computing cluster for large-scale applications. This task will require some additional user effort compared to the toy example, but the scripts for each module are documented and quite intuitive.



## Authors

   - [Matteo Sesia](https://msesia.github.io/) (Stanford University).

## Contributors

   - [Stephen Bates](https://stephenbates19.github.io/) (Stanford University).

## License

This software is distributed under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Further references

Read more about the broader framework of [knockoffs](https://web.stanford.edu/group/candes/knockoffs/).