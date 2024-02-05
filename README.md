# SedSeqQuant - RNA Fraction Quantification

Quantify the distribution of RNAs from cellular fractions.

------------------------------------------------------------------------

# Introduction

This package employs a bayesian statistical model to quantify the distribution of mRNA

------------------------------------------------------------------------

## 1.1 Required input files and their format

Users need to provide a count matrix and a samplesheet file. Those files need to meet the conditions below:

**1. Both the count matrix and samplesheet file need to be in the ".tsv" format**

**2. The count matrix should have a column of transcript_IDs then a column of counts for every Sample_ID. It is highly recommended to first filter this for verified ORFs. Counts can be integers or fraction estimates (e.g. est_counts from kallisto).**

| transcript_ID | F01~snake_230331 | F02~snake_230331 |
|---------------|------------------|------------------|
| YAL068C       | 1.08             | 3.41             |
| YAL067C       | 17               | 15               |
| YAL064W       | 16               | 24               |

**3. In samplesheet file, there should be at least four columns "Sample_ID","Condition","Rep",and "Fraction". Every Sample_ID must be unique. For every Condition and Rep, there should be a Sample_ID corresponding to the three fractions: Total, Supernatant, and Pellet **


| Sample_ID       | Condition                      | Rep    | Fraction |
|-----------------|--------------------------------|--------|----------|
| F01~snake_230331| hairpinReporters~30C~NA~none   | yHG010 | Total    |
| F02~snake_230331| hairpinReporters~30C~NA~none   | yHG005 | Total    |

**4. For every line of comment, there should be a tag "\#" at the begining of the text**

# Installation

We recommend you to use package "devtools" for dowloading this package from GitHub. Please refer [devtools installation instructions](https://www.r-project.org/nosvn/pandoc/devtools.html) for more information.

```         
install.packages("devtools")
library(devtools)
install_github("jabard89/SedSeqQuant")
```

Then load SedSeqQuant as a standard package:

```         
library(SedSeqQuant)
```

------------------------------------------------------------------------

## Quick guide

```         
count_data <- load_counts(count_file=file.path("count_matrix.tsv.gz"),
                          samplesheet = file.path("samplesheet.tsv"))
pd_noreps <- prepare_data_noreps(count_data,min_counts=20)
stanfit_noreps <- model_fit_noreps(pd_noreps)
sum_noreps <- get_stan_summary_noreps(stanfit_noreps,pd_noreps)
write_stan_summary_noreps(sum_noreps,file.path(wrk.dir,"output/no_reps"))

pd_withreps <- prepare_data_with_reps(count_data,min_counts=20)
stanfit_withreps <- model_fit_reps(pd_withreps,iter=8000,chains=4)
sum_withreps <- get_stan_summary_reps(stanfit_withreps,pd_withreps)
write_stan_summary_with_reps(sum_withreps,file.path(wrk.dir,"output/with_reps"))
```
