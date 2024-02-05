## code to prepare `DATASET` dataset goes here
library(tidyverse)
generated_count_data <- read_tsv("src/generated_count_data.tsv")
usethis::use_data(generated_count_data, overwrite = TRUE)

generated_data_input_pSups <- read_tsv("src/generated_data_input_pSups.tsv")
usethis::use_data(generated_data_input_pSups, overwrite = TRUE)

generated_data_samplesheet <- read_tsv("src/generated_data_samplesheet.tsv")
usethis::use_data(generated_data_samplesheet, overwrite = TRUE)
