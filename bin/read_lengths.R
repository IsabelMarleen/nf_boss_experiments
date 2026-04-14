#!/usr/bin/env Rscript

library(readr)
library(dplyr)

father <- read_delim(
    "/hps/software/users/goldman/ipoetzsch/boss_exp/barcoded_hg002_hg003_hg004_chr21/analyse/read_lengths/read_length_control_5_bc7_hschr21.csv",
    delim = " ",
    col_names = c("length", "n"))%>%
    slice_head(n=-1)

mother <- read_delim(
    "/hps/software/users/goldman/ipoetzsch/boss_exp/barcoded_hg002_hg003_hg004_chr21/analyse/read_lengths/read_length_control_5_bc30_hschr21.csv",
    delim = " ",
    col_names = c("length", "n"))%>%
    slice_head(n=-1)

son <- read_delim(
    "/hps/software/users/goldman/ipoetzsch/boss_exp/barcoded_hg002_hg003_hg004_chr21/analyse/read_lengths/read_length_control_5_bc15_hschr21.csv",
    delim = " ",
    col_names = c("length", "n"))%>%
    slice_head(n=-1)

read_lengths <- bind_rows("father" = father, "mother" = mother, "son"=son, .id = "person") %>%
    group_by(person) %>%
    mutate(length = as.numeric(length), n = as.numeric(n))%>%
    mutate(intermed_mean = length*n, intermed_var = (length - mean(intermed_mean))^2*n)%>%
    summarise(mean = mean(intermed_mean), med = median (intermed_mean), var = sum(intermed_var)/n())
    %>% print()
