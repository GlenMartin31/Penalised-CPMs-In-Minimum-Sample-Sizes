library(tidyverse)

df <- list.files(pattern = ".RDS") %>% map_dfr(read_rds)

write_rds(df, path = "./simulation_results_all.RDS")
