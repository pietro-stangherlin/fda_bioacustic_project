# subsample gufi due to computational issues
# (too many observations) -> memory problem while 
# computing ANOVA se

load("data/data_1.RData")

table(gufi$Climate_zone)

set.seed(123)

# keep only 150

gufi_ids_exclude = sample(which(gufi$Climate_zone == "Temperate"),
                               496 - 150)

write.table(gufi_ids_exclude,
            "data/gufi_ids_exclude.txt",
            row.names = FALSE)

