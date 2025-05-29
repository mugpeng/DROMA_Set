#!/usr/bin/env Rscript

# Using DROMA with SQLite Database
# This example demonstrates how to use the DROMA package with SQLite database
# to efficiently load and query large omics datasets

library(DROMA.Set)

# ---------- Setup ----------
# First-time setup: Create the database from .Rda files
# This step only needs to be done once
cat("Creating DROMA database (this may take a few minutes)...\n")
createDROMADatabase(db_path = "sql_db/droma.sqlite", rda_dir = "data/Rda")

# Connect to the database
db_conn <- connectDROMADatabase()

# ---------- List Available Tables ----------

# View all tables in the database
all_tables <- listDROMADatabaseTables()
cat("Total tables in database:", nrow(all_tables), "\n")

# View just mRNA tables
mrna_tables <- listDROMADatabaseTables(pattern = "_mRNA$")
print(mrna_tables)

# ---------- Example 1: Retrieve Specific mRNA Data ----------

# Retrieve BRCA1 mRNA expression from CCLE and GDSC datasets only
brca1_mRNA <- getFeatureFromDatabase(
  select_feas_type = "mRNA",
  select_feas = c("BRCA1"),
  data_sources = c("CCLE","gCSI")
)

# View the number of samples in each dataset
cat("\nNumber of samples with BRCA1 expression data:\n")
sapply(brca1_mRNA, length)

# ---------- Example 2: Filter by Data Type ----------

# Get a drug's response data only from PDX models
drug_pdx <- getFeatureFromDatabase(
  select_feas_type = "drug",
  select_feas = "Tamoxifen",
  data_type = "PDX"
)

# ---------- Example 3: Filter by Tumor Type ----------

# Get mutation data for TP53 in breast cancer cell lines
tp53_mut_breast <- getFeatureFromDatabase(
  select_feas_type = "mutation_gene",
  select_feas = "TP53",
  data_type = "CellLine",
  tumor_type = "breast cancer"
)

# ---------- Clean Up ----------

# Close the database connection when done
closeDROMADatabase()

cat("\nExample completed successfully\n")
