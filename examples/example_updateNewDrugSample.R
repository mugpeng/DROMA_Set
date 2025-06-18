#!/usr/bin/env Rscript

# Example script for DROMA name harmonization functions
# This example demonstrates how to check and harmonize sample/drug names using the DROMA database

# Load required libraries
library(DBI)
library(RSQLite)
library(DROMA.Set)

# If you have the DROMA packages loaded:
# library(DROMA.R)  # or source the functions directly

######################################
# Setup: Connect to DROMA Database
######################################

# Set path to your DROMA database
# db_path <- "path/to/your/droma.sqlite"  # Update this path
db_path <- "250513-DROMA_R/sql_db/droma.sqlite"

# Connect to the database
cat("Connecting to DROMA database...\n")
con <- connectDROMADatabase(db_path)

######################################
# Example 1: Sample Name Harmonization
######################################

cat("\n=== Example 1: Sample Name Harmonization ===\n")

# Example sample names that might come from a new dataset
# These represent various naming conventions and potential issues
example_sample_names <- c(
  "MCF7",           # Should match existing cell line
  "MCF-7",          # Variant with hyphen
  "MCF7 [ATCC]",    # With annotation brackets
  "22RV1",          # Another existing cell line
  "PC3",            # Existing cell line
  "PC-3",           # Variant with hyphen
  "HeLa",           # Common cell line
  "HeLa_S3",        # Variant with suffix
  "NewCellLine1",   # Completely new sample
  "SAMPLE_001",     # Generic new sample
  "Patient_001_T",  # Patient sample tumor
  "Patient_001_N"   # Patient sample normal
)

cat("Original sample names to harmonize:\n")
print(example_sample_names)

# Step 1: Check and harmonize sample names
cat("\nStep 1: Checking sample names against database...\n")
sample_mapping <- checkDROMASampleNames(
  sample_names = example_sample_names,
  connection = con,
  max_distance = 0.2,  # Allow 20% character differences for fuzzy matching
  min_name_length = 5   # Minimum length for partial matching
)

# View the mapping results
cat("\nSample name mapping results:\n")
print(sample_mapping)

# Step 2: Update sample annotations in the database
cat("\nStep 2: Updating sample annotations in database...\n")
updateDROMAAnnotation(
  anno_type = "sample",
  name_mapping = sample_mapping,
  project_name = "MyNewProject",
  data_type = "PDO",  # "CellLine", "PDX", "PDO", "PDC" as appropriate
  tumor_type = "mixed",    # or specific tumor type like "breast cancer"
  connection = con
)

updateDROMAAnnotation("sample", sample_mapping,
                      project_name = "MyProject3",
                      data_type = "CellLine",
                      tumor_type = "breast cancer",
                      Gender = "Female",
                      Age = 1888,
                      FullEthnicity = "European",
                      SimpleEthnicity = "Caucasian")

######################################
# Example 2: Drug Name Harmonization
######################################

cat("\n=== Example 2: Drug Name Harmonization ===\n")

# Example drug names that might come from a new dataset
example_drug_names <- c(
  "Paclitaxel",                    # Should match existing drug
  "paclitaxel",                    # Case variation
  "Taxol",                         # Brand name for paclitaxel
  "Cisplatin",                     # Common chemotherapy drug
  "cis-Diamminedichloroplatinum",  # Chemical name for cisplatin
  "5-FU",                          # Common abbreviation
  "5-Fluorouracil",                # Full name
  "Doxorubicin HCl",               # With salt specification
  "Adriamycin",                    # Brand name for doxorubicin
  "NewDrug_001",                   # Completely new drug
  "Experimental_Compound_X",       # Another new drug
  "((2-(3-Ethoxy-3-oxopropyl)-1,1-dioxido-3-oxo-2,3-dihydro-1,2-benzisothiazol-4-yl)amino)(oxo)acetic acid"  # Very long name
)

cat("Original drug names to harmonize:\n")
print(example_drug_names)

# Step 1: Check and harmonize drug names
cat("\nStep 1: Checking drug names against database...\n")
drug_mapping <- checkDROMADrugNames(
  drug_names = example_drug_names,
  connection = con,
  max_distance = 0.2,              # Allow 20% character differences
  min_name_length = 5,             # Minimum length for partial matching
  keep_long_names_threshold = 18   # Keep very long names as original
)

# View the mapping results
cat("\nDrug name mapping results:\n")
print(drug_mapping)

# Step 2: Update drug annotations in the database
cat("\nStep 2: Updating drug annotations in database...\n")
updateDROMAAnnotation(
  anno_type = "drug",
  name_mapping = drug_mapping,
  project_name = "MyNewProject",
  data_type = "CellLine",  # or "PDX", "PDO", "PDC" as appropriate
  connection = con
)

######################################
# Example 3: Working with Real Data Matrices
######################################

cat("\n=== Example 3: Real Data Matrix Example ===\n")

# Simulate loading a real data matrix
# In practice, you would load your actual data like:
# my_expression_data <- read.csv("my_expression_data.csv", row.names = 1)
# my_drug_response_data <- read.csv("my_drug_response.csv", row.names = 1)

# For this example, create simulated data
set.seed(123)
simulated_expression <- matrix(
  rnorm(100),
  nrow = 10,
  ncol = 10,
  dimnames = list(
    paste0("GENE_", 1:10),  # Gene names as rownames
    c("MCF7", "MCF-7_replicate", "HeLa", "PC3", "22RV1",
      "NewCell_A", "NewCell_B", "Sample_X", "Sample_Y", "Control_1")  # Sample names as colnames
  )
)

simulated_drug_response <- matrix(
  runif(60, min = -3, max = 1),
  nrow = 6,
  ncol = 10,
  dimnames = list(
    c("Paclitaxel", "cisplatin", "5-FU", "Doxorubicin", "NewDrug_A", "NewDrug_B"),  # Drug names as rownames
    c("MCF7", "MCF-7_replicate", "HeLa", "PC3", "22RV1",
      "NewCell_A", "NewCell_B", "Sample_X", "Sample_Y", "Control_1")   # Sample names as colnames
  )
)

cat("Simulated expression data dimensions:", dim(simulated_expression), "\n")
cat("Sample names in expression data:", colnames(simulated_expression), "\n")

cat("Simulated drug response data dimensions:", dim(simulated_drug_response), "\n")
cat("Drug names in drug response data:", rownames(simulated_drug_response), "\n")

# Harmonize sample names from the expression data
cat("\nHarmonizing sample names from expression data...\n")
expr_sample_mapping <- checkDROMASampleNames(
  sample_names = colnames(simulated_expression),
  connection = con
)

# Update sample annotations
updateDROMAAnnotation(
  anno_type = "sample",
  name_mapping = expr_sample_mapping,
  project_name = "ExpressionStudy",
  data_type = "CellLine",
  tumor_type = "mixed",
  connection = con
)

# Harmonize drug names from the drug response data
cat("\nHarmonizing drug names from drug response data...\n")
drug_response_mapping <- checkDROMADrugNames(
  drug_names = rownames(simulated_drug_response),
  connection = con
)

# Update drug annotations
updateDROMAAnnotation(
  anno_type = "drug",
  name_mapping = drug_response_mapping,
  project_name = "DrugScreenStudy",
  data_type = "CellLine",
  connection = con
)

# Now you could rename your actual data matrices using the harmonized names
# colnames(simulated_expression) <- expr_sample_mapping$harmonized_name[match(colnames(simulated_expression), expr_sample_mapping$original_name)]
# rownames(simulated_drug_response) <- drug_response_mapping$harmonized_name[match(rownames(simulated_drug_response), drug_response_mapping$original_name)]
# colnames(simulated_drug_response) <- expr_sample_mapping$harmonized_name[match(colnames(simulated_drug_response), expr_sample_mapping$original_name)]

######################################
# Example 4: Quality Control and Review
######################################

cat("\n=== Example 4: Quality Control and Review ===\n")

# Check for samples/drugs that need manual review
cat("Checking for mappings that need manual review...\n")

# Low confidence sample matches
low_confidence_samples <- sample_mapping[sample_mapping$match_confidence %in% c("medium", "low"), ]
if (nrow(low_confidence_samples) > 0) {
  cat("\nSample mappings requiring manual review:\n")
  print(low_confidence_samples[, c("original_name", "harmonized_name", "match_type", "match_confidence")])
} else {
  cat("No sample mappings require manual review.\n")
}

# Low confidence drug matches
low_confidence_drugs <- drug_mapping[drug_mapping$match_confidence %in% c("medium", "low"), ]
if (nrow(low_confidence_drugs) > 0) {
  cat("\nDrug mappings requiring manual review:\n")
  print(low_confidence_drugs[, c("original_name", "harmonized_name", "match_type", "match_confidence")])
} else {
  cat("No drug mappings require manual review.\n")
}

# Check what was actually added to the database
cat("\nChecking updated annotation tables...\n")

# Get sample annotations for our new project
new_sample_annotations <- getDROMAAnnotation(
  anno_type = "sample",
  project_name = "MyNewProject",
  connection = con
)
cat("New sample annotations added:\n")
print(new_sample_annotations[, c("SampleID", "ProjectID", "DataType", "ProjectRawName")])

# Get drug annotations for our new project
new_drug_annotations <- getDROMAAnnotation(
  anno_type = "drug",
  project_name = "MyNewProject",
  connection = con
)
cat("New drug annotations added:\n")
print(new_drug_annotations[, c("DrugName", "ProjectID", "DataType", "ProjectRawName")])

######################################
# Example 5: Advanced Usage - Batch Processing
######################################

cat("\n=== Example 5: Batch Processing Multiple Projects ===\n")

# Example of processing multiple datasets with different characteristics
projects_info <- list(
  list(
    name = "CellLineProject",
    sample_names = c("MCF7", "HeLa", "A549", "NewLine1"),
    drug_names = c("Paclitaxel", "Cisplatin", "NewDrug1"),
    data_type = "CellLine",
    tumor_type = "mixed"
  ),
  list(
    name = "PDXProject",
    sample_names = c("PDX_001", "PDX_002", "PDX_003"),
    drug_names = c("5-FU", "Doxorubicin", "TargetedTherapy1"),
    data_type = "PDX",
    tumor_type = "breast cancer"
  ),
  list(
    name = "PDOProject",
    sample_names = c("Organoid_A", "Organoid_B", "Organoid_C"),
    drug_names = c("Immunotherapy1", "CombinationDrug"),
    data_type = "PDO",
    tumor_type = "colorectal cancer"
  )
)

# Process each project
for (project in projects_info) {
  cat("\nProcessing project:", project$name, "\n")

  # Harmonize sample names
  sample_map <- checkDROMASampleNames(project$sample_names, connection = con)
  updateDROMAAnnotation("sample", sample_map, project$name, project$data_type, project$tumor_type, con)

  # Harmonize drug names
  drug_map <- checkDROMADrugNames(project$drug_names, connection = con)
  updateDROMAAnnotation("drug", drug_map, project$name, project$data_type, connection = con)

  cat("Completed processing", project$name, "\n")
}

######################################
# Cleanup
######################################

cat("\n=== Cleanup ===\n")

# Close database connection
closeDROMADatabase(con)

cat("\nName harmonization example completed successfully!\n")
cat("\nKey takeaways:\n")
cat("1. Always check sample/drug names before adding new data\n")
cat("2. Review medium/low confidence matches manually\n")
cat("3. Use consistent project names and data types\n")
cat("4. The harmonized names can be used to rename your data matrices\n")
cat("5. All original names are preserved in ProjectRawName columns\n")
cat("6. Batch processing allows efficient handling of multiple projects\n")

######################################
# Additional Notes
######################################

cat("\n=== Additional Usage Notes ===\n")
cat("
# To use these functions with your own data:

# 1. Load your data
my_expression_data <- read.csv('my_expression_data.csv', row.names = 1)
my_drug_data <- read.csv('my_drug_data.csv', row.names = 1)

# 2. Connect to database
con <- connectDROMADatabase('path/to/droma.sqlite')

# 3. Harmonize names
sample_mapping <- checkDROMASampleNames(colnames(my_expression_data), con)
drug_mapping <- checkDROMADrugNames(rownames(my_drug_data), con)

# 4. Review mappings and update database
updateDROMAAnnotation('sample', sample_mapping, 'MyProject', 'CellLine', 'breast cancer', con)
updateDROMAAnnotation('drug', drug_mapping, 'MyProject', 'CellLine', con)

# 5. Apply harmonized names to your data
colnames(my_expression_data) <- sample_mapping$harmonized_name[match(colnames(my_expression_data), sample_mapping$original_name)]
rownames(my_drug_data) <- drug_mapping$harmonized_name[match(rownames(my_drug_data), drug_mapping$original_name)]

# 6. Now you can add the harmonized data to the database
updateDROMADatabase(my_expression_data, 'MyProject_mRNA', overwrite = TRUE)
updateDROMADatabase(my_drug_data, 'MyProject_drug', overwrite = TRUE)

# 7. Close connection
closeDROMADatabase(con)
")
