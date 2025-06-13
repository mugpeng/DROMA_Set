#!/usr/bin/env Rscript

# Example: Using updateDROMAProjects to manage project metadata
# This example demonstrates how to use the updateDROMAProjects function
# to automatically manage project metadata in the DROMA database

library(DROMA.Set)

# Example workflow for updating project metadata
# ============================================

# Step 1: Connect to database
# ---------------------------
# First, connect to your DROMA database
db_path <- "path/to/your/droma.sqlite"
con <- connectDROMADatabase(db_path)

# Step 2: Add data tables using updateDROMADatabase
# ------------------------------------------------
# Create some example data (in real usage, you'd load your actual data)

# Example mRNA expression data
# set.seed(123)
# mrna_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
# rownames(mrna_data) <- paste0("Gene_", 1:100)
# colnames(mrna_data) <- paste0("Sample_", 1:10)
#
# # Example CNV data
# cnv_data <- matrix(sample(c(-2, -1, 0, 1, 2), 500, replace = TRUE), nrow = 50, ncol = 10)
# rownames(cnv_data) <- paste0("Gene_", 1:50)
# colnames(cnv_data) <- paste0("Sample_", 1:10)
#
# # Example drug response data
# drug_data <- matrix(runif(50, 0, 10), nrow = 5, ncol = 10)
# rownames(drug_data) <- paste0("Drug_", 1:5)
# colnames(drug_data) <- paste0("Sample_", 1:10)
#
# # Add data to database for project "gCSI"
# updateDROMADatabase(mrna_data, "gCSI_mRNA", overwrite = TRUE)
# updateDROMADatabase(cnv_data, "gCSI_cnv", overwrite = TRUE)
# updateDROMADatabase(drug_data, "gCSI_drug", overwrite = TRUE)
#
# # Add data to database for project "CCLE"
# updateDROMADatabase(mrna_data, "CCLE_mRNA", overwrite = TRUE)
# updateDROMADatabase(drug_data, "CCLE_drug", overwrite = TRUE)

# Step 3: Update project metadata
# ------------------------------

# Option 1: Update metadata for a specific project
# cat("Updating metadata for gCSI project...\n")
# updateDROMAProjects("gCSI")
#
# # Option 2: Update metadata for all projects in the database
# cat("\nUpdating metadata for all projects...\n")
# updateDROMAProjects()
#
# # Step 4: View the updated projects table
# # ---------------------------------------
# cat("\nCurrent projects in database:\n")
# projects_info <- listDROMAProjects()
# print(projects_info)

# Step 5: Add more data and update again
# -------------------------------------
# Add protein data for gCSI
protein_data <- matrix(rnorm(300), nrow = 30, ncol = 10)
rownames(protein_data) <- paste0("Protein_", 1:30)
colnames(protein_data) <- paste0("Sample_", 1:10)

updateDROMADatabase(protein_data, "test_proteinrppa", overwrite = TRUE)

# Update project metadata again
cat("\nUpdating gCSI metadata after adding protein data...\n")
updateDROMAProjects("test")

# View updated information
cat("\nUpdated projects information:\n")
projects_info_updated <- listDROMAProjects()
print(projects_info_updated)

# Step 6: Demonstrate error handling
# ---------------------------------
cat("\nTesting error handling...\n")

# Try to update a project that doesn't exist in the database
tryCatch({
  updateDROMAProjects("NonExistentProject")
}, warning = function(w) {
  cat("Warning caught:", w$message, "\n")
})

# Step 7: Clean up
# ---------------
closeDROMADatabase()

cat("\nExample completed successfully!\n")

# Additional Usage Examples
# ========================

# Example with sample annotation data
# ----------------------------------
if (FALSE) {  # Set to TRUE to run this part

  # Create sample annotation data
  sample_anno <- data.frame(
    SampleID = paste0("Sample_", 1:10),
    ProjectID = "gCSI",
    DataType = "CellLine",
    TumorType = c(rep("breast cancer", 5), rep("lung cancer", 5)),
    stringsAsFactors = FALSE
  )

  # Add sample annotation to database
  updateDROMADatabase(sample_anno, "sample_anno", overwrite = TRUE)

  # Update project metadata (now with dataset type information)
  updateDROMAProjects("gCSI")

  # View the enhanced project information
  projects_with_anno <- listDROMAProjects()
  print(projects_with_anno)
}

# Example of batch updating multiple projects
# ------------------------------------------
if (FALSE) {  # Set to TRUE to run this part

  # Add data for multiple projects
  projects <- c("GDSC", "CTRPv2", "PRISM")

  for (proj in projects) {
    # Add some example data
    example_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(example_data) <- paste0("Feature_", 1:10)
    colnames(example_data) <- paste0("Sample_", 1:10)

    updateDROMADatabase(example_data, paste0(proj, "_mRNA"), overwrite = TRUE)
    updateDROMADatabase(example_data, paste0(proj, "_drug"), overwrite = TRUE)
  }

  # Update all projects at once
  updateDROMAProjects()

  # View all projects
  all_projects <- listDROMAProjects()
  print(all_projects)
}

# Tips for using updateDROMAProjects:
# ==================================
# 1. Always run updateDROMAProjects after adding new data tables
# 2. Use updateDROMAProjects() without arguments to update all projects
# 3. The function automatically creates the projects table if it doesn't exist
# 4. Project metadata includes: data types, sample count, drug count, timestamps
# 5. The function handles both new projects and updates to existing ones
# 6. Sample and drug counts are automatically calculated from the data tables
# 7. Dataset type is inferred from sample_anno table if available
