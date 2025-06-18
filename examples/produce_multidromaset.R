#!/usr/bin/env Rscript

# Using the MultiDromaSet Class for Cross-Project DROMA Analysis
# This example demonstrates how to use the MultiDromaSet class for
# multi-project data access, sample overlap analysis, and cross-project comparisons

library(DROMA.Set)

# Connect to the database
connectDROMADatabase("data/droma.sqlite")
all_tables <- listDROMADatabaseTables()
all_projects <- listDROMAProjects()
cat("Total tables in database:", nrow(all_tables), "\n")
cat("Available projects:", paste(listDROMAProjects(show_names_only = TRUE), collapse = ", "), "\n\n")

# ========== METHOD 1: Create MultiDromaSet from Database ==========

cat("========== Creating MultiDromaSet from Database ==========\n")

# Create MultiDromaSet with specific projects
multi_set1 <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE", "GDSC", "Xeva"),
  db_path = "data/droma.sqlite",
  load_metadata = TRUE
)

# Display information about the MultiDromaSet
print(multi_set1)

# Check available projects
available_projs <- availableProjects(multi_set1)
cat("\nProjects in MultiDromaSet:", paste(available_projs, collapse = ", "), "\n")

# ========== METHOD 2: Create MultiDromaSet from Existing DromaSet Objects ==========

cat("\n========== Creating MultiDromaSet from DromaSet Objects ==========\n")

# First create individual DromaSet objects
gCSI <- createDromaSetFromDatabase("gCSI", db_path = "data/droma.sqlite")
CCLE <- createDromaSetFromDatabase("CCLE", db_path = "data/droma.sqlite")

# Method 2a: Create from list of DromaSet objects
multi_set2 <- createMultiDromaSetFromObjects(list(gCSI = gCSI, CCLE = CCLE))

# Method 2b: Create from individual DromaSet objects as arguments
multi_set2 <- createMultiDromaSetFromObjects(gCSI, CCLE)

cat("Created MultiDromaSet with", length(availableProjects(multi_set2)), "projects\n")

# ========== METHOD 3: Create from All Available Projects ==========

cat("\n========== Creating MultiDromaSet from All Projects ==========\n")

# Create MultiDromaSet with all projects in database
tryCatch({
  multi_set_all <- createMultiDromaSetFromAllProjects(
    db_path = "data/droma.sqlite",
    load_metadata = TRUE
  )
  print(multi_set_all)
}, error = function(e) {
  cat("Note: Could not create from all projects:", e$message, "\n")
  multi_set_all <- multi_set1  # Use the first one as fallback
})

# ========== Sample Overlap Analysis ==========

cat("\n========== Sample Overlap Analysis ==========\n")

# Find overlapping samples between projects
overlap_info <- getOverlappingSamples(multi_set1)
cat("Projects analyzed:", paste(overlap_info$projects, collapse = " vs "), "\n")
cat("Number of overlapping samples:", overlap_info$overlap_count, "\n")
cat("Sample counts per project:\n")
for (i in seq_along(overlap_info$project_sample_counts)) {
  cat("  ", names(overlap_info$project_sample_counts)[i], ":",
      overlap_info$project_sample_counts[i], "samples\n")
}

if (length(overlap_info$overlapping_samples) > 0) {
  cat("First 10 overlapping samples:",
      paste(head(overlap_info$overlapping_samples, 10), collapse = ", "), "\n")
}

# Demonstrate the new sample and drug overlap display in the show method
cat("\nMultiDromaSet object display (showing improved overlap information):\n")
print(multi_set1)

# ========== Drug Overlap Analysis ==========

cat("\n========== Drug Overlap Analysis ==========\n")

# Get drug lists for each project
drug_lists <- list()
for (proj in availableProjects(multi_set1)) {
  ds <- getDromaSet(multi_set1, proj)
  if (nrow(ds@treatmentMetadata) > 0) {
    drug_lists[[proj]] <- ds@treatmentMetadata$DrugName
    cat("Project", proj, "has", length(drug_lists[[proj]]), "drugs\n")
  }
}

# Find overlapping drugs manually for detailed analysis
if (length(drug_lists) >= 2) {
  all_drugs <- unique(unlist(drug_lists))
  common_drugs <- all_drugs

  for (proj in names(drug_lists)) {
    common_drugs <- intersect(common_drugs, drug_lists[[proj]])
  }

  cat("\nDetailed drug overlap analysis:\n")
  cat("Total unique drugs across all projects:", length(all_drugs), "\n")
  cat("Drugs common to all projects:", length(common_drugs), "\n")

  if (length(common_drugs) > 0) {
    cat("Examples of common drugs:", paste(head(common_drugs, 10), collapse = ", "), "\n")
  }

  # Show pairwise drug overlaps
  if (length(drug_lists) == 2) {
    proj_names <- names(drug_lists)
    overlap_drugs <- intersect(drug_lists[[proj_names[1]]], drug_lists[[proj_names[2]]])
    cat("Drugs overlapping between", proj_names[1], "and", proj_names[2], ":", length(overlap_drugs), "\n")
  }
}

# ========== Metadata Analysis ==========

cat("\n========== Metadata Analysis ==========\n")

# Demonstrate the improved metadata with ProjectID labels
cat("Sample metadata now includes ProjectID labels:\n")
if (nrow(multi_set1@sampleMetadata) > 0) {
  sample_counts_by_project <- table(multi_set1@sampleMetadata$ProjectID)
  cat("Sample entries by project:\n")
  for (proj in names(sample_counts_by_project)) {
    cat("  ", proj, ":", sample_counts_by_project[proj], "entries\n")
  }

  # Show first few rows with ProjectID
  cat("\nFirst 5 sample metadata entries:\n")
  print(head(multi_set1@sampleMetadata[, c("SampleID", "ProjectID",
                                          names(multi_set1@sampleMetadata)[1:3])], 5))
}

cat("\nTreatment metadata now includes ProjectID labels:\n")
if (nrow(multi_set1@treatmentMetadata) > 0) {
  drug_counts_by_project <- table(multi_set1@treatmentMetadata$ProjectID)
  cat("Drug entries by project:\n")
  for (proj in names(drug_counts_by_project)) {
    cat("  ", proj, ":", drug_counts_by_project[proj], "entries\n")
  }

  # Show first few rows with ProjectID
  cat("\nFirst 5 treatment metadata entries:\n")
  print(head(multi_set1@treatmentMetadata[, c("DrugName", "ProjectID",
                                             names(multi_set1@treatmentMetadata)[1:3])], 5))
}

# ========== Cross-Project Data Loading ==========

cat("\n========== Cross-Project Data Loading ==========\n")

# Load molecular profiles across projects for overlapping samples
cat("Loading mRNA data across projects...\n")

mRNA_multi2 <- loadMultiProjectMolecularProfiles(
  multi_set1,
  molecular_type = "mRNA",
  data_type = "CellLine"
)

mRNA_multi <- tryCatch({
  loadMultiProjectMolecularProfiles(
    multi_set1,
    molecular_type = "mRNA",
    features = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
    overlap_only = FALSE
  )
}, error = function(e) {
  cat("Could not load mRNA data:", e$message, "\n")
  NULL
})

if (!is.null(mRNA_multi)) {
  cat("Loaded mRNA data from", length(mRNA_multi), "projects:\n")
  for (proj in names(mRNA_multi)) {
    if (is.matrix(mRNA_multi[[proj]]) || is.data.frame(mRNA_multi[[proj]])) {
      cat("  ", proj, ":", nrow(mRNA_multi[[proj]]), "features x",
          ncol(mRNA_multi[[proj]]), "samples\n")
    }
  }
}

# ========== NEW: Load ALL Molecular Profiles ==========

cat("\n========== Loading ALL Molecular Profiles ==========\n")

# Example 1: Load all molecular profiles for a single DromaSet
cat("Loading all molecular profiles for gCSI project...\n")
gCSI_all_mol <- tryCatch({
  loadMolecularProfiles(
    gCSI,
    molecular_type = "all",
    features = c("BRCA1", "BRCA2", "TP53"),  # Limit features for faster loading
    return_data = TRUE
  )
}, error = function(e) {
  cat("Could not load all molecular data for gCSI:", e$message, "\n")
  NULL
})

if (!is.null(gCSI_all_mol)) {
  cat("Loaded", length(gCSI_all_mol), "molecular profile types for gCSI:\n")
  for (mol_type in names(gCSI_all_mol)) {
    data <- gCSI_all_mol[[mol_type]]
    if (is.matrix(data) || is.data.frame(data)) {
      cat("  ", mol_type, ":", nrow(data), "features x", ncol(data), "samples\n")
    } else {
      cat("  ", mol_type, ": data loaded\n")
    }
  }
}

# Example 2: Load all molecular profiles across multiple projects
cat("\nLoading all molecular profiles across projects (overlapping samples only)...\n")
all_mol_multi <- tryCatch({
  loadMultiProjectMolecularProfiles(
    multi_set1,
    molecular_type = "all",
    features = c("BRCA1", "TP53"),  # Limit features for demonstration
    # overlap_only = FALSE
  )
}, error = function(e) {
  cat("Could not load all molecular data across projects:", e$message, "\n")
  NULL
})

if (!is.null(all_mol_multi)) {
  cat("Loaded", length(all_mol_multi), "molecular profile types across projects:\n")
  for (mol_type in names(all_mol_multi)) {
    cat("  Molecular type:", mol_type, "\n")
    proj_data <- all_mol_multi[[mol_type]]
    for (proj in names(proj_data)) {
      data <- proj_data[[proj]]
      if (is.matrix(data) || is.data.frame(data)) {
        cat("    ", proj, ":", nrow(data), "features x", ncol(data), "samples\n")
      } else {
        cat("    ", proj, ": data loaded\n")
      }
    }
  }
}

# Load treatment response data across projects
cat("\nLoading drug response data across projects...\n")
drug_multi <- tryCatch({
  loadMultiProjectTreatmentResponse(
    multi_set1,
    data_type = "treatment",
    drugs = c("Tamoxifen", "Cisplatin", "Paclitaxel"),
    overlap_only = FALSE
  )
}, error = function(e) {
  cat("Could not load drug data:", e$message, "\n")
  NULL
})

if (!is.null(drug_multi)) {
  cat("Loaded drug response data from", length(drug_multi), "projects:\n")
  for (proj in names(drug_multi)) {
    if (is.matrix(drug_multi[[proj]]) || is.data.frame(drug_multi[[proj]])) {
      cat("  ", proj, ":", nrow(drug_multi[[proj]]), "drugs x",
          ncol(drug_multi[[proj]]), "samples\n")
    }
  }
}

# ========== Access Individual DromaSet Objects ==========

cat("\n========== Accessing Individual DromaSet Objects ==========\n")

# Get specific DromaSet from MultiDromaSet
gCSI_from_multi <- getDromaSet(multi_set1, "gCSI")
cat("Retrieved gCSI DromaSet from MultiDromaSet\n")
cat("gCSI project name:", gCSI_from_multi@name, "\n")

# Load data for individual project within MultiDromaSet
gCSI_from_multi <- loadTreatmentResponse(gCSI_from_multi, return_data = FALSE)
cat("Loaded treatment response data for gCSI\n")

# ========== Cross-Project Analysis Example ==========

cat("\n========== Cross-Project Analysis Example ==========\n")

# Example: Compare BRCA1 expression patterns between projects
if (!is.null(mRNA_multi) && length(mRNA_multi) >= 2) {
  cat("Comparing BRCA1 expression between projects...\n")

  # Extract BRCA1 data from each project
  brca1_data <- list()
  for (proj in names(mRNA_multi)) {
    if ("BRCA1" %in% rownames(mRNA_multi[[proj]])) {
      brca1_data[[proj]] <- mRNA_multi[[proj]]["BRCA1", ]
    }
  }

  if (length(brca1_data) >= 2) {
    # Calculate summary statistics for each project
    cat("\nBRCA1 expression summary:\n")
    for (proj in names(brca1_data)) {
      values <- brca1_data[[proj]][!is.na(brca1_data[[proj]])]
      if (length(values) > 0) {
        cat("  ", proj, ": mean =", round(mean(values), 3),
            ", sd =", round(sd(values), 3),
            ", n =", length(values), "\n")
      }
    }

    # Compare between first two projects if possible
    proj_names <- names(brca1_data)[1:2]
    common_samples <- intersect(names(brca1_data[[proj_names[1]]]),
                               names(brca1_data[[proj_names[2]]]))

    if (length(common_samples) >= 5) {
      values1 <- brca1_data[[proj_names[1]]][common_samples]
      values2 <- brca1_data[[proj_names[2]]][common_samples]

      # Remove NA values
      valid_idx <- !is.na(values1) & !is.na(values2)
      if (sum(valid_idx) >= 5) {
        correlation <- cor.test(values1[valid_idx], values2[valid_idx])
        cat("\nCorrelation of BRCA1 expression between", proj_names[1], "and", proj_names[2], ":\n")
        cat("  r =", round(correlation$estimate, 3),
            ", p-value =", format(correlation$p.value, scientific = TRUE, digits = 3), "\n")
        cat("  Based on", sum(valid_idx), "overlapping samples\n")
      }
    }
  }
}

# ========== Drug Response Comparison ==========

if (!is.null(drug_multi) && length(drug_multi) >= 2) {
  cat("\n========== Drug Response Comparison ==========\n")

  # Find common drugs across projects
  all_drugs <- unique(unlist(lapply(drug_multi, rownames)))
  common_drugs <- all_drugs
  for (proj in names(drug_multi)) {
    common_drugs <- intersect(common_drugs, rownames(drug_multi[[proj]]))
  }

  cat("Common drugs across projects:", length(common_drugs), "\n")
  if (length(common_drugs) > 0) {
    cat("Examples:", paste(head(common_drugs, 5), collapse = ", "), "\n")

    # Analyze first common drug if available
    if (length(common_drugs) >= 1) {
      drug_name <- common_drugs[1]
      cat("\nAnalyzing", drug_name, "response across projects:\n")

      drug_responses <- list()
      for (proj in names(drug_multi)) {
        if (drug_name %in% rownames(drug_multi[[proj]])) {
          drug_responses[[proj]] <- drug_multi[[proj]][drug_name, ]
        }
      }

      # Summary statistics
      for (proj in names(drug_responses)) {
        values <- drug_responses[[proj]][!is.na(drug_responses[[proj]])]
        if (length(values) > 0) {
          cat("  ", proj, ": mean =", round(mean(values), 3),
              ", sd =", round(sd(values), 3),
              ", n =", length(values), "\n")
        }
      }
    }
  }
}

# ========== MultiDromaSet Management ==========

cat("\n========== MultiDromaSet Management ==========\n")

# Add a new DromaSet to existing MultiDromaSet
tryCatch({
  gdsc <- createDromaSetFromDatabase("gdsc", db_path = "data/droma.sqlite")
  multi_set_expanded <- addDromaSetToMulti(multi_set1, gdsc, project_name = "GDSC")
  cat("Added GDSC project to MultiDromaSet\n")
  cat("New project count:", length(availableProjects(multi_set_expanded)), "\n")

  # Remove a project
  multi_set_reduced <- removeDromaSetFromMulti(multi_set_expanded, "GDSC")
  cat("Removed GDSC project from MultiDromaSet\n")
  cat("Final project count:", length(availableProjects(multi_set_reduced)), "\n")

}, error = function(e) {
  cat("Note: Could not demonstrate project management:", e$message, "\n")
})

# ========== Practical Analysis Workflow ==========

cat("\n========== Practical Analysis Workflow ==========\n")

# Workflow: Find drug-gene associations across projects
cat("Example workflow: Finding drug-gene associations across projects\n")

# 1. Load specific molecular and drug data
workflow_features <- c("TP53", "BRCA1", "EGFR")
workflow_drugs <- c("Tamoxifen", "Cisplatin")

cat("Loading data for analysis...\n")
mol_data <- tryCatch({
  loadMultiProjectMolecularProfiles(
    multi_set1,
    molecular_type = "mRNA",
    features = workflow_features,
    overlap_only = FALSE
  )
}, error = function(e) {
  cat("Could not load molecular data for workflow\n")
  NULL
})

drug_data <- tryCatch({
  loadMultiProjectTreatmentResponse(
    multi_set1,
    data_type = "treatment",
    drugs = workflow_drugs,
    overlap_only = FALSE
  )
}, error = function(e) {
  cat("Could not load drug data for workflow\n")
  NULL
})

# 2. Perform correlation analysis if data is available
if (!is.null(mol_data) && !is.null(drug_data) &&
    length(mol_data) > 0 && length(drug_data) > 0) {

  cat("\nPerforming correlation analysis...\n")

  for (proj in names(mol_data)) {
    if (proj %in% names(drug_data)) {
      mol_proj <- mol_data[[proj]]
      drug_proj <- drug_data[[proj]]

      if (is.matrix(mol_proj) && is.matrix(drug_proj) &&
          ncol(mol_proj) > 0 && ncol(drug_proj) > 0) {

        # Find common samples
        common_samples <- intersect(colnames(mol_proj), colnames(drug_proj))

        if (length(common_samples) >= 10) {
          cat("  Project", proj, ": analyzing", length(common_samples), "common samples\n")

          # Example: correlate first gene with first drug
          if (nrow(mol_proj) > 0 && nrow(drug_proj) > 0) {
            gene_values <- mol_proj[1, common_samples]
            drug_values <- drug_proj[1, common_samples]

            # Remove missing values
            valid_idx <- !is.na(gene_values) & !is.na(drug_values)
            if (sum(valid_idx) >= 5) {
              cor_result <- cor.test(gene_values[valid_idx], drug_values[valid_idx])
              cat("    ", rownames(mol_proj)[1], "vs", rownames(drug_proj)[1],
                  ": r =", round(cor_result$estimate, 3),
                  ", p =", format(cor_result$p.value, scientific = TRUE, digits = 2), "\n")
            }
          }
        }
      }
    }
  }
}

# ========== Summary ==========

cat("\n========== Analysis Summary ==========\n")
cat("MultiDromaSet analysis completed successfully!\n")
cat("Key capabilities demonstrated:\n")
cat("1. Creating MultiDromaSet objects from database and existing DromaSet objects\n")
cat("2. Sample overlap detection and analysis\n")
cat("3. Cross-project data loading with automatic filtering\n")
cat("4. Individual DromaSet access within MultiDromaSet\n")
cat("5. Cross-project correlation and comparison analysis\n")
cat("6. MultiDromaSet management (adding/removing projects)\n")
cat("7. Practical analysis workflows for drug-gene associations\n")
cat("8. NEW: Loading ALL molecular profiles using molecular_type='all'\n")
cat("9. Enhanced metadata handling with ProjectID labels\n")
cat("10. Comprehensive drug overlap analysis\n")

# ========== Clean Up ==========

# Close the database connection when done
closeDROMADatabase()

cat("\nMultiDromaSet example completed successfully!\n")
