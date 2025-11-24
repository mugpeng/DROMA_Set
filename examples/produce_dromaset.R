#!/usr/bin/env Rscript

# Using the DromaSet Class for DROMA
# This example demonstrates how to use the new project-oriented DromaSet class
# for more efficient data access and organization

library(DROMA.Set)

# Connect to the database
connectDROMADatabase("data/droma.sqlite")
all_tables <- listDROMADatabaseTables()
all_projects <- listDROMAProjects()
cat("Total tables in database:", nrow(all_tables), "\n")

# ---------- Create DromaSet Objects ----------

# Get just the project names
project_names <- listDROMAProjects(show_names_only = TRUE)

# Create a DromaSet object for gCSI data
gCSI <- createDromaSetFromDatabase("gCSI", db_path = "data/droma.sqlite")

# load data include drug and molecule
# gCSI_2 <- createDromaSetFromDatabase("gCSI", "data/droma.sqlite", auto_load = TRUE)

# Display information about the DromaSet
print(gCSI)

# Load drug response data
gCSI <- loadTreatmentResponse(gCSI)
# directly return data
gCSI_drug <- loadTreatmentResponse(gCSI, return_data = T)

# ---------- Check Annotation ----------
# Get all sample annotations
sample_anno <- getDROMAAnnotation("sample")

# Get sample annotations for gCSI project only
gCSI_samples <- getDROMAAnnotation("sample", project_name = "gCSI")

# Get annotations for specific samples
specific_samples <- getDROMAAnnotation("sample", ids = c("22RV1", "2313287"))

# Get only cell line samples
cell_lines <- getDROMAAnnotation("sample", data_type = "CellLine")

# Get only breast cancer samples
breast_samples <- getDROMAAnnotation("sample", tumor_type = "breast cancer")

# Get all drug annotations
drug_anno <- getDROMAAnnotation("drug")

# Get drug annotations for gCSI project
gCSI_drugs <- getDROMAAnnotation("drug", project_name = "gCSI")

# Get annotations for specific drugs
specific_drugs <- getDROMAAnnotation("drug", ids = c("Tamoxifen", "Cisplatin"))

# Combine filters - get first 100 cell line samples from gCSI
gCSI_cell_lines <- getDROMAAnnotation("sample",
                                      project_name = "gCSI",
                                      data_type = "CellLine",
                                      limit = 100)

# ---------- Check Specific Project, data type, feature, drug, sample ----------
# Check what molecular profiles are available
profiles <- availableMolecularProfiles(gCSI)
cat("\nAvailable molecular profiles for gCSI:\n")
print(profiles)

# Check what treatment response data is available
responses <- availableTreatmentResponses(gCSI)
cat("\nAvailable treatment response data for gCSI:\n")
print(responses)

# Get data types for gCSI project
gCSI_data_types <- listDROMAProjects(projects = "gCSI")

# get anno data
# check samples from anno data
gCSI_anno_sample <- gCSI@sampleMetadata
gCSI_anno_drug <- gCSI@treatmentMetadata

# List all genes in gCSI mRNA data
genes <- listDROMAFeatures("gCSI", "cnv")

# List all drugs in gCSI drug response data
drugs <- listDROMAFeatures("gCSI", "drug")

# Find BRCA genes
brca_genes <- listDROMAFeatures("gCSI", "cnv", pattern = "^BRCA")

# Get first 100 features
top_genes <- listDROMAFeatures("gCSI", "mRNA", limit = 100)

# Select samples
samples_sel <- listDROMASamples("gCSI")

# ---------- Load all Data ----------
CCLE <- createDromaSetFromDatabase("CCLE", db_path = "data/droma.sqlite")
CCLE <- loadMolecularProfiles(
  CCLE, molecular_type = "mRNA"
)
CCLE <- loadTreatmentResponse(
  CCLE
)

# ---------- Load Specific Data ----------
# Partially Load Copy Number data for specific genes
gCSI_cnv_sel <- loadMolecularProfiles(
  gCSI,
  molecular_type = "cnv",
  features = c("BRCA1", "BRCA2", "TP53", "cf999"),
  samples = c("2313287", "NCIH522", "22RV1"),
  return_data = T
)

gCSI_mut_sel <- loadMolecularProfiles(
  gCSI,
  molecular_type = "mutation_gene",
  features = c("BRCA1"),
  # samples = c("2313287", "NCIH522", "22RV1"),
  return_data = T
)

# Partially Load drug for samples
gCSI_drug_sel <- loadTreatmentResponse(
  gCSI,
  tumor_type = "thyroid cancer",
  drugs = c("Azacitidine", "Idelalisib", "sd"),
  # data_type = "PDO",
  return_data = T)

# ---------- Example Analysis ----------

# Analyze relationship between BRCA1 copy number and drug response
if ("cnv" %in% names(gCSI@molecularProfiles) && "drug" %in% names(gCSI@treatmentResponse)) {
  # Get BRCA1 CNV data
  brca1_cnv <- gCSI@molecularProfiles$cnv["BRCA1", ]

  # Get response to a specific drug (e.g., first drug in the matrix)
  drug_name <- rownames(gCSI@treatmentResponse$drug)[1]
  drug_response <- gCSI@treatmentResponse$drug[drug_name, ]

  # Find common samples
  common_samples <- intersect(names(brca1_cnv), names(drug_response))

  # Calculate correlation if there are enough samples
  if (length(common_samples) >= 10) {
    correlation <- cor.test(
      brca1_cnv[common_samples],
      drug_response[common_samples],
      method = "pearson"
    )

    cat("\nCorrelation between BRCA1 CNV and", drug_name, "response:\n")
    cat("r =", round(correlation$estimate, 3),
        ", p-value =", round(correlation$p.value, 5), "\n")
  } else {
    cat("\nNot enough common samples for correlation analysis\n")
  }
}

# ---------- Multiple Projects ----------

# Create DromaSet objects for multiple projects
ccle <- createDromaSetFromDatabase("ccle")
gdsc <- createDromaSetFromDatabase("gdsc")

# Load the same molecular data from each project
projects <- list(gCSI = gCSI, CCLE = ccle, GDSC = gdsc)

# For each project, load TP53 mutation data
for (proj_name in names(projects)) {
  tryCatch({
    projects[[proj_name]] <- loadMolecularProfiles(
      projects[[proj_name]],
      molecular_type = "mutation_gene",
      features = "TP53"
    )
    cat("\nLoaded TP53 mutation data for", proj_name, "\n")
  }, error = function(e) {
    cat("\nCould not load TP53 mutation data for", proj_name, ":", e$message, "\n")
  })
}

# ---------- Clean Up ----------

# Close the database connection when done
closeDROMADatabase()

cat("\nExample completed successfully\n")
