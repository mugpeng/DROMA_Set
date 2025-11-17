# DROMA_Set: Drug Response and Omics Multi-project Analysis Set

[![Website](https://img.shields.io/website?url=https%3A//droma01.github.io/DROMA/)](https://droma01.github.io/DROMA/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: mpl-2-0](https://img.shields.io/badge/MPL-2.0-yellow.svg)](https://opensource.org/licenses/mpl-2-0)

## Overview

**DROMA_Set** is a comprehensive R package for managing and analyzing drug response and omics data across multiple projects. It provides a robust framework for handling complex multi-omics datasets with integrated drug sensitivity information, enabling seamless cross-project comparisons and analyses.

It is a part of [DROMA project](https://github.com/mugpeng/DROMA). Visit the [official DROMA website](https://droma01.github.io/DROMA/) for comprehensive documentation and interactive examples.

### Key Features

- **🔬 Multi-omics Data Management**: Support for various molecular profile types (mRNA, CNV, mutations, methylation, proteomics)
- **💊 Drug Response Integration**: Comprehensive treatment response data handling and analysis
- **🔗 Cross-Project Analysis**: Advanced tools for comparing and analyzing data across multiple projects
- **📊 Sample Overlap Detection**: Automatic identification and analysis of overlapping samples between projects
- **🗄️ Database Integration**: Robust SQLite database connectivity with efficient data storage and retrieval
- **📈 Flexible Data Loading**: Smart data loading with filtering by data type, tumor type, and specific features
- **🎯 Metadata Management**: Comprehensive sample and treatment metadata handling with ProjectID tracking

## Installation

### From GitHub (Recommended)

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install DROMA_Set
devtools::install_github("mugpeng/DROMA_Set")
```

### Dependencies

**Required packages:**
- `DBI` (>= 1.1.0)
- `RSQLite` (>= 2.2.0)
- `methods`

**Suggested packages for enhanced functionality:**
- `data.table`: For efficient large dataset processing
- `parallel`: For parallel processing of multiple molecular types (Unix/Linux/macOS)

These will be automatically installed when you install DROMA_Set.

## Quick Start

### 1. Load the Package

```r
library(DROMA.Set)
```

### 2. Connect to Database

```r
# Connect to your DROMA database
connectDROMADatabase("path/to/your/droma.sqlite")

# List available projects
projects <- listDROMAProjects()
print(projects)
```

### 3. Create DromaSet Objects

```r
# Create a single DromaSet for one project
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Create a DromaSet with automatic data loading
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite", auto_load = TRUE)

# Create a MultiDromaSet for multiple projects
multi_set <- createMultiDromaSetFromDatabase(
    project_names = c("gCSI", "CCLE"),
    db_path = "path/to/droma.sqlite"
)

# Create MultiDromaSet with specific dataset types
multi_set <- createMultiDromaSetFromDatabase(
    project_names = c("gCSI", "PDX_data"),
    db_path = "path/to/droma.sqlite",
    dataset_types = c("CellLine", "PDX")
)
```

### 4. Load and Analyze Data

```r
# Load molecular profiles
gCSI <- loadMolecularProfiles(gCSI, molecular_type = "mRNA", 
                             features = c("BRCA1", "BRCA2", "TP53"))

# Load molecular profiles with advanced filtering
gCSI <- loadMolecularProfiles(gCSI, molecular_type = "mRNA",
                             data_type = "CellLine", 
                             tumor_type = "breast cancer",
                             chunk_size = 100000,
                             validate_features = TRUE)

# Load treatment response data
gCSI <- loadTreatmentResponse(gCSI, drugs = c("Tamoxifen", "Cisplatin"))

# Load treatment response with filtering
gCSI <- loadTreatmentResponse(gCSI, drugs = c("Tamoxifen", "Cisplatin"),
                             data_type = "CellLine", 
                             tumor_type = "breast cancer")

# Cross-project molecular analysis
mRNA_data <- loadMultiProjectMolecularProfiles(multi_set, 
                                              molecular_type = "mRNA",
                                              overlap_only = FALSE)

# Cross-project treatment response analysis
drug_data <- loadMultiProjectTreatmentResponse(multi_set,
                                              drugs = c("Tamoxifen", "Cisplatin"),
                                              overlap_only = FALSE)
```

## Core Classes

### DromaSet Class

The `DromaSet` class represents a single project's drug response and omics data:

```r
# Create DromaSet
dataset <- createDromaSetFromDatabase("project_name", "database.sqlite")

# Load all molecular profiles
dataset <- loadMolecularProfiles(dataset, molecular_type = "all")

# Check available data types
availableMolecularProfiles(dataset)
availableTreatmentResponses(dataset)
```

**Key Methods:**
- `loadMolecularProfiles()`: Load omics data (mRNA, CNV, mutations, etc.) with advanced filtering options
- `loadTreatmentResponse()`: Load drug sensitivity data with filtering by data type and tumor type
- `availableMolecularProfiles()`: List available molecular data types
- `availableTreatmentResponses()`: List available treatment response types

### MultiDromaSet Class

The `MultiDromaSet` class manages multiple projects for cross-project analysis:

```r
# Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "database.sqlite")

# Create from existing DromaSet objects
multi_set <- createMultiDromaSetFromObjects(gCSI, CCLE)

# Add new DromaSet to existing MultiDromaSet
multi_set <- addDromaSetToMulti(multi_set, new_dromaset)

# Remove DromaSet from MultiDromaSet
multi_set <- removeDromaSetFromMulti(multi_set, "CCLE")

# Create subset of MultiDromaSet
subset_multi <- subset(multi_set, projects = c("gCSI"))

# Find overlapping samples
overlap_info <- getOverlappingSamples(multi_set)

# Load molecular data across projects
mRNA_data <- loadMultiProjectMolecularProfiles(multi_set, 
                                              molecular_type = "mRNA")

# Load treatment response data across projects
drug_data <- loadMultiProjectTreatmentResponse(multi_set,
                                              drugs = c("Tamoxifen", "Cisplatin"))
```

**Key Methods:**
- `getOverlappingSamples()`: Identify samples present in multiple projects
- `loadMultiProjectMolecularProfiles()`: Load molecular data across multiple projects with filtering
- `loadMultiProjectTreatmentResponse()`: Load treatment response data across multiple projects with filtering
- `getDromaSet()`: Extract individual DromaSet from MultiDromaSet
- `availableProjects()`: List available projects
- `createMultiDromaSetFromObjects()`: Create from existing DromaSet objects
- `addDromaSetToMulti()`: Add new DromaSet to existing MultiDromaSet
- `removeDromaSetFromMulti()`: Remove DromaSet from MultiDromaSet
- `subset()`: Create subset with specific projects

## Advanced Features

### 1. Advanced Database Operations

```r
# Check and harmonize sample names
sample_mapping <- checkDROMASampleNames(colnames(my_data))

# Update sample annotations with harmonized names
updateDROMAAnnotation("sample", sample_mapping, project_name = "MyProject",
                     data_type = "CellLine", tumor_type = "breast cancer")

# Check and harmonize drug names  
drug_mapping <- checkDROMADrugNames(rownames(my_drug_data))

# Update drug annotations
updateDROMAAnnotation("drug", drug_mapping, project_name = "MyProject")

# Get feature data with advanced filtering
feature_data <- getFeatureFromDatabase("mRNA", "BRCA1", 
                                      data_sources = c("gCSI", "CCLE"),
                                      data_type = "CellLine")

# Create MultiDromaSet from all available projects
multi_all <- createMultiDromaSetFromAllProjects("droma.sqlite",
                                               exclude_projects = "test_data")
```

### 2. Matrix Database Storage (FuncMatrixDatabase.R)

```r
# Store matrix data directly in database
storeMatricesInDatabase("my_database.sqlite", expression_matrix, "experiment1_mRNA")

# Retrieve matrix data
retrieved_matrix <- retrieveMatrixFromDatabase("my_database.sqlite", "experiment1_mRNA")

# Retrieve specific features only
subset_matrix <- retrieveMatrixFromDatabase("my_database.sqlite", "experiment1_mRNA",
                                          features = c("BRCA1", "TP53", "EGFR"))

# List all matrix tables in database
matrix_tables <- listMatrixTables("my_database.sqlite")
```

### 3. Load All Molecular Profiles

```r
# Load all available molecular profile types
all_data <- loadMolecularProfiles(dataset, molecular_type = "all")

# Cross-project loading of all molecular types
all_cross_data <- loadMultiProjectMolecularProfiles(multi_set,
                                                   molecular_type = "all")
```

### 2. Sample and Data Filtering

```r
# Filter by data type and tumor type
filtered_data <- loadMolecularProfiles(dataset,
                                      molecular_type = "mRNA",
                                      data_type = "CellLine",
                                      tumor_type = "breast cancer")

# Load specific features and samples
specific_data <- loadMolecularProfiles(dataset,
                                      molecular_type = "mRNA",
                                      features = c("BRCA1", "TP53"),
                                      samples = c("sample1", "sample2"))

# Cross-project filtering by data type and tumor type
filtered_cross_data <- loadMultiProjectMolecularProfiles(multi_set,
                                                        molecular_type = "mRNA",
                                                        data_type = "CellLine",
                                                        tumor_type = "breast cancer",
                                                        overlap_only = FALSE)
```

### 3. Database Management

```r
# Connect to database
connectDROMADatabase("droma.sqlite")

# Add new data to database
updateDROMADatabase(expression_matrix, "new_project_mRNA")

# List all tables with metadata
tables <- listDROMADatabaseTables()

# List available projects
projects <- listDROMAProjects()

# Update project metadata
updateDROMAProjects("gCSI", dataset_type = "CellLine")

# List features for a specific project and data type
features <- listDROMAFeatures("gCSI", "mRNA", limit = 100)

# List samples for a project
samples <- listDROMASamples("gCSI", data_type = "CellLine")

# Get annotation data
sample_anno <- getDROMAAnnotation("sample", project_name = "gCSI")
drug_anno <- getDROMAAnnotation("drug", project_name = "gCSI")

# Close connection
closeDROMADatabase()
```

### 4. Cross-Project Analysis Workflow

```r
# 1. Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))

# 2. Find overlapping samples
overlaps <- getOverlappingSamples(multi_set)
cat("Found", overlaps$overlap_count, "overlapping samples")

# 3. Load molecular data for overlapping samples
mRNA_data <- loadMultiProjectMolecularProfiles(multi_set,
                                              molecular_type = "mRNA",
                                              features = c("BRCA1", "BRCA2"),
                                              overlap_only = FALSE,
                                              data_type = "CellLine")

# 4. Load drug response data for overlapping samples
drug_data <- loadMultiProjectTreatmentResponse(multi_set,
                                              drugs = c("Tamoxifen", "Cisplatin"),
                                              overlap_only = FALSE,
                                              data_type = "CellLine")

# 5. Perform correlation analysis
for (project in names(mRNA_data)) {
    if (project %in% names(drug_data)) {
        # Analyze correlations between gene expression and drug response
        # Your analysis code here
    }
}
```

## Data Types Supported

### Molecular Profiles
- **mRNA**: Gene expression data
- **cnv**: Copy number variation data
- **mutation_gene**: Gene-level mutation data
- **mutation_site**: Site-specific mutation data
- **fusion**: Gene fusion data
- **meth**: DNA methylation data
- **proteinrppa**: Reverse-phase protein array data
- **proteinms**: Mass spectrometry proteomics data

### Treatment Response
- **drug**: Drug sensitivity/response data

## Database Structure

The DROMA database uses a standardized table naming convention:
- `{project}_{datatype}`: Data tables (e.g., `gCSI_mRNA`, `CCLE_drug`)
- `sample_anno`: Sample metadata with ProjectID tracking
- `drug_anno`: Drug/treatment metadata with ProjectID tracking
- `projects`: Project summary information

## Database Utility Functions

### Connection Management
- `connectDROMADatabase()`: Connect to DROMA database
- `closeDROMADatabase()`: Close database connection
- `connectCTRDatabase()`: Connect to Clinical Trial Response Database
- `closeCTRDatabase()`: Close CTRDB connection

### Data Management
- `updateDROMADatabase()`: Add/update data tables
- `updateDROMAProjects()`: Update project metadata
- `updateDROMAAnnotation()`: Update sample/drug annotations with harmonized names

### Query Functions
- `listDROMAProjects()`: List available projects
- `listDROMADatabaseTables()`: List all data tables with metadata
- `listDROMAFeatures()`: List features for specific project/data type
- `listDROMASamples()`: List samples with filtering options
- `getDROMAAnnotation()`: Get annotation data
- `getFeatureFromDatabase()`: Get feature data with complex filtering

### Name Harmonization
- `checkDROMASampleNames()`: Check and harmonize sample names
- `checkDROMADrugNames()`: Check and harmonize drug names

### Matrix Storage (FuncMatrixDatabase.R)
- `storeMatricesInDatabase()`: Store matrix data in SQLite database
- `retrieveMatrixFromDatabase()`: Retrieve matrix data from database
- `listMatrixTables()`: List matrix tables with metadata

## Examples

Comprehensive examples are provided in the `examples/` directory:

- `examples/produce_dromaset.R`: Basic DromaSet usage
- `examples/produce_multidromaset.R`: MultiDromaSet cross-project analysis  
- `examples/produce_droma_database.R`: Database creation and management

### Complete Workflow Example

```r
# 1. Connect to database
library(DROMA.Set)
con <- connectDROMADatabase("path/to/droma.sqlite")

# 2. List available projects and data types
projects <- listDROMAProjects()
print(projects)

# 3. Create DromaSet with automatic loading
gCSI <- createDromaSetFromDatabase("gCSI", auto_load = TRUE)

# 4. Load specific molecular profiles with filtering
gCSI <- loadMolecularProfiles(gCSI, 
                             molecular_type = "mRNA",
                             features = c("BRCA1", "BRCA2", "TP53"),
                             data_type = "CellLine",
                             tumor_type = "breast cancer")

# 5. Create MultiDromaSet for cross-project analysis
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))

# 6. Find overlapping samples
overlaps <- getOverlappingSamples(multi_set)
print(paste("Found", overlaps$overlap_count, "overlapping samples"))

# 7. Load cross-project data
cross_mRNA <- loadMultiProjectMolecularProfiles(multi_set,
                                               molecular_type = "mRNA",
                                               overlap_only = TRUE)

# 8. Clean up
closeDROMADatabase()
```

## Performance Tips

1. **Use `overlap_only = TRUE`** when loading cross-project data to focus on overlapping samples
2. **Specify `features` parameter** to load only genes/drugs of interest
3. **Use `return_data = TRUE`** when you only need the data without updating the object
4. **Filter by `data_type` and `tumor_type`** to reduce data loading time and focus on specific sample types
5. **Load molecular profiles incrementally** rather than using `molecular_type = "all"` for large datasets
6. **Use `chunk_size` parameter** for large datasets to optimize memory usage (default: 100,000 rows)
7. **Set `validate_features = FALSE`** to skip feature validation for faster loading when you're confident features exist
8. **Use parallel processing** - the package automatically uses parallel processing for loading multiple molecular types on Unix-like systems
9. **Leverage database indexing** - the package creates indexes on feature_id columns for faster queries
10. **Use `limit` parameter** in list functions to preview data before loading full datasets

## Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use DROMA_Set in your research, please cite:

```
Li, S., Peng, Y., Chen, M. et al. Facilitating integrative and personalized oncology omics analysis with UCSCXenaShiny. Commun Biol 7, 1200 (2024). https://doi.org/10.1038/s42003-024-06891-2
```

## License

This project is licensed under the MPL-2 License - see the [LICENSE](LICENSE) file for details.

## Support

- 📧 **Email**: yc47680@um.edu.mo
- 🐛 **Issues**: [GitHub Issues](https://github.com/mugpeng/DROMA_Set/issues)
- 📖 **Documentation**: [Package Documentation](https://mugpeng.github.io/DROMA_Set/)

## Changelog
### Version 0.4.6
Current stable release with comprehensive documentation updates and enhanced functionality.

### Version 0.4.4
Refactor updateDROMADatabase and updateDROMAProjects functions to improve project tracking and metadata handling; enhance listDROMADatabaseTables to filter out backup tables and include created/updated dates; update documentation for new parameters in updateDROMAAnnotation function to support vector inputs for age, data type, and other attributes.

Enhancements Made:
✅ Removed projects table auto-updates from updateDROMADatabase
✅ Added _mutation_raw table exclusion across all relevant functions
✅ Added dataset_type parameter to updateDROMAProjects
✅ Enhanced updateDROMAAnnotation with vector support and created_date logic
✅ Improved parameter validation and documentation

### Version 0.4.3
Add updateDROMAProjects function to manage project metadata in DROMA database; enhance listDROMADatabaseTables with feature and sample counts; minor adjustments in example script.

### Version 0.4.1
- Initial release
- DromaSet and MultiDromaSet classes
- Database integration and management
- Cross-project analysis capabilities
- Comprehensive molecular profile support
- Sample overlap detection and analysis
- Enhanced metadata management with ProjectID tracking
- Support for loading all molecular profile types with `molecular_type = "all"`
- Split cross-project data loading into specialized functions:
  - `loadMultiProjectMolecularProfiles()` for molecular data
  - `loadMultiProjectTreatmentResponse()` for treatment response data
- Added `data_type` and `tumor_type` filtering parameters for enhanced sample selection

---

**DROMA_Set** - Empowering multi-project drug response and omics analysis 🧬💊
