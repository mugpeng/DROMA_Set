# DROMA_Set: Drug Response and Omics Multi-project Analysis Set

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: mpl-2-0](https://img.shields.io/badge/MPL-2.0-yellow.svg)](https://opensource.org/licenses/mpl-2-0)

## Overview

**DROMA_Set** is a comprehensive R package for managing and analyzing drug response and omics data across multiple projects. It provides a robust framework for handling complex multi-omics datasets with integrated drug sensitivity information, enabling seamless cross-project comparisons and analyses.

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

The package requires the following R packages:
- `DBI` (>= 1.1.0)
- `RSQLite` (>= 2.2.0)
- `methods`

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

# Create a MultiDromaSet for multiple projects
multi_set <- createMultiDromaSetFromDatabase(
    project_names = c("gCSI", "CCLE"),
    db_path = "path/to/droma.sqlite"
)
```

### 4. Load and Analyze Data

```r
# Load molecular profiles
gCSI <- loadMolecularProfiles(gCSI, molecular_type = "mRNA", 
                             features = c("BRCA1", "BRCA2", "TP53"))

# Load treatment response data
gCSI <- loadTreatmentResponse(gCSI, drugs = c("Tamoxifen", "Cisplatin"))

# Cross-project molecular analysis
mRNA_data <- loadMultiProjectMolecularProfiles(multi_set, 
                                              molecular_type = "mRNA",
                                              overlap_only = FLASE)

# Cross-project treatment response analysis
drug_data <- loadMultiProjectTreatmentResponse(multi_set,
                                              drugs = c("Tamoxifen", "Cisplatin"),
                                              overlap_only = FLASE)
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
- `loadMolecularProfiles()`: Load omics data (mRNA, CNV, mutations, etc.)
- `loadTreatmentResponse()`: Load drug sensitivity data
- `availableMolecularProfiles()`: List available molecular data types
- `availableTreatmentResponses()`: List available treatment response types

### MultiDromaSet Class

The `MultiDromaSet` class manages multiple projects for cross-project analysis:

```r
# Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "database.sqlite")

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
- `loadMultiProjectMolecularProfiles()`: Load molecular data across multiple projects
- `loadMultiProjectTreatmentResponse()`: Load treatment response data across multiple projects
- `getDromaSet()`: Extract individual DromaSet from MultiDromaSet
- `availableProjects()`: List available projects

## Advanced Features

### 1. Load All Molecular Profiles

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
                                                        overlap_only = FLASE)
```

### 3. Database Management

```r
# Connect to database
connectDROMADatabase("droma.sqlite")

# Add new data to database
updateDROMADatabase(expression_matrix, "new_project_mRNA")

# List all tables
tables <- listDROMADatabaseTables()

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
                                              overlap_only = FLASE,
                                              data_type = "CellLine")

# 4. Load drug response data for overlapping samples
drug_data <- loadMultiProjectTreatmentResponse(multi_set,
                                              drugs = c("Tamoxifen", "Cisplatin"),
                                              overlap_only = FLASE,
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
- **drug_raw**: Raw drug response measurements

## Database Structure

The DROMA database uses a standardized table naming convention:
- `{project}_{datatype}`: Data tables (e.g., `gCSI_mRNA`, `CCLE_drug`)
- `sample_anno`: Sample metadata with ProjectID tracking
- `drug_anno`: Drug/treatment metadata with ProjectID tracking
- `projects`: Project summary information

## Examples

Comprehensive examples are provided in the `examples/` directory:

- `examples/produce_dromaset.R`: Basic DromaSet usage
- `examples/produce_multidromaset.R`: MultiDromaSet cross-project analysis
- `examples/produce_droma_database.R`: Database creation and management

## Performance Tips

1. You may **Use `overlap_only = TRUE`** when loading cross-project data to focus on same samples
2. **Specify `features` parameter** to load only genes/drugs of interest
3. **Use `return_data = TRUE`** when you only need the data without updating the object
4. **Filter by `data_type` and `tumor_type`** to reduce data loading time and focus on specific sample types
5. **Load molecular profiles incrementally** rather than using `molecular_type = "all"` for large datasets

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
