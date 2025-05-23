# DROMASet: Drug Response and Omics Multi-project Analysis Set

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**DROMASet** is a comprehensive R package for managing and analyzing drug response and omics data across multiple projects. It provides a robust framework for handling complex multi-omics datasets with integrated drug sensitivity information, enabling seamless cross-project comparisons and analyses.

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

# Install DROMASet
devtools::install_github("mugpeng/DROMASet")
```

### Dependencies

The package requires the following R packages:
- `DBI` (>= 1.1.0)
- `RSQLite` (>= 2.2.0)
- `methods`

These will be automatically installed when you install DROMASet.

## Quick Start

### 1. Load the Package

```r
library(DROMASet)
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

# Cross-project analysis
mRNA_data <- loadMultiProjectData(multi_set, 
                                 data_type = "molecular",
                                 molecular_type = "mRNA",
                                 overlap_only = TRUE)
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

# Load data across projects
cross_data <- loadMultiProjectData(multi_set, 
                                  data_type = "molecular",
                                  molecular_type = "mRNA")
```

**Key Methods:**
- `getOverlappingSamples()`: Identify samples present in multiple projects
- `loadMultiProjectData()`: Load data across multiple projects
- `getDromaSet()`: Extract individual DromaSet from MultiDromaSet
- `availableProjects()`: List available projects

## Advanced Features

### 1. Load All Molecular Profiles

```r
# Load all available molecular profile types
all_data <- loadMolecularProfiles(dataset, molecular_type = "all")

# Cross-project loading of all molecular types
all_cross_data <- loadMultiProjectData(multi_set,
                                      data_type = "molecular", 
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
mRNA_data <- loadMultiProjectData(multi_set,
                                 data_type = "molecular",
                                 molecular_type = "mRNA",
                                 features = c("BRCA1", "BRCA2"),
                                 overlap_only = TRUE)

# 4. Load drug response data
drug_data <- loadMultiProjectData(multi_set,
                                 data_type = "treatment",
                                 drugs = c("Tamoxifen", "Cisplatin"),
                                 overlap_only = TRUE)

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

1. **Use `overlap_only = TRUE`** when loading cross-project data to focus on comparable samples
2. **Specify `features` parameter** to load only genes/drugs of interest
3. **Use `return_data = TRUE`** when you only need the data without updating the object
4. **Filter by `data_type` and `tumor_type`** to reduce data loading time
5. **Load molecular profiles incrementally** rather than using `molecular_type = "all"` for large datasets

## Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use DROMASet in your research, please cite:

```
Zhong, P.Y. (2024). DROMASet: Drug Response and Omics Multi-project Analysis Set. 
R package version 0.9.0. https://github.com/mugpeng/DROMASet
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- 📧 **Email**: yc47680@um.edu.mo
- 🐛 **Issues**: [GitHub Issues](https://github.com/mugpeng/DROMASet/issues)
- 📖 **Documentation**: [Package Documentation](https://mugpeng.github.io/DROMASet/)

## Changelog

### Version 0.9.0
- Initial release
- DromaSet and MultiDromaSet classes
- Database integration and management
- Cross-project analysis capabilities
- Comprehensive molecular profile support
- Sample overlap detection and analysis
- Enhanced metadata management with ProjectID tracking
- Support for loading all molecular profile types with `molecular_type = "all"`

---

**DROMASet** - Empowering multi-project drug response and omics analysis 🧬💊
