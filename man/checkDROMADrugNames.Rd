% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DROMA_SQLManager.R
\name{checkDROMADrugNames}
\alias{checkDROMADrugNames}
\title{Check and Harmonize Drug Names Against DROMA Database}
\usage{
checkDROMADrugNames(
  drug_names,
  connection = NULL,
  max_distance = 0.2,
  min_name_length = 5,
  keep_long_names_threshold = 17
)
}
\arguments{
\item{drug_names}{Character vector of drug names to check and harmonize}

\item{connection}{Optional database connection object. If NULL, uses global connection}

\item{max_distance}{Numeric, maximum distance for fuzzy matching (default: 0.2)}

\item{min_name_length}{Integer, minimum name length for partial matching (default: 5)}

\item{keep_long_names_threshold}{Integer, names longer than this will be kept as original (default: 17)}
}
\value{
A data frame with columns: original_name, cleaned_name, harmonized_name, match_type, match_confidence, new_name
}
\description{
Checks drug names (row names) against the drug_anno table in the DROMA database
and provides harmonized mappings. Uses fuzzy matching and name cleaning approach.
}
\examples{
\dontrun{
# Connect to database
con <- connectDROMADatabase("path/to/droma.sqlite")

# Check drug names from a drug response matrix
drug_names <- rownames(my_drug_matrix)
name_mapping <- checkDROMADrugNames(drug_names)

# View the mapping results
print(name_mapping)
}
}
