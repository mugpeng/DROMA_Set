% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiDromaSet-database.R
\name{createMultiDromaSetFromDatabase}
\alias{createMultiDromaSetFromDatabase}
\title{Create a MultiDromaSet from Database}
\usage{
createMultiDromaSetFromDatabase(
  project_names,
  db_path = file.path(path.expand("~"), "droma.sqlite"),
  db_groups = NULL,
  load_metadata = TRUE,
  dataset_types = NULL,
  auto_load = FALSE
)
}
\arguments{
\item{project_names}{Character vector, the names of the projects/datasets to include (e.g., c("gCSI", "CCLE"))}

\item{db_path}{Path to the SQLite database}

\item{db_groups}{Optional character vector of group names in the database, if different from project_names}

\item{load_metadata}{Logical, whether to load sample and treatment metadata (default: TRUE)}

\item{dataset_types}{Optional character vector of dataset types for each project (e.g., c("CellLine", "PDX"))}

\item{auto_load}{Logical, whether to automatically load treatment response and molecular profiles (default: FALSE)}
}
\value{
A MultiDromaSet object linked to the database
}
\description{
Creates a MultiDromaSet object containing multiple projects from a SQLite database
}
\examples{
\dontrun{
# Create a MultiDromaSet for gCSI and CCLE data from database
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "~/droma.sqlite")

# Create a MultiDromaSet and automatically load all data
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "~/droma.sqlite", auto_load = TRUE)

# Create with specific dataset types
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "PDX_data"),
  db_path = "~/droma.sqlite",
  dataset_types = c("CellLine", "PDX")
)
}
}
