% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DromaSet-class.R
\name{loadTreatmentResponse}
\alias{loadTreatmentResponse}
\alias{loadTreatmentResponse,DromaSet-method}
\title{Load Treatment Response from Database}
\usage{
loadTreatmentResponse(
  object,
  drugs = NULL,
  samples = NULL,
  return_data = FALSE,
  data_type = "all",
  tumor_type = "all"
)

\S4method{loadTreatmentResponse}{DromaSet}(
  object,
  drugs = NULL,
  samples = NULL,
  return_data = FALSE,
  data_type = "all",
  tumor_type = "all"
)
}
\arguments{
\item{object}{A DromaSet object}

\item{drugs}{Optional vector of drug names to load. If NULL, loads all drugs.}

\item{samples}{Optional vector of sample IDs to load. If NULL, loads all samples.}

\item{return_data}{Logical, if TRUE returns the loaded data directly instead of updating the object (default: FALSE)}

\item{data_type}{Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"}

\item{tumor_type}{Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")}
}
\value{
Updated DromaSet object with loaded drug response data or the loaded data directly if return_data=TRUE
}
\description{
Loads drug response data from the database into a DromaSet object
}
