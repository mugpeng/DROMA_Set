#' Create a DromaSet from Database
#'
#' @description Creates a DromaSet object linked to data in a SQLite database
#' @param project_name The name of the project/dataset (e.g., "gCSI", "CCLE")
#' @param db_path Path to the SQLite database
#' @param db_group Optional group name in the database, if different from project_name
#' @param load_metadata Logical, whether to load sample and treatment metadata (default: TRUE)
#' @param dataset_type Optional dataset type (e.g., "CellLine", "PDX", "PDO")
#' @param auto_load Logical, whether to automatically load treatment response and molecular profiles (default: FALSE)
#' @return A DromaSet object linked to the database
#' @export
#' @examples
#' \dontrun{
#' # Create a DromaSet for gCSI data from database
#' gCSI <- createDromaSetFromDatabase("gCSI", "~/droma.sqlite")
#'
#' # Create a DromaSet and automatically load drug response data
#' gCSI <- createDromaSetFromDatabase("gCSI", "~/droma.sqlite", auto_load = TRUE)
#' }
createDromaSetFromDatabase <- function(project_name, db_path = file.path(path.expand("~"), "droma.sqlite"),
                                    db_group = NULL, load_metadata = TRUE, dataset_type = NULL, auto_load = FALSE) {

  if (!file.exists(db_path)) {
    stop("Database file not found: ", db_path)
  }

  # Set db_group to project_name if not specified
  if (is.null(db_group)) {
    db_group <- project_name
  }

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if tables exist for this project
  all_tables <- DBI::dbListTables(con)
  project_tables <- grep(paste0("^", db_group, "_"), all_tables, value = TRUE)

  if (length(project_tables) == 0) {
    stop("No tables found for project '", project_name, "' with group prefix '", db_group, "'")
  }

  # Load metadata if requested
  sample_metadata <- data.frame()
  treatment_metadata <- data.frame()

  if (load_metadata) {
    # Try to load sample metadata
    if ("sample_anno" %in% all_tables) {
      # Get sample IDs from project tables - handle both continuous and discrete data
      sample_ids <- character()

      for (table in project_tables) {
        tryCatch({
          # Check if this is a discrete data table (mutation, fusion)
          table_info <- DBI::dbGetQuery(con, paste0("SELECT name FROM pragma_table_info('", table, "')"))
          column_names <- table_info$name

          if ("samples" %in% column_names || "cells" %in% column_names) {
            # Discrete data - samples are in a column
            sample_col <- ifelse("samples" %in% column_names, "samples", "cells")
            sample_query_discrete <- paste0("SELECT DISTINCT ", sample_col, " FROM ", table)
            discrete_samples <- DBI::dbGetQuery(con, sample_query_discrete)[[sample_col]]
            sample_ids <- c(sample_ids, discrete_samples)
          } else {
            # Continuous data - samples are column names (excluding feature_id)
            continuous_samples <- setdiff(column_names, "feature_id")
            sample_ids <- c(sample_ids, continuous_samples)
          }
        }, error = function(e) {
          warning("Problem with getting sample IDs from table '", table, "': ", e$message)
        })
      }

      # Get unique sample IDs
      sample_ids <- unique(sample_ids)

      if (length(sample_ids) > 0) {
        # Construct query to get sample metadata for these samples
        sample_ids_str <- paste0("'", sample_ids, "'", collapse = ",")
        sample_query <- paste0("SELECT * FROM sample_anno WHERE SampleID IN (", sample_ids_str, ")")

        tryCatch({
          sample_metadata <- DBI::dbGetQuery(con, sample_query)
          # Filter by project if ProjectID column exists
          if ("ProjectID" %in% colnames(sample_metadata)) {
            sample_metadata <- sample_metadata[sample_metadata$ProjectID %in% project_name, ]
          }
          sample_metadata <- unique(sample_metadata)
        }, error = function(e) {
          warning("Problem with loading sample metadata: ", e$message)
        })
      }
    }

    # Try to load treatment metadata
    if ("drug_anno" %in% all_tables) {
      # Construct query to get drugs for this project
      drug_query <- paste0(
        "SELECT * FROM drug_anno WHERE DrugName IN (",
        "SELECT feature_id FROM ", db_group, "_drug)"
      )

      tryCatch({
        treatment_metadata <- DBI::dbGetQuery(con, drug_query)
        treatment_metadata <- treatment_metadata[treatment_metadata$ProjectID %in% project_name,]
        treatment_metadata <- unique(treatment_metadata)
      }, error = function(e) {
        warning("Problem with loading treatment metadata: ", e$message)
      })
    }
  }

  # If dataset_type is not provided, try to infer from sample metadata
  if (is.null(dataset_type) && nrow(sample_metadata) > 0 && "DataType" %in% colnames(sample_metadata)) {
    dataset_type <- unique(sample_metadata$DataType)[1]
  }

  # Create the DromaSet object
  object <- DromaSet(
    name = project_name,
    sampleMetadata = sample_metadata,
    treatmentMetadata = treatment_metadata,
    datasetType = ifelse(is.null(dataset_type), NA_character_, dataset_type),
    db_info = list(
      db_path = db_path,
      db_group = db_group
    )
  )

  # Auto-load treatment response and molecular profiles if requested
  if (auto_load) {
    # Load treatment response data if available
    if (paste0(db_group, "_drug") %in% all_tables) {
      tryCatch({
        object <- loadTreatmentResponse(object)
        message("Loaded treatment response data for project '", project_name, "'")
      }, error = function(e) {
        warning("Problem with loading treatment response data: ", e$message)
      })
    }

    # Get molecular profile types available for this project
    profile_types <- unique(sapply(project_tables, function(t) {
      sub(paste0("^", db_group, "_"), "", t)
    }))
    # Remove drug table from profile types
    profile_types <- setdiff(profile_types, "drug")

    if (length(profile_types) > 0) {
      for (profile_type in profile_types) {
        tryCatch({
          object <- loadMolecularProfiles(object, molecular_type = profile_type)
          message("Loaded '", profile_type, "' molecular profile data for project '", project_name, "'")
        }, error = function(e) {
          warning("Problem with loading '", profile_type, "' molecular profile data: ", e$message)
        })
      }
    }
  }

  return(object)
}
