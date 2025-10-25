#!/usr/bin/env Rscript

#' Z-score normalization for omics data
#'
#' @description Apply Z-score normalization to omics data at the gene level
#' @param mat A matrix or data frame with genes in rows and samples in columns
#' @return A normalized matrix with Z-scores
#' @export
zscoreNormalize <- function(mat) {
  # Apply z-score normalization to each row (gene)
  normalized <- t(scale(t(mat)))
  # Handle any potential NaN values (e.g., if SD was 0)
  normalized[is.nan(normalized)] <- 0
  return(normalized)
}

#' DromaSet Class
#'
#' @description A class to represent a DROMA project dataset with drug response and omics data.
#' Each DromaSet contains treatment response data (drug sensitivity) and multiple molecular profiles
#' (omics data) for a specific project.
#'
#' @slot name Character string, the name of the dataset (e.g., "gCSI", "CCLE")
#' @slot treatmentResponse List containing drug response data matrices
#' @slot molecularProfiles List containing different types of molecular data (omics)
#' @slot sampleMetadata Data frame with sample annotations
#' @slot treatmentMetadata Data frame with drug annotations
#' @slot datasetType Character string, the type of dataset (e.g., "CellLine", "PDX", "PDO")
#' @slot db_info List containing database connection information
#' @export
setClass("DromaSet",
  slots = c(
    name = "character",
    treatmentResponse = "list",
    molecularProfiles = "list",
    sampleMetadata = "data.frame",
    treatmentMetadata = "data.frame",
    datasetType = "character",
    db_info = "list"
  ),
  prototype = list(
    name = NA_character_,
    treatmentResponse = list(),
    molecularProfiles = list(),
    sampleMetadata = data.frame(),
    treatmentMetadata = data.frame(),
    datasetType = NA_character_,
    db_info = list()
  )
)

#' Create a DromaSet Object
#'
#' @description Creates a DromaSet object to store drug response and omics data for a specific project
#' @param name Character string, the name of the dataset (e.g., "gCSI", "CCLE")
#' @param treatmentResponse List containing drug response data matrices
#' @param molecularProfiles List containing different types of molecular data (omics)
#' @param sampleMetadata Data frame with sample annotations
#' @param treatmentMetadata Data frame with drug annotations
#' @param datasetType Character string, the type of dataset (e.g., "CellLine", "PDX", "PDO")
#' @param db_info List containing database connection information
#' @return A DromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Create a DromaSet object for gCSI data
#' gCSI <- DromaSet(name = "gCSI",
#'                  datasetType = "CellLine",
#'                  db_info = list(db_path = "~/droma.sqlite", db_group = "gCSI"))
#' }
DromaSet <- function(name,
                   treatmentResponse = list(),
                   molecularProfiles = list(),
                   sampleMetadata = data.frame(),
                   treatmentMetadata = data.frame(),
                   datasetType = NA_character_,
                   db_info = list()) {

  # Create new DromaSet object
  object <- new("DromaSet",
              name = name,
              treatmentResponse = treatmentResponse,
              molecularProfiles = molecularProfiles,
              sampleMetadata = sampleMetadata,
              treatmentMetadata = treatmentMetadata,
              datasetType = datasetType,
              db_info = db_info)

  # Return the object
  return(object)
}

#' Show Method for DromaSet objects
#'
#' @description Displays information about a DromaSet object
#' @param object A DromaSet object
#' @return NULL, prints information to console
#' @export
setMethod("show", "DromaSet", function(object) {
  cat("DromaSet Object:", object@name, "\n")
  cat("Dataset Type:", object@datasetType, "\n\n")

  # Treatment Response
  cat("Treatment Response Data:\n")
  if (length(object@treatmentResponse) > 0) {
    tr_info <- lapply(object@treatmentResponse, function(x) {
      if(is.matrix(x) || is.data.frame(x)) {
        paste0("  (", nrow(x), " drugs x ", ncol(x), " samples)")
      } else {
        "  (data in database)"
      }
    })
    cat(paste(names(object@treatmentResponse), tr_info, sep = ": ", collapse = "\n"), "\n")
  } else {
    cat("  None loaded (may be available in database)\n")
  }

  # Molecular Profiles
  cat("\nMolecular Profiles:\n")
  if (length(object@molecularProfiles) > 0) {
    mp_info <- lapply(object@molecularProfiles, function(x) {
      if(is.matrix(x) || is.data.frame(x)) {
        paste0("  (", nrow(x), " features x ", ncol(x), " samples)")
      } else {
        "  (data in database)"
      }
    })
    cat(paste(names(object@molecularProfiles), mp_info, sep = ": ", collapse = "\n"), "\n")
  } else {
    cat("  None loaded (may be available in database)\n")
  }

  # Database information
  cat("\nDatabase Connection Information:\n")
  if (length(object@db_info) > 0) {
    cat("  Path:", ifelse(is.null(object@db_info$db_path), "Not specified", object@db_info$db_path), "\n")
    cat("  Group:", ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group), "\n")
  } else {
    cat("  No database information available\n")
  }

  invisible(NULL)
})

#' Get Available Molecular Profile Types
#'
#' @description Returns the types of molecular profiles available for a DromaSet
#' @param object A DromaSet object
#' @param include_db Logical, whether to include profiles available in the database but not loaded (default: TRUE)
#' @return Character vector of available molecular profile types
#' @export
setGeneric("availableMolecularProfiles", function(object, include_db = TRUE) standardGeneric("availableMolecularProfiles"))

#' @rdname availableMolecularProfiles
#' @export
setMethod("availableMolecularProfiles", "DromaSet", function(object, include_db = TRUE) {
  # Get profiles already loaded in the object
  loaded_profiles <- names(object@molecularProfiles)

  # If requested, check database for additional profiles
  if (include_db && length(object@db_info) > 0 && !is.null(object@db_info$db_path)) {
    if (file.exists(object@db_info$db_path)) {
      # Connect to database
      con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      # Get group prefix (if specified) or use dataset name
      group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

      # List all tables with this prefix
      all_tables <- DBI::dbListTables(con)
      project_tables <- grep(paste0("^", group_prefix, "_"), all_tables, value = TRUE)

      # Extract molecular profile types from table names
      if (length(project_tables) > 0) {
        # Remove prefix to get profile types
        db_profiles <- unique(sub(paste0("^", group_prefix, "_"), "", project_tables))
        # Remove 'drug' and 'drug_raw' which are treatment responses
        db_profiles <- db_profiles[!db_profiles %in% c("drug", "drug_raw")]

        # Combine with loaded profiles
        return(unique(c(loaded_profiles, db_profiles)))
      }
    }
  }

  return(loaded_profiles)
})

#' Get Available Treatment Response Types
#'
#' @description Returns the types of treatment response data available for a DromaSet
#' @param object A DromaSet object
#' @param include_db Logical, whether to include data available in the database but not loaded (default: TRUE)
#' @return Character vector of available treatment response types
#' @export
setGeneric("availableTreatmentResponses", function(object, include_db = TRUE) standardGeneric("availableTreatmentResponses"))

#' @rdname availableTreatmentResponses
#' @export
setMethod("availableTreatmentResponses", "DromaSet", function(object, include_db = TRUE) {
  # Get responses already loaded in the object
  loaded_responses <- names(object@treatmentResponse)

  # If requested, check database for additional responses
  if (include_db && length(object@db_info) > 0 && !is.null(object@db_info$db_path)) {
    if (file.exists(object@db_info$db_path)) {
      # Connect to database
      con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      # Get group prefix (if specified) or use dataset name
      group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

      # Look for drug and drug_raw tables
      all_tables <- DBI::dbListTables(con)
      drug_tables <- grep(paste0("^", group_prefix, "_(drug|drug_raw)$"), all_tables, value = TRUE)

      # Extract response types from table names
      if (length(drug_tables) > 0) {
        # Remove prefix to get response types
        db_responses <- unique(sub(paste0("^", group_prefix, "_"), "", drug_tables))

        # Combine with loaded responses
        return(unique(c(loaded_responses, db_responses)))
      }
    }
  }

  return(loaded_responses)
})

#' Load Molecular Profiles from Database
#'
#' @description Loads specific molecular profile data from the database into a DromaSet object
#' @param object A DromaSet object
#' @param molecular_type The type of molecular data to load (e.g., "mRNA", "cnv", "mutation_gene") or "all" to load all available types
#' @param features Optional vector of feature names to load. If NULL, loads all features.
#' @param samples Optional vector of sample IDs to load. If NULL, loads all samples.
#' @param return_data Logical, if TRUE returns the loaded data directly instead of updating the object (default: FALSE)
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")
#' @param chunk_size Integer, number of rows to process at a time for large datasets (default: 100000)
#' @param validate_features Logical, whether to validate that specified features exist in the database (default: TRUE)
#' @param zscore Logical, whether to apply z-score normalization to continuous data (default: FALSE)
#' @param format Character, format of returned discrete data: "wide" (features as rows, samples as columns) or "long" (default, original database format with features and samples columns)
#' @return Updated DromaSet object with loaded molecular data or the loaded data directly if return_data=TRUE
#' @export
setGeneric("loadMolecularProfiles", function(object, molecular_type, features = NULL, samples = NULL, return_data = FALSE, data_type = "all", tumor_type = "all", chunk_size = 100000, validate_features = TRUE, zscore = FALSE, format = "long")
  standardGeneric("loadMolecularProfiles"))

#' @rdname loadMolecularProfiles
#' @export
setMethod("loadMolecularProfiles", "DromaSet", function(object, molecular_type, features = NULL, samples = NULL, return_data = FALSE, data_type = "all", tumor_type = "all", chunk_size = 100000, validate_features = TRUE, zscore = FALSE, format = "long") {
  # Handle "all" molecular_type option with parallel processing
  if (molecular_type == "all") {
    # Get all available molecular profile types
    available_types <- availableMolecularProfiles(object, include_db = TRUE)

    if (length(available_types) == 0) {
      warning("No molecular profile types available for dataset '", object@name, "'")
      if (return_data) {
        return(list())
      } else {
        return(object)
      }
    }

    # Load molecular profile types in parallel if available
    if (requireNamespace("parallel", quietly = TRUE) && length(available_types) > 1 && Sys.info()["sysname"] != "Windows") {
      # Use mclapply for Unix-like systems
      all_data <- parallel::mclapply(available_types, function(mol_type) {
        tryCatch({
          loadMolecularProfiles(
            object = object,
            molecular_type = mol_type,
            features = features,
            samples = samples,
            return_data = TRUE,
            data_type = data_type,
            tumor_type = tumor_type,
            chunk_size = chunk_size,
            validate_features = validate_features,
            zscore = zscore,
            format = format
          )
        }, error = function(e) {
          warning("Failed to load molecular profile '", mol_type, "': ", e$message)
          if (mol_type %in% c("mutation_gene", "mutation_site", "fusion")) {
            return(data.frame())
          } else {
            return(matrix(nrow = 0, ncol = 0))
          }
        })
      }, mc.cores = min(length(available_types), parallel::detectCores() - 1))

      names(all_data) <- available_types
    } else {
      # Sequential loading for Windows or single core
      all_data <- list()
      for (mol_type in available_types) {
        tryCatch({
          mol_data <- loadMolecularProfiles(
            object = object,
            molecular_type = mol_type,
            features = features,
            samples = samples,
            return_data = TRUE,
            data_type = data_type,
            tumor_type = tumor_type,
            chunk_size = chunk_size,
            validate_features = validate_features,
            zscore = zscore,
            format = format
          )
          all_data[[mol_type]] <- mol_data
          message("Loaded molecular profile: ", mol_type, " (",
                 ifelse(is.matrix(mol_data) || is.data.frame(mol_data),
                       paste(nrow(mol_data), "features x", ncol(mol_data), "samples"),
                       "data loaded"), ")")
        }, error = function(e) {
          warning("Failed to load molecular profile '", mol_type, "': ", e$message)
        })
      }
    }

    # Update object if not returning data
    if (!return_data) {
      for (mol_type in names(all_data)) {
        object@molecularProfiles[[mol_type]] <- all_data[[mol_type]]
      }
      message("Loaded ", length(all_data), " molecular profile types for dataset '", object@name, "'")
      return(object)
    } else {
      return(all_data)
    }
  }

  # Verify we have database connection info
  if (length(object@db_info) == 0 || is.null(object@db_info$db_path)) {
    stop("No database connection information available")
  }

  if (!file.exists(object@db_info$db_path)) {
    stop("Database file not found: ", object@db_info$db_path)
  }

  # Get group prefix
  group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

  # Connect to database with optimized settings
  con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Optimize SQLite settings for performance
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "PRAGMA cache_size = -10000")  # 10MB cache

  # Check table existence
  table_name <- paste0(group_prefix, "_", molecular_type)
  if (!table_name %in% DBI::dbListTables(con)) {
    stop("Molecular profile type '", molecular_type, "' not found for dataset '", object@name, "'")
  }

  # Pre-filter samples based on data_type and tumor_type
  if (data_type != "all" || tumor_type != "all") {
    if (nrow(object@sampleMetadata) == 0) {
      warning("Sample metadata not available for filtering by data_type or tumor_type")
    } else {
      filtered_samples <- object@sampleMetadata$SampleID

      if (data_type != "all") {
        filtered_samples <- object@sampleMetadata$SampleID[object@sampleMetadata$DataType == data_type]
        if (length(filtered_samples) == 0) {
          return(empty_result(molecular_type, return_data, object))
        }
      }

      if (tumor_type != "all") {
        tumor_samples <- object@sampleMetadata$SampleID[object@sampleMetadata$TumorType == tumor_type]
        if (length(tumor_samples) == 0) {
          return(empty_result(molecular_type, return_data, object))
        }
        filtered_samples <- intersect(filtered_samples, tumor_samples)
      }

      samples <- if (is.null(samples)) filtered_samples else intersect(samples, filtered_samples)
      if (length(samples) == 0) {
        return(empty_result(molecular_type, return_data, object))
      }
    }
  }

  # Load data based on molecular type
  if (molecular_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms")) {
    # For matrix data - use optimized query construction
    result <- load_matrix_data(con, table_name, molecular_type, features, samples,
                              chunk_size, validate_features, return_data)

    # Apply z-score normalization if requested
    if (zscore && !is.null(result) && is.matrix(result) && nrow(result) > 0 && ncol(result) > 0) {
      result <- zscoreNormalize(result)
      attr(result, "zscore_normalized") <- TRUE
    }

    if (!return_data) {
      object@molecularProfiles[[molecular_type]] <- result
      return(object)
    } else {
      return(result)
    }

  } else if (molecular_type %in% c("mutation_gene", "mutation_site", "fusion")) {
    # For dataframe data
    result <- load_dataframe_data(con, table_name, molecular_type, features, samples,
                                 validate_features, return_data, format)

    # No z-score normalization for discrete data
    if (zscore) {
      warning("Z-score normalization not applicable for molecular type: ", molecular_type)
    }

    if (!return_data) {
      object@molecularProfiles[[molecular_type]] <- result
      return(object)
    } else {
      return(result)
    }
  } else {
    warning("Unrecognized molecular profile type: ", molecular_type)
    if (return_data) return(NULL)
  }
})

# Helper function for empty results
empty_result <- function(molecular_type, return_data, object) {
  if (return_data) return(matrix(nrow = 0, ncol = 0))
  object@molecularProfiles[[molecular_type]] <- matrix(nrow = 0, ncol = 0)
  return(object)
}

# Helper function to load matrix data with optimizations
load_matrix_data <- function(con, table_name, molecular_type, features, samples,
                            chunk_size, validate_features, return_data) {
  # Build optimized SQL query
  query_parts <- c("SELECT * FROM", table_name)
  where_clauses <- c()

  # Add feature filter
  if (!is.null(features)) {
    if (validate_features) {
      # Single query to check and filter features
      feature_query <- paste0("SELECT DISTINCT feature_id FROM ", table_name,
                             " WHERE feature_id IN (", paste0("'", features, "'", collapse = ","), ")")
      existing_features <- DBI::dbGetQuery(con, feature_query)$feature_id

      missing_features <- setdiff(features, existing_features)
      if (length(missing_features) > 0) {
        warning("The following features do not exist in the ", molecular_type, " data: ",
                paste(missing_features, collapse = ", "))
      }

      if (length(existing_features) == 0) {
        warning("None of the specified features exist in the ", molecular_type, " data")
        return(matrix(nrow = 0, ncol = 0))
      }

      features <- existing_features
    }

    where_clauses <- c(where_clauses, paste0("feature_id IN (", paste0("'", features, "'", collapse = ","), ")"))
  }

  # Add sample filter if specified
  if (!is.null(samples)) {
    # Get all column names except feature_id
    col_info <- DBI::dbGetQuery(con, paste0("PRAGMA table_info(", table_name, ")"))
    sample_cols <- setdiff(col_info$name, "feature_id")
    existing_samples <- intersect(samples, sample_cols)

    if (length(existing_samples) == 0) {
      warning("No samples match the specified filter")
      return(matrix(nrow = 0, ncol = 0))
    }

    if (length(existing_samples) < length(samples)) {
      missing_samples <- setdiff(samples, existing_samples)
      warning("The following samples do not exist: ", paste(missing_samples, collapse = ", "))
    }

    # Sample filtering will be done after data retrieval
    samples <- existing_samples
  }

  # Combine WHERE clauses
  if (length(where_clauses) > 0) {
    query_parts <- c(query_parts, "WHERE", paste(where_clauses, collapse = " AND "))
  }

  query <- paste(query_parts, collapse = " ")

  # Check table size for chunked loading
  size_query <- paste0("SELECT COUNT(*) FROM ", table_name)
  total_rows <- DBI::dbGetQuery(con, size_query)[1, 1]

  if (total_rows > chunk_size && requireNamespace("data.table", quietly = TRUE)) {
    # Use chunked loading with data.table
    return(load_data_chunks(con, query, table_name, features, samples, chunk_size))
  } else {
    # Standard loading
    data <- DBI::dbGetQuery(con, query)

    if (nrow(data) == 0) {
      warning("No data found for molecular profile type: ", molecular_type)
      return(matrix(nrow = 0, ncol = 0))
    }

    # Efficient conversion to matrix
    feature_ids <- data$feature_id
    data$feature_id <- NULL

    if (requireNamespace("data.table", quietly = TRUE)) {
      mat <- as.matrix(data.table::as.data.table(data))
    } else {
      mat <- as.matrix(data)
    }
    rownames(mat) <- feature_ids

    # Filter samples if specified
    if (!is.null(samples)) {
      mat <- mat[, samples, drop = FALSE]
    }

    return(mat)
  }
}

# Helper function to load dataframe data
load_dataframe_data <- function(con, table_name, molecular_type, features, samples, validate_features, return_data, format = "long") {
  # Check table structure
  col_info <- DBI::dbGetQuery(con, paste0("PRAGMA table_info(", table_name, ")"))
  col_names <- col_info$name

  # Standard format: features in 'features' column, samples in 'samples' column
  # Get all unique samples from database if not specified
  if (is.null(samples) && "samples" %in% col_names) {
    all_samples_query <- paste0("SELECT DISTINCT samples FROM ", table_name)
    samples <- DBI::dbGetQuery(con, all_samples_query)$samples
  }
  
  query_parts <- c("SELECT * FROM", table_name)
  where_clauses <- c()

  # Add feature filter
  if (!is.null(features)) {
    if (validate_features) {
      feature_query <- paste0("SELECT DISTINCT features FROM ", table_name,
                          " WHERE features IN (", paste0("'", features, "'", collapse = ","), ")")
      existing_features <- DBI::dbGetQuery(con, feature_query)$features

      missing_features <- setdiff(features, existing_features)
      if (length(missing_features) > 0) {
        warning("The following features do not exist: ", paste(missing_features, collapse = ", "))
      }

      if (length(existing_features) == 0) {
        warning("None of the specified features exist")
        return(data.frame())
      }

      features <- existing_features
    }

    where_clauses <- c(where_clauses, paste0("features IN (", paste0("'", features, "'", collapse = ","), ")"))
  }

  # Add sample filter
  if (!is.null(samples) && "samples" %in% col_names) {
    # Verify which samples exist in the database
    sample_check_query <- paste0("SELECT DISTINCT samples FROM ", table_name)
    all_db_samples <- DBI::dbGetQuery(con, sample_check_query)$samples
    
    existing_samples <- intersect(samples, all_db_samples)
    
    if (length(existing_samples) == 0) {
      warning("No matching samples found")
      return(data.frame())
    }
    
    if (length(existing_samples) < length(samples)) {
      missing_samples <- setdiff(samples, existing_samples)
      warning("The following samples do not exist: ", paste(missing_samples, collapse = ", "))
    }
    
    where_clauses <- c(where_clauses, paste0("samples IN (", paste0("'", existing_samples, "'", collapse = ","), ")"))
    # Update samples variable for later use in wide format conversion
    samples <- existing_samples
  }

  # Combine WHERE clauses
  if (length(where_clauses) > 0) {
    query_parts <- c(query_parts, "WHERE", paste(where_clauses, collapse = " AND "))
  }

  query <- paste(query_parts, collapse = " ")
  data <- DBI::dbGetQuery(con, query)

  if (nrow(data) == 0) {
    warning("No data found for molecular profile type: ", molecular_type)
    return(data.frame())
  }

  # If format is "long", return the data directly without conversion
  if (format == "long") {
    return(data)
  }

  # Reshape long format to wide format for consistency
  if ("features" %in% colnames(data) && "samples" %in% colnames(data)) {
    # Convert from long to wide format
    message("Converting long discrete data format to wide format...")

    # Get unique features and samples
    unique_features <- unique(data$features)
    # Use all samples (from database query or specified parameter)
    # This ensures samples without mutations are included and marked as 0
    unique_samples <- samples

    if (length(unique_features) == 0 || length(unique_samples) == 0) {
      return(data.frame())
    }

    # Create wide format matrix initialized with 0 (no mutation)
    wide_data <- matrix(0, nrow = length(unique_features), ncol = length(unique_samples))
    rownames(wide_data) <- unique_features
    colnames(wide_data) <- unique_samples

    # Fill in values using vectorized indexing (mark mutations with 1)
    indices <- cbind(
      match(data$features, unique_features),
      match(data$samples, unique_samples)
    )
    valid_idx <- !is.na(indices[,1]) & !is.na(indices[,2])
    if (any(valid_idx)) {
      wide_data[indices[valid_idx, , drop = FALSE]] <- 1
    }

    # Return as matrix with feature names as row names
    return(wide_data)
  }

  # If we reach here with data that couldn't be converted,
  # return empty matrix
  return(matrix(nrow = 0, ncol = 0))
}

# Helper function for chunked loading
load_data_chunks <- function(con, base_query, table_name, features, samples, chunk_size) {
  # Get total count for progress
  count_query <- gsub("SELECT \\*", "SELECT COUNT(*)", base_query)
  total_count <- DBI::dbGetQuery(con, count_query)[1, 1]

  # Get sample columns
  col_info <- DBI::dbGetQuery(con, paste0("PRAGMA table_info(", table_name, ")"))
  sample_cols <- setdiff(col_info$name, "feature_id")

  # Initialize empty matrix
  result_matrix <- matrix(nrow = 0, ncol = length(sample_cols))
  colnames(result_matrix) <- sample_cols

  # Process in chunks
  offset <- 0
  processed <- 0

  while (processed < total_count) {
    chunk_query <- paste0(base_query, " LIMIT ", chunk_size, " OFFSET ", offset)
    chunk_data <- DBI::dbGetQuery(con, chunk_query)

    if (nrow(chunk_data) > 0) {
      feature_ids <- chunk_data$feature_id
      chunk_data$feature_id <- NULL

      chunk_matrix <- as.matrix(chunk_data)
      rownames(chunk_matrix) <- feature_ids

      # Combine with results
      if (nrow(result_matrix) == 0) {
        result_matrix <- chunk_matrix
      } else {
        # Use rbind efficiently
        result_matrix <- rbind(result_matrix, chunk_matrix)
      }
    }

    processed <- processed + nrow(chunk_data)
    offset <- offset + chunk_size

    # Progress message for large datasets
    if (total_count > chunk_size * 10) {
      message(sprintf("Processed %d/%d rows (%.1f%%)", processed, total_count,
                      (processed/total_count)*100))
    }
  }

  # Final filtering if needed
  if (!is.null(samples)) {
    common_samples <- intersect(colnames(result_matrix), samples)
    if (length(common_samples) > 0) {
      result_matrix <- result_matrix[, common_samples, drop = FALSE]
    }
  }

  return(result_matrix)
}

#' Load Treatment Response from Database
#'
#' @description Loads drug response data from the database into a DromaSet object
#' @param object A DromaSet object
#' @param drugs Optional vector of drug names to load. If NULL, loads all drugs.
#' @param samples Optional vector of sample IDs to load. If NULL, loads all samples.
#' @param return_data Logical, if TRUE returns the loaded data directly instead of updating the object (default: FALSE)
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")
#' @param zscore Logical, whether to apply z-score normalization (default: FALSE)
#' @return Updated DromaSet object with loaded drug response data or the loaded data directly if return_data=TRUE
#' @export
setGeneric("loadTreatmentResponse", function(object, drugs = NULL, samples = NULL, return_data = FALSE, data_type = "all", tumor_type = "all", zscore = FALSE)
  standardGeneric("loadTreatmentResponse"))

#' @rdname loadTreatmentResponse
#' @export
setMethod("loadTreatmentResponse", "DromaSet", function(object, drugs = NULL, samples = NULL, return_data = FALSE, data_type = "all", tumor_type = "all", zscore = FALSE) {
  # Verify we have database connection info
  if (length(object@db_info) == 0 || is.null(object@db_info$db_path)) {
    stop("No database connection information available")
  }

  if (!file.exists(object@db_info$db_path)) {
    stop("Database file not found: ", object@db_info$db_path)
  }

  # Get group prefix (if specified) or use dataset name
  group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if the drug table exists
  table_name <- paste0(group_prefix, "_drug")
  all_tables <- DBI::dbListTables(con)

  if (!table_name %in% all_tables) {
    stop("Treatment response data not found for dataset '", object@name, "'")
  }

  # Filter samples by data_type and tumor_type if specified
  if (data_type != "all" || tumor_type != "all") {
    if (nrow(object@sampleMetadata) == 0) {
      warning("Sample metadata not available for filtering by data_type or tumor_type")
    } else {
      filtered_samples <- object@sampleMetadata$SampleID

      if (data_type != "all") {
        filtered_samples <- object@sampleMetadata$SampleID[object@sampleMetadata$DataType == data_type]
        if (length(filtered_samples) == 0) {
          warning("No samples match the specified data_type: ", data_type)
          if (return_data) {
            return(matrix(nrow = 0, ncol = 0))
          } else {
            object@treatmentResponse[["drug"]] <- matrix(nrow = 0, ncol = 0)
            return(object)
          }
        }
      }

      if (tumor_type != "all") {
        tumor_samples <- object@sampleMetadata$SampleID[object@sampleMetadata$TumorType == tumor_type]
        if (length(tumor_samples) == 0) {
          warning("No samples match the specified tumor_type: ", tumor_type)
          if (return_data) {
            return(matrix(nrow = 0, ncol = 0))
          } else {
            object@treatmentResponse[["drug"]] <- matrix(nrow = 0, ncol = 0)
            return(object)
          }
        }
        filtered_samples <- intersect(filtered_samples, tumor_samples)
      }

      # Update samples parameter with filtered samples
      if (is.null(samples)) {
        samples <- filtered_samples
      } else {
        samples <- intersect(samples, filtered_samples)
        if (length(samples) == 0) {
          warning("No samples match both the specified filters and sample list")
          if (return_data) {
            return(matrix(nrow = 0, ncol = 0))
          } else {
            object@treatmentResponse[["drug"]] <- matrix(nrow = 0, ncol = 0)
            return(object)
          }
        }
      }
    }
  }

  # Construct query
  query <- paste0("SELECT * FROM ", table_name)

  # Add drug filter if specified
  if (!is.null(drugs)) {
    # First check if all drugs exist in the database
    drug_check_query <- paste0("SELECT DISTINCT feature_id FROM ", table_name)
    all_drugs <- DBI::dbGetQuery(con, drug_check_query)$feature_id

    missing_drugs <- setdiff(drugs, all_drugs)
    if (length(missing_drugs) > 0) {
      warning("The following drugs do not exist in the treatment response data: ",
              paste(missing_drugs, collapse = ", "))
    }

    # Only query for drugs that exist
    existing_drugs <- intersect(drugs, all_drugs)
    if (length(existing_drugs) == 0) {
      warning("None of the specified drugs exist in the treatment response data")
      if (return_data) {
        return(matrix(nrow = 0, ncol = 0))
      } else {
        object@treatmentResponse[["drug"]] <- matrix(nrow = 0, ncol = 0)
        return(object)
      }
    }

    drugs_str <- paste0("'", existing_drugs, "'", collapse = ",")
    query <- paste0(query, " WHERE feature_id IN (", drugs_str, ")")
  }

  # Execute query
  data <- DBI::dbGetQuery(con, query)

  # Variable to store the loaded data
  loaded_data <- NULL

  # Reshape to matrix format
  if (nrow(data) > 0) {
    # Extract feature_id column (drug names)
    drug_ids <- data$feature_id
    data$feature_id <- NULL

    # Convert to matrix
    mat <- as.matrix(data)
    rownames(mat) <- drug_ids

    # Filter by samples if needed
    if (!is.null(samples)) {
      # Check for missing samples
      missing_samples <- setdiff(samples, colnames(mat))
      if (length(missing_samples) > 0) {
        warning("The following samples do not exist in the treatment response data: ",
                paste(missing_samples, collapse = ", "))
      }

      common_samples <- intersect(colnames(mat), samples)
      if (length(common_samples) == 0) {
        warning("No samples match the specified filter")
        mat <- matrix(nrow = 0, ncol = 0)
      } else {
        mat <- mat[, common_samples, drop = FALSE]
      }
    }

    # Store the loaded data
    loaded_data <- mat

    # Apply z-score normalization if requested
    if (zscore && !is.null(loaded_data) && is.matrix(loaded_data) && nrow(loaded_data) > 0 && ncol(loaded_data) > 0) {
      loaded_data <- zscoreNormalize(loaded_data)
      attr(loaded_data, "zscore_normalized") <- TRUE
    }

    # Store in object if not returning data directly
    if (!return_data) {
      object@treatmentResponse[["drug"]] <- loaded_data
    }
  } else {
    warning("No data found for treatment response")
    loaded_data <- matrix(nrow = 0, ncol = 0)

    if (!return_data) {
      object@treatmentResponse[["drug"]] <- matrix(nrow = 0, ncol = 0)
    }
  }

  # Return either the updated object or the loaded data
  if (return_data) {
    return(loaded_data)
  } else {
    return(object)
  }
})
