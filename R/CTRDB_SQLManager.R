# CTRDB SQL Manager Functions ----
# Connection helpers and patient-level accessors for the Clinical Trial Response Database.
# Path-based matrix storage: storeMatricesInDatabase(), listMatrixTables()
# in DROMA_SQLManager.R; multi-feature reads: getFeatureFromDatabase().

#' Get patient expression data from CTRDB database
#'
#' @description Retrieves expression data for a specific patient from CTRDB database.
#' Automatically detects and applies log transformation if needed based on data range.
#' @param patient_id Patient identifier (will be sanitized for database table name)
#' @param connection Database connection object
#' @param auto_log Logical, whether to automatically detect and apply log transformation (default: TRUE)
#' @param log_threshold Maximum value threshold to determine if data is already log-transformed (default: 50)
#' @return Expression matrix with genes as rows and samples as columns, or NULL if retrieval fails
#' @export
#' @examples
#' \dontrun{
#' # Connect to CTRDB database
#' con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#'
#' # Get expression data for a patient
#' expr_data <- getPatientExpressionData("Patient-001", con)
#' }
getPatientExpressionData <- function(patient_id,
                                    connection,
                                    auto_log = TRUE,
                                    log_threshold = 50) {

  # Validate inputs
  if (missing(patient_id) || is.null(patient_id) || patient_id == "") {
    stop("patient_id must be specified")
  }
  if (is.null(connection)) {
    stop("Database connection must be provided")
  }

  # Replace "-" with "_" for database table name
  table_name_db <- gsub("-", "_", patient_id)

  # Query expression matrix
  expr_data <- tryCatch({
    expr_query <- paste0("SELECT * FROM ", table_name_db)
    full_expr <- DBI::dbGetQuery(connection, expr_query)

    # Check if data was retrieved
    if (is.null(full_expr) || nrow(full_expr) == 0) {
      warning("No data found in table ", table_name_db, " for patient ", patient_id)
      return(NULL)
    }

    # Convert to matrix with feature_id as rownames
    if (!"feature_id" %in% colnames(full_expr)) {
      stop("feature_id column not found in expression data for patient ", patient_id)
    }

    rownames(full_expr) <- full_expr$feature_id
    full_expr$feature_id <- NULL

    expr_matrix <- as.matrix(full_expr)

    # Check if log transformation is needed
    if (auto_log && nrow(expr_matrix) > 0 && ncol(expr_matrix) > 0) {
      # Calculate statistics to determine if data is log-transformed
      max_value <- max(expr_matrix, na.rm = TRUE)
      min_value <- min(expr_matrix, na.rm = TRUE)
      median_value <- median(expr_matrix, na.rm = TRUE)

      # If min value is negative, skip log transformation
      if (min_value < 0) {
        message("Expression data for patient ", patient_id,
                " contains negative values (min = ", round(min_value, 2),
                "). Skipping log transformation.")
      } else if (max_value > log_threshold || median_value > 20) {
        # If max value exceeds threshold, likely not log-transformed
        # Also check if median is suspiciously high (> 20 typically indicates non-log data)
        message("Expression data for patient ", patient_id,
                " appears to be non-log-transformed (max = ", round(max_value, 2),
                ", median = ", round(median_value, 2),
                "). Applying log2(x + 1) transformation...")
        expr_matrix <- log2(expr_matrix + 1)
      } else {
        message("Expression data for patient ", patient_id,
                " appears to be log-transformed (max = ", round(max_value, 2),
                ", median = ", round(median_value, 2), ")")
      }
    }

    return(expr_matrix)

  }, error = function(e) {
    warning("Failed to retrieve expression data for patient ", patient_id, ": ", e$message)
    return(NULL)
  })

  return(expr_data)
}

#' Connect to CTRDB Database
#'
#' @description Establishes a connection to the Clinical Trial Response Database (CTRDB) SQLite database
#' @param db_path Path to the CTRDB SQLite database file
#' @return A database connection object
#' @export
#' @examples
#' \dontrun{
#' # Connect to CTRDB database
#' con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#' }
connectCTRDatabase <- function(db_path = file.path(path.expand("sql_db"), "ctrdb.sqlite")) {
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Please install them with install.packages(c('RSQLite', 'DBI'))")
  }

  if (!file.exists(db_path)) {
    stop("CTRDB database file not found at: ", db_path)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  message("Connected to CTRDB database at ", db_path)

  # Store the connection in the package environment with a different name
  assign("ctrdb_connection", con, envir = .GlobalEnv)

  # Register an exit handler to close the connection when R exits
  reg.finalizer(.GlobalEnv, function(e) {
    if (exists("ctrdb_connection", envir = e) &&
        inherits(get("ctrdb_connection", envir = e), "DBIConnection")) {
      DBI::dbDisconnect(get("ctrdb_connection", envir = e))
    }
  }, onexit = TRUE)

  return(con)
}

#' Close CTRDB Database Connection
#'
#' @description Closes the connection to the CTRDB database
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return TRUE if successfully disconnected
#' @export
#' @examples
#' \dontrun{
#' # Close CTRDB database connection
#' closeCTRDatabase()
#' }
closeCTRDatabase <- function(connection = NULL) {
  if (is.null(connection)) {
    if (!exists("ctrdb_connection", envir = .GlobalEnv)) {
      message("No open CTRDB database connection found")
      return(invisible(FALSE))
    }
    connection <- get("ctrdb_connection", envir = .GlobalEnv)
  }

  DBI::dbDisconnect(connection)

  if (exists("ctrdb_connection", envir = .GlobalEnv)) {
    rm("ctrdb_connection", envir = .GlobalEnv)
  }

  message("CTRDB database connection closed")
  invisible(TRUE)
}
