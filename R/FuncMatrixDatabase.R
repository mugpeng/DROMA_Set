# Matrix Database Storage Functions ----

#' Store Matrix in SQLite Database
#'
#' @description Converts and stores a matrix or data.frame in a SQLite database with efficient indexing.
#' Row names are preserved as feature_id column.
#' @param db_path Path where the SQLite database file should be created or updated
#' @param matrix Matrix or data.frame to store
#' @param table_name Name of the table to create in database
#' @return Invisibly returns the path to the database
#' @export
#' @examples
#' \dontrun{
#' # Store matrix
#' exp_matrix <- matrix(rnorm(1000), nrow=100, ncol=10)
#' rownames(exp_matrix) <- paste0("Gene_", 1:100)
#' colnames(exp_matrix) <- paste0("Sample_", 1:10)
#' 
#' storeMatricesInDatabase(
#'   db_path = "expression_data.sqlite",
#'   matrix = exp_matrix,
#'   table_name = "expression"
#' )
#' 
#' # Store data.frame
#' df <- data.frame(matrix(runif(1000), nrow=100, ncol=10))
#' rownames(df) <- paste0("Gene_", 1:100)
#' storeMatricesInDatabase("data.sqlite", df, "methylation")
#' }
storeMatricesInDatabase <- function(db_path,
                                   matrix,
                                   table_name) {
  
  # First principles validation: check fundamental requirements
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Install with: install.packages(c('RSQLite', 'DBI'))")
  }
  
  # Validate inputs
  if (missing(db_path) || !is.character(db_path) || length(db_path) != 1) {
    stop("db_path must be a single character string specifying the database file path")
  }
  
  if (missing(matrix)) {
    stop("matrix is required and must be a matrix or data.frame")
  }
  
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    stop("matrix must be a matrix or data.frame")
  }
  
  if (missing(table_name) || !is.character(table_name) || length(table_name) != 1) {
    stop("table_name must be a single character string")
  }
  
  # Validate table name (SQLite requirements)
  if (grepl("^[0-9]|[^a-zA-Z0-9_]", table_name)) {
    stop("Invalid table name '", table_name, 
         "'. Table names must start with letter/underscore and contain only letters, numbers, underscores")
  }
  
  # Create database directory if needed
  db_dir <- dirname(db_path)
  if (!dir.exists(db_dir)) {
    dir.create(db_dir, recursive = TRUE)
    message("Created directory: ", db_dir)
  }
  
  # Connect to database with proper cleanup
  message("Connecting to database: ", db_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Process matrix
  message("Processing matrix: ", table_name)
  
  # Validate matrix has row and column names
  if (is.null(rownames(matrix))) {
    rownames(matrix) <- paste0("feature_", seq_len(nrow(matrix)))
    warning("Matrix has no row names. Generated generic names.")
  }
  
  if (is.null(colnames(matrix))) {
    colnames(matrix) <- paste0("sample_", seq_len(ncol(matrix)))
    warning("Matrix has no column names. Generated generic names.")
  }
  
  # Convert matrix to data frame with feature_id
  tryCatch({
    # Convert to data frame preserving all data types
    df <- as.data.frame(matrix, stringsAsFactors = FALSE)
    
    # Add feature_id column with row names
    df$feature_id <- rownames(matrix)
    
    # Reorder columns to put feature_id first
    df <- df[, c("feature_id", setdiff(names(df), "feature_id"))]
    
    # Store dimensions for logging
    n_features <- nrow(df)
    n_samples <- ncol(df) - 1  # Subtract 1 for feature_id column
    
    # Write to database (overwrite if exists)
    DBI::dbWriteTable(con, table_name, df, overwrite = TRUE)
    
    # Create index on feature_id for efficient queries
    index_name <- paste0("idx_", table_name, "_feature_id")
    index_sql <- paste0("CREATE INDEX IF NOT EXISTS ", index_name, 
                       " ON ", table_name, " (feature_id)")
    DBI::dbExecute(con, index_sql)
    
    message("Stored ", n_features, " features × ", n_samples, 
            " samples with feature_id index")
    
  }, error = function(e) {
    stop("Failed to process matrix '", table_name, "': ", e$message)
  })
  
  # Final summary
  message("Database storage complete. Table '", table_name, "' created.")
  
  invisible(db_path)
}

#' Retrieve Matrix from SQLite Database
#'
#' @description Loads a matrix from SQLite database, reconstructing the original matrix format
#' with row names from feature_id column.
#' @param db_path Path to the SQLite database file
#' @param table_name Name of the table containing the matrix data
#' @param features Optional vector of specific feature IDs to retrieve. If NULL, retrieves all
#' @return Matrix object with feature_id values as row names
#' @export
#' @examples
#' \dontrun{
#' # Retrieve entire matrix
#' exp_matrix <- retrieveMatrixFromDatabase("expression_data.sqlite", "expression")
#' 
#' # Retrieve specific features
#' subset_matrix <- retrieveMatrixFromDatabase(
#'   "expression_data.sqlite", 
#'   "expression", 
#'   features = c("Gene_1", "Gene_5", "Gene_10")
#' )
#' }
retrieveMatrixFromDatabase <- function(db_path,
                                     table_name,
                                     features = NULL) {
  
  # Validate inputs
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Install with: install.packages(c('RSQLite', 'DBI'))")
  }
  
  if (!file.exists(db_path)) {
    stop("Database file does not exist: ", db_path)
  }
  
  if (missing(table_name) || !is.character(table_name) || length(table_name) != 1) {
    stop("table_name must be a single character string")
  }
  
  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Check if table exists
  if (!table_name %in% DBI::dbListTables(con)) {
    stop("Table '", table_name, "' not found in database. Available tables: ",
         paste(DBI::dbListTables(con), collapse = ", "))
  }
  
  # Build query
  if (is.null(features)) {
    query <- paste0("SELECT * FROM ", table_name)
  } else {
    if (!is.character(features)) {
      stop("features must be a character vector of feature IDs")
    }
    # Use parameterized query to prevent SQL injection
    feature_placeholders <- paste0("'", features, "'", collapse = ", ")
    query <- paste0("SELECT * FROM ", table_name, 
                   " WHERE feature_id IN (", feature_placeholders, ")")
  }
  
  # Execute query
  tryCatch({
    df <- DBI::dbGetQuery(con, query)
    
    if (nrow(df) == 0) {
      warning("No data retrieved for table '", table_name, "'")
      return(NULL)
    }
    
    # Check if feature_id column exists
    if (!"feature_id" %in% colnames(df)) {
      stop("Table '", table_name, "' does not have a feature_id column")
    }
    
    # Convert back to matrix
    feature_ids <- df$feature_id
    df$feature_id <- NULL
    
    # Convert to matrix
    mat <- as.matrix(df)
    rownames(mat) <- feature_ids
    
    message("Retrieved matrix: ", nrow(mat), " features × ", ncol(mat), " samples")
    
    return(mat)
    
  }, error = function(e) {
    stop("Failed to retrieve matrix from table '", table_name, "': ", e$message)
  })
}

#' List Matrix Tables in Database
#'
#' @description Lists all matrix tables stored in the database along with metadata
#' @param db_path Path to the SQLite database file
#' @return Data frame with table information including dimensions and creation dates
#' @export
#' @examples
#' \dontrun{
#' # List all matrices in database
#' table_info <- listMatrixTables("expression_data.sqlite")
#' print(table_info)
#' }
listMatrixTables <- function(db_path) {
  
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Install with: install.packages(c('RSQLite', 'DBI'))")
  }
  
  if (!file.exists(db_path)) {
    stop("Database file does not exist: ", db_path)
  }
  
  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  all_tables <- DBI::dbListTables(con)
  
  if (length(all_tables) == 0) {
    message("No tables found in database")
    return(data.frame())
  }
  
  # Check if metadata table exists
  if ("matrix_metadata" %in% all_tables) {
    metadata <- DBI::dbGetQuery(con, "SELECT * FROM matrix_metadata")
    return(metadata)
  } else {
    # Generate metadata on the fly
    matrix_tables <- setdiff(all_tables, c("matrix_metadata", "sqlite_sequence"))
    
    if (length(matrix_tables) == 0) {
      message("No matrix tables found in database")
      return(data.frame())
    }
    
    metadata <- data.frame(
      table_name = matrix_tables,
      stringsAsFactors = FALSE
    )
    
    # Try to get dimensions for each table
    for (i in seq_len(nrow(metadata))) {
      table_name <- metadata$table_name[i]
      tryCatch({
        count_query <- paste0("SELECT COUNT(*) FROM ", table_name)
        n_features <- DBI::dbGetQuery(con, count_query)[1,1]
        
        # Get column count (subtract 1 for feature_id)
        cols_query <- paste0("PRAGMA table_info(", table_name, ")")
        col_info <- DBI::dbGetQuery(con, cols_query)
        n_samples <- nrow(col_info) - 1
        
        metadata$n_features[i] <- n_features
        metadata$n_samples[i] <- n_samples
      }, error = function(e) {
        metadata$n_features[i] <- NA
        metadata$n_samples[i] <- NA
      })
    }
    
    return(metadata)
  }
}
