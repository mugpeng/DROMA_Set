#' Connect to DROMA Database
#'
#' @description Establishes a connection to the DROMA SQLite database
#' @param db_path Path to the SQLite database file
#' @return A database connection object
#' @export
connectDROMADatabase <- function(db_path = file.path(path.expand("sql_db"), "droma.sqlite")) {
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Please install them with install.packages(c('RSQLite', 'DBI'))")
  }

  if (!file.exists(db_path)) {
    stop("Database file not found. Create the database first with createDROMADatabase()")
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  message("Connected to DROMA database at ", db_path)

  # Store the connection in the package environment
  assign("droma_db_connection", con, envir = .GlobalEnv)

  # Register an exit handler to close the connection when R exits
  reg.finalizer(.GlobalEnv, function(e) {
    if (exists("droma_db_connection", envir = e) &&
        inherits(get("droma_db_connection", envir = e), "DBIConnection")) {
      DBI::dbDisconnect(get("droma_db_connection", envir = e))
    }
  }, onexit = TRUE)

  return(con)
}

#' Update DROMA Database with New Object
#'
#' @description Adds or updates a table in the DROMA database with a new object
#' @param obj The object to add to the database (matrix or data.frame)
#' @param table_name The name to use for the table in the database
#' @param overwrite Logical, whether to overwrite if table already exists (default: FALSE)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return Invisibly returns TRUE if successful
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Create a matrix of gene expression data
#' expr_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(expr_data) <- paste0("gene", 1:10)
#' colnames(expr_data) <- paste0("sample", 1:10)
#'
#' # Add to database as "myproject_mRNA"
#' updateDROMADatabase(expr_data, "myproject_mRNA", overwrite = TRUE)
#' }
updateDROMADatabase <- function(obj, table_name, overwrite = FALSE, connection = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Packages 'DBI' and 'RSQLite' are required. Please install them with install.packages(c('DBI', 'RSQLite'))")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Check if table already exists
  all_tables <- DBI::dbListTables(connection)
  if (table_name %in% all_tables && !overwrite) {
    stop("Table '", table_name, "' already exists. Set overwrite = TRUE to replace it.")
  } else if (table_name %in% all_tables && overwrite) {
    message("Overwriting existing table '", table_name, "'")
  }

  # Process the object based on its type
  if (is.matrix(obj)) {
    # Convert matrix to data frame with feature_id column
    df <- as.data.frame(obj)
    df$feature_id <- rownames(df)

    # Write to database
    DBI::dbWriteTable(connection, table_name, df, overwrite = TRUE)

    # Create index on feature_id for faster lookups
    DBI::dbExecute(connection, paste0("CREATE INDEX IF NOT EXISTS idx_", table_name, "_feature_id ON ",
                                     table_name, " (feature_id)"))

    message("Added matrix to database as '", table_name, "' with ", nrow(df), " rows and ",
            ncol(df) - 1, " columns")
  } else if (is.data.frame(obj)) {
    obj <- as.data.frame(obj)
    # Check if the data frame has rownames and preserve them
    if (!is.null(rownames(obj)) && !all(rownames(obj) == seq_len(nrow(obj)))) {
      # If data frame has meaningful rownames, add them as feature_id column
      df <- obj
      df$feature_id <- rownames(df)

      # Write to database
      DBI::dbWriteTable(connection, table_name, df, overwrite = TRUE)

      # Create index on feature_id for faster lookups
      DBI::dbExecute(connection, paste0("CREATE INDEX IF NOT EXISTS idx_", table_name, "_feature_id ON ",
                                       table_name, " (feature_id)"))

      message("Added data frame to database as '", table_name, "' with ", nrow(df), " rows and ",
              ncol(df) - 1, " columns")
    } else {
      # If no meaningful rownames, write directly
      DBI::dbWriteTable(connection, table_name, obj, overwrite = TRUE)
      message("Added data frame to database as '", table_name, "' with ", nrow(obj), " rows and ",
              ncol(obj), " columns")
    }
  } else {
    stop("Object must be a matrix or data frame")
  }

  # Update projects table if it exists
  if ("projects" %in% all_tables) {
    # Extract project name from table name
    parts <- strsplit(table_name, "_")[[1]]
    if (length(parts) >= 2) {
      project_name <- parts[1]
      data_type <- paste(parts[-1], collapse = "_")

      # Check if project exists in projects table
      projects_df <- DBI::dbReadTable(connection, "projects")

      if (project_name %in% projects_df$project_name) {
        # Update existing project entry
        project_row <- projects_df[projects_df$project_name == project_name, ]

        # Update data_types field
        current_types <- unlist(strsplit(project_row$data_types, ","))
        if (!data_type %in% current_types) {
          new_types <- sort(c(current_types, data_type))
          project_row$data_types <- paste(new_types, collapse = ",")

          # Update the row in the projects table
          DBI::dbExecute(connection, paste0(
            "UPDATE projects SET data_types = '", project_row$data_types,
            "' WHERE project_name = '", project_name, "'"
          ))

          message("Updated project '", project_name, "' with new data type '", data_type, "'")
        }
      } else {
        # Create new project entry
        new_project <- data.frame(
          project_name = project_name,
          dataset_type = NA_character_,
          data_types = data_type,
          sample_count = length(setdiff(colnames(obj), "feature_id")),
          drug_count = if(data_type == "drug") nrow(obj) else 0,
          stringsAsFactors = FALSE
        )

        # Add to projects table
        DBI::dbWriteTable(connection, "projects", new_project, append = TRUE)
        message("Added new project '", project_name, "' to projects table")
      }
    }
  }

  invisible(TRUE)
}

#' Retrieve Feature Data from DROMA Database
#'
#' @description Fetches specific feature data from the DROMA database based on selection criteria
#' @param select_feas_type The type of feature to select (e.g., "mRNA", "cnv", "drug")
#' @param select_feas The specific feature to select within the feature type
#' @param data_sources Vector of data sources to select from (e.g., c("ccle", "gdsc"))
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO", "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or specific tumor type
#' @param connection Optional database connection object. If NULL, uses global connection.
#' @return A list of selected features from specified data sources
#' @export
#' @note This function is provided for backward compatibility. For new code, consider
#' using the DromaSet object approach instead.
getFeatureFromDatabase <- function(select_feas_type, select_feas,
                                 data_sources = "all",
                                 data_type = "all", tumor_type = "all",
                                 connection = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Get data source tables that match the feature type
  all_tables <- DBI::dbListTables(connection)
  pattern <- paste0("_", select_feas_type, "$")
  feature_tables <- grep(pattern, all_tables, value = TRUE)

  if (length(feature_tables) == 0) {
    stop("No tables found for feature type: ", select_feas_type)
  }

  # Filter by data sources if specified
  if (!identical(data_sources, "all")) {
    feature_tables <- grep(paste0("^(", paste(data_sources, collapse = "|"), ")_"),
                           feature_tables, value = TRUE)
  }

  if (length(feature_tables) == 0) {
    stop("No matching tables found for the specified data sources")
  }

  # Get sample IDs from sample_anno based on data_type and tumor_type
  filtered_samples <- NULL
  if (data_type != "all" || tumor_type != "all") {
    # Construct SQL query for sample filtering
    sample_query <- "SELECT SampleID FROM sample_anno WHERE 1=1"

    if (data_type != "all") {
      sample_query <- paste0(sample_query, " AND DataType = '", data_type, "'")
    }

    if (tumor_type != "all") {
      sample_query <- paste0(sample_query, " AND TumorType = '", tumor_type, "'")
    }

    # Execute the query
    filtered_samples <- DBI::dbGetQuery(connection, sample_query)$SampleID

    if (length(filtered_samples) == 0) {
      stop("No samples match the specified data_type and tumor_type criteria")
    }
  }

  # Retrieve data for each table
  result_list <- list()

  for (table in feature_tables) {
    # Extract data source name from table name
    data_source <- sub(paste0("_", select_feas_type, "$"), "", table)

    # Query for the specified feature
    if (select_feas_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
      # For continuous data, get the row for the feature
      query <- paste0("SELECT * FROM ", table, " WHERE feature_id = '", select_feas, "'")
      feature_data <- DBI::dbGetQuery(connection, query)

      if (nrow(feature_data) == 0) {
        next  # Skip if feature not found
      }

      # Convert to vector format (excluding feature_id column)
      feature_vector <- as.numeric(as.vector(feature_data[1, -which(names(feature_data) == "feature_id")]))
      names(feature_vector) <- colnames(feature_data)[-which(names(feature_data) == "feature_id")]
    } else {
      # For discrete data like mutations, get sample IDs where feature is present
      query <- paste0("SELECT cells FROM ", table, " WHERE gene = '", select_feas, "'")
      feature_data <- DBI::dbGetQuery(connection, query)

      if (nrow(feature_data) == 0) {
        next  # Skip if feature not found
      }

      feature_vector <- feature_data$cells
    }

    # Filter by samples if needed
    if (!is.null(filtered_samples)) {
      if (select_feas_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
        common_samples <- intersect(names(feature_vector), filtered_samples)
        if (length(common_samples) == 0) {
          next  # Skip if no samples match the filter
        }
        feature_vector <- feature_vector[common_samples]
      } else {
        feature_vector <- intersect(feature_vector, filtered_samples)
        if (length(feature_vector) == 0) {
          next  # Skip if no samples match the filter
        }
      }
    }

    # Add to result list
    result_list[[data_source]] <- feature_vector
  }

  if (length(result_list) == 0) {
    stop("No data found for feature '", select_feas, "' with the specified criteria")
  }

  return(result_list)
}

#' Close DROMA Database Connection
#'
#' @description Closes the connection to the DROMA database
#' @param connection Optional database connection object. If NULL, uses global connection.
#' @return TRUE if successfully disconnected
#' @export
closeDROMADatabase <- function(connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      message("No open database connection found")
      return(invisible(FALSE))
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  DBI::dbDisconnect(connection)

  if (exists("droma_db_connection", envir = .GlobalEnv)) {
    rm("droma_db_connection", envir = .GlobalEnv)
  }

  message("Database connection closed")
  invisible(TRUE)
}

#' List Available Tables in DROMA Database
#'
#' @description Provides information about tables available in the DROMA database
#' @param pattern Optional regex pattern to filter table names
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return A data frame with table information
#' @export
listDROMADatabaseTables <- function(pattern = NULL, connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Get table list
  tables <- DBI::dbListTables(connection)

  # Filter by pattern if provided
  if (!is.null(pattern)) {
    tables <- grep(pattern, tables, value = TRUE)
  }

  # Get metadata for each table
  if ("droma_metadata" %in% DBI::dbListTables(connection)) {
    metadata <- DBI::dbReadTable(connection, "droma_metadata")
    result <- metadata[metadata$table_name %in% tables, ]
  } else {
    # If metadata table doesn't exist, create basic info
    result <- data.frame(
      table_name = tables,
      row_count = sapply(tables, function(t)
        DBI::dbGetQuery(connection, paste0("SELECT COUNT(*) FROM ", t))[1,1]),
      stringsAsFactors = FALSE
    )
  }

  # Add categorization by data type
  result$data_type <- sub("_.*$", "", result$table_name)
  result$feature_type <- sub("^.*_", "", result$table_name)

  return(result)
}

#' List Available Projects in DROMA Database
#'
#' @description Lists all projects available in the DROMA database
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param show_names_only Logical, if TRUE returns only a character vector of project names
#' @param project_data_types Character, project name to get specific data types for
#' @return A data frame with project information or a character vector of project names or data types
#' @export
listDROMAProjects <- function(connection = NULL, show_names_only = FALSE, project_data_types = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Check if projects table exists
  if ("projects" %in% DBI::dbListTables(connection)) {
    projects_df <- DBI::dbReadTable(connection, "projects")

    # If user just wants project names, return them
    if (show_names_only) {
      return(projects_df$project_name)
    }

    # If user wants data types for a specific project
    if (!is.null(project_data_types)) {
      if (project_data_types %in% projects_df$project_name) {
        project_row <- projects_df[projects_df$project_name == project_data_types, ]
        # Split the data_types string by comma and return as a vector
        return(unlist(strsplit(project_row$data_types, ",")))
      } else {
        warning("Project '", project_data_types, "' not found")
        return(character(0))
      }
    }

    return(projects_df)
  }

  # Otherwise, try to infer projects from table names
  tables <- DBI::dbListTables(connection)

  # Extract project names from table prefixes
  project_names <- unique(sapply(tables, function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 2 && !t %in% c("sample_anno", "drug_anno", "droma_metadata", "search_vectors")) {
      return(parts[1])
    } else {
      return(NA)
    }
  }))
  project_names <- project_names[!is.na(project_names)]

  if (length(project_names) == 0) {
    message("No projects found in database")
    if (show_names_only) {
      return(character(0))
    }
    return(data.frame())
  }

  # Return just project names if requested
  if (show_names_only) {
    return(project_names)
  }

  # If user wants data types for a specific project
  if (!is.null(project_data_types)) {
    if (project_data_types %in% project_names) {
      # Get all tables for this project
      project_tables <- grep(paste0("^", project_data_types, "_"), tables, value = TRUE)
      # Extract data types from table names
      data_types <- unique(sapply(project_tables, function(t) {
        sub(paste0("^", project_data_types, "_"), "", t)
      }))
      return(data_types)
    } else {
      warning("Project '", project_data_types, "' not found")
      return(character(0))
    }
  }

  # Create a data frame with project information
  result <- data.frame(
    project_name = project_names,
    stringsAsFactors = FALSE
  )

  return(result)
}

#' List Available Features in DROMA Database
#'
#' @description Lists all available features (genes, drugs, etc.) for a specific project and data type
#' @param project_name Character, the name of the project (e.g., "gCSI", "CCLE")
#' @param data_sources Character, the type of data to query (e.g., "mRNA", "cnv", "drug", "mutation_gene")
#' @param data_type Character, filter by data type: "all" (default), "CellLine", "PDO", "PDC", or "PDX"
#' @param tumor_type Character, filter by tumor type: "all" (default) or specific tumor type
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param limit Integer, maximum number of features to return (default: NULL for all features)
#' @param pattern Character, optional regex pattern to filter feature names
#' @return A character vector of available feature names
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # List all genes available in gCSI mRNA data
#' genes <- listDROMAFeatures("gCSI", "mRNA")
#'
#' # List all drugs available in gCSI drug response data
#' drugs <- listDROMAFeatures("gCSI", "drug")
#'
#' # List genes matching a pattern
#' brca_genes <- listDROMAFeatures("gCSI", "mRNA", pattern = "^BRCA")
#'
#' # List genes for cell lines only
#' cell_line_genes <- listDROMAFeatures("gCSI", "mRNA", data_type = "CellLine")
#'
#' # List genes for breast cancer samples only
#' breast_genes <- listDROMAFeatures("gCSI", "mRNA", tumor_type = "breast cancer")
#'
#' # List first 100 features
#' top_genes <- listDROMAFeatures("gCSI", "mRNA", limit = 100)
#' }
listDROMAFeatures <- function(project_name, data_sources, data_type = "all", tumor_type = "all",
                             connection = NULL, limit = NULL, pattern = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Construct table name
  table_name <- paste0(project_name, "_", data_sources)

  # Check if table exists
  all_tables <- DBI::dbListTables(connection)
  if (!table_name %in% all_tables) {
    stop("Table '", table_name, "' not found. Available tables: ",
         paste(grep(paste0("^", project_name, "_"), all_tables, value = TRUE), collapse = ", "))
  }

  # Get filtered sample IDs if data_type or tumor_type filters are specified
  filtered_samples <- NULL
  if (data_type != "all" || tumor_type != "all") {
    # Check if sample_anno table exists
    if (!"sample_anno" %in% all_tables) {
      warning("Sample annotation table 'sample_anno' not found. Ignoring data_type and tumor_type filters.")
    } else {
      # Construct query for sample filtering
      sample_query <- "SELECT DISTINCT SampleID FROM sample_anno WHERE ProjectID = ?"
      sample_params <- list(project_name)

      if (data_type != "all") {
        sample_query <- paste0(sample_query, " AND DataType = ?")
        sample_params <- append(sample_params, data_type)
      }

      if (tumor_type != "all") {
        sample_query <- paste0(sample_query, " AND TumorType = ?")
        sample_params <- append(sample_params, tumor_type)
      }

      # Execute the query
      sample_result <- DBI::dbGetQuery(connection, sample_query, params = sample_params)
      filtered_samples <- sample_result$SampleID

      if (length(filtered_samples) == 0) {
        filter_parts <- c()
        if(data_type != "all") filter_parts <- c(filter_parts, paste0("data_type='", data_type, "'"))
        if(tumor_type != "all") filter_parts <- c(filter_parts, paste0("tumor_type='", tumor_type, "'"))
        filter_desc <- paste0(" with ", paste(filter_parts, collapse = " and "))

        message("No samples found for project '", project_name, "'", filter_desc)
        return(character(0))
      }
    }
  }

  # Determine the column name based on data type
  if (data_sources %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
    # For continuous data, features are in feature_id column
    feature_column <- "feature_id"
  } else if (data_sources %in% c("mutation_gene", "mutation_site", "fusion")) {
    # For discrete data, features are in genes column
    feature_column <- "genes"
  } else {
    # Try to detect the column automatically
    # Get column names from the table
    columns_query <- paste0("PRAGMA table_info(", table_name, ")")
    columns_info <- DBI::dbGetQuery(connection, columns_query)

    if ("feature_id" %in% columns_info$name) {
      feature_column <- "feature_id"
    } else if ("genes" %in% columns_info$name) {
      feature_column <- "genes"
    } else {
      stop("Cannot determine feature column for data type '", data_sources,
           "'. Available columns: ", paste(columns_info$name, collapse = ", "))
    }
  }

  # For continuous data types, we need to check which features have data for the filtered samples
  if (!is.null(filtered_samples) && data_sources %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
    # Get column names from the table to see which samples have data
    columns_query <- paste0("PRAGMA table_info(", table_name, ")")
    columns_info <- DBI::dbGetQuery(connection, columns_query)
    available_columns <- columns_info$name[columns_info$name != "feature_id"]

    # Find intersection of filtered samples and available columns
    common_samples <- intersect(filtered_samples, available_columns)

    if (length(common_samples) == 0) {
      message("No data available for the specified sample filters in ", table_name)
      return(character(0))
    }

    # For continuous data, we'll get all features since sample filtering is conceptual
    # (all features exist, but only some samples would be used in analysis)
    query <- paste0("SELECT DISTINCT ", feature_column, " FROM ", table_name,
                    " WHERE ", feature_column, " IS NOT NULL")
  } else {
    # Construct query to get distinct features
    query <- paste0("SELECT DISTINCT ", feature_column, " FROM ", table_name,
                    " WHERE ", feature_column, " IS NOT NULL")
  }

  # Add pattern filter if specified
  if (!is.null(pattern)) {
    # Convert regex pattern to SQL LIKE pattern for basic matching
    # For simple patterns like "^BRCA", convert to "BRCA%"
    # For patterns with *, convert to %
    like_pattern <- pattern

    # Handle common regex patterns
    if (grepl("^\\^", pattern)) {
      # Pattern starts with ^, remove ^ and add % at end
      like_pattern <- paste0(sub("^\\^", "", pattern), "%")
    } else if (grepl("\\$$", pattern)) {
      # Pattern ends with $, remove $ and add % at start
      like_pattern <- paste0("%", sub("\\$$", "", pattern))
    } else if (grepl("^\\^.*\\$$", pattern)) {
      # Pattern has both ^ and $, remove both (exact match)
      like_pattern <- gsub("^\\^|\\$$", "", pattern)
    } else {
      # Default: add % on both sides for contains matching
      like_pattern <- paste0("%", pattern, "%")
    }

    # Replace common regex characters with SQL LIKE equivalents
    like_pattern <- gsub("\\*", "%", like_pattern)
    like_pattern <- gsub("\\.", "_", like_pattern)

    query <- paste0(query, " AND ", feature_column, " LIKE '", like_pattern, "'")
  }

  # Add ordering
  query <- paste0(query, " ORDER BY ", feature_column)

  # Add limit if specified
  if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    query <- paste0(query, " LIMIT ", as.integer(limit))
  }

  # Execute query
  tryCatch({
    result <- DBI::dbGetQuery(connection, query)
    features <- result[[feature_column]]

    if (length(features) == 0) {
      filter_parts <- c()
      if (!is.null(pattern)) filter_parts <- c(filter_parts, paste0("pattern='", pattern, "'"))
      if (data_type != "all") filter_parts <- c(filter_parts, paste0("data_type='", data_type, "'"))
      if (tumor_type != "all") filter_parts <- c(filter_parts, paste0("tumor_type='", tumor_type, "'"))

      filter_desc <- if(length(filter_parts) > 0) paste0(" with ", paste(filter_parts, collapse = " and ")) else ""

      message("No features found in ", table_name, filter_desc)
      return(character(0))
    }

    # Print summary information
    total_query <- paste0("SELECT COUNT(DISTINCT ", feature_column, ") as total FROM ", table_name,
                         " WHERE ", feature_column, " IS NOT NULL")
    total_result <- DBI::dbGetQuery(connection, total_query)
    total_features <- total_result$total

    filter_desc <- ""
    filters <- c()
    if (!is.null(pattern)) filters <- c(filters, paste0("pattern='", pattern, "'"))
    if (data_type != "all") filters <- c(filters, paste0("data_type='", data_type, "'"))
    if (tumor_type != "all") filters <- c(filters, paste0("tumor_type='", tumor_type, "'"))

    if (length(filters) > 0) {
      filter_desc <- paste0(" (filtered by ", paste(filters, collapse = " and "), ")")
    }

    if (!is.null(limit)) {
      message("Showing first ", length(features), " features out of ", total_features,
              " total features in ", table_name, filter_desc)
    } else if (!is.null(pattern)) {
      message("Found ", length(features), " features matching pattern '", pattern,
              "' out of ", total_features, " total features in ", table_name, filter_desc)
    } else {
      message("Found ", length(features), " features in ", table_name, filter_desc)
    }

    return(features)

  }, error = function(e) {
    stop("Error querying features from ", table_name, ": ", e$message)
  })
}

#' List Available Samples in DROMA Database
#'
#' @description Lists all available samples for a specific project, optionally filtered by data type or tumor type
#' @param project_name Character, the name of the project (e.g., "gCSI", "CCLE")
#' @param data_type Character, filter by data type: "all" (default), "CellLine", "PDO", "PDC", or "PDX"
#' @param tumor_type Character, filter by tumor type: "all" (default) or specific tumor type
#' @param data_sources Character, filter by data sources: "all" (default) or specific data type (e.g., "mRNA", "cnv", "drug")
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param limit Integer, maximum number of samples to return (default: NULL for all samples)
#' @param pattern Character, optional regex pattern to filter sample names
#' @return A character vector of available sample IDs
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # List all samples for gCSI project
#' samples <- listDROMASamples("gCSI")
#'
#' # List only cell line samples
#' cell_lines <- listDROMASamples("gCSI", data_type = "CellLine")
#'
#' # List only breast cancer samples
#' breast_samples <- listDROMASamples("gCSI", tumor_type = "breast cancer")
#'
#' # List samples with mRNA data
#' mrna_samples <- listDROMASamples("gCSI", data_sources = "mRNA")
#'
#' # List first 50 samples
#' top_samples <- listDROMASamples("gCSI", limit = 50)
#'
#' # List samples matching a pattern
#' mcf_samples <- listDROMASamples("gCSI", pattern = "^MCF")
#' }
listDROMASamples <- function(project_name, data_sources = "all", data_type = "all", tumor_type = "all",
                            connection = NULL, limit = NULL, pattern = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Check if sample_anno table exists
  all_tables <- DBI::dbListTables(connection)
  if (!"sample_anno" %in% all_tables) {
    stop("Sample annotation table 'sample_anno' not found in database")
  }

  # Get samples with data_sources filter if specified
  filtered_samples_by_data <- NULL
  if (data_sources != "all") {
    # Construct table name for the data source
    data_table_name <- paste0(project_name, "_", data_sources)

    # Check if the data source table exists
    if (!data_table_name %in% all_tables) {
      stop("Data source table '", data_table_name, "' not found. Available tables: ",
           paste(grep(paste0("^", project_name, "_"), all_tables, value = TRUE), collapse = ", "))
    }

    # Get samples that have data in this data source
    if (data_sources %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
      # For continuous data, get column names (excluding feature_id)
      columns_query <- paste0("PRAGMA table_info(", data_table_name, ")")
      columns_info <- DBI::dbGetQuery(connection, columns_query)
      filtered_samples_by_data <- columns_info$name[columns_info$name != "feature_id"]
    } else if (data_sources %in% c("mutation_gene", "mutation_site", "fusion")) {
      # For discrete data, get unique values from cells column
      discrete_query <- paste0("SELECT DISTINCT cells FROM ", data_table_name, " WHERE cells IS NOT NULL")
      discrete_result <- DBI::dbGetQuery(connection, discrete_query)
      filtered_samples_by_data <- discrete_result$cells
    } else {
      # Try to detect automatically
      columns_query <- paste0("PRAGMA table_info(", data_table_name, ")")
      columns_info <- DBI::dbGetQuery(connection, columns_query)

      if ("cells" %in% columns_info$name) {
        # Discrete data
        discrete_query <- paste0("SELECT DISTINCT cells FROM ", data_table_name, " WHERE cells IS NOT NULL")
        discrete_result <- DBI::dbGetQuery(connection, discrete_query)
        filtered_samples_by_data <- discrete_result$cells
      } else {
        # Continuous data
        filtered_samples_by_data <- columns_info$name[columns_info$name != "feature_id"]
      }
    }

    if (length(filtered_samples_by_data) == 0) {
      message("No samples found with data in '", data_sources, "' for project '", project_name, "'")
      return(character(0))
    }
  }

  # Construct query
  query <- "SELECT DISTINCT SampleID FROM sample_anno WHERE ProjectID = ?"
  params <- list(project_name)

  # Add data type filter
  if (data_type != "all") {
    query <- paste0(query, " AND DataType = ?")
    params <- append(params, data_type)
  }

  # Add tumor type filter
  if (tumor_type != "all") {
    query <- paste0(query, " AND TumorType = ?")
    params <- append(params, tumor_type)
  }

  # Add data sources filter by restricting to samples with data
  if (!is.null(filtered_samples_by_data)) {
    if (length(filtered_samples_by_data) > 0) {
      placeholders <- paste(rep("?", length(filtered_samples_by_data)), collapse = ",")
      query <- paste0(query, " AND SampleID IN (", placeholders, ")")
      params <- append(params, as.list(filtered_samples_by_data))
    } else {
      # No samples with data in this data source
      return(character(0))
    }
  }

  # Add pattern filter if specified
  if (!is.null(pattern)) {
    # Convert regex pattern to SQL LIKE pattern for basic matching
    # For simple patterns like "^MCF", convert to "MCF%"
    # For patterns with *, convert to %
    like_pattern <- pattern

    # Handle common regex patterns
    if (grepl("^\\^", pattern)) {
      # Pattern starts with ^, remove ^ and add % at end
      like_pattern <- paste0(sub("^\\^", "", pattern), "%")
    } else if (grepl("\\$$", pattern)) {
      # Pattern ends with $, remove $ and add % at start
      like_pattern <- paste0("%", sub("\\$$", "", pattern))
    } else if (grepl("^\\^.*\\$$", pattern)) {
      # Pattern has both ^ and $, remove both (exact match)
      like_pattern <- gsub("^\\^|\\$$", "", pattern)
    } else {
      # Default: add % on both sides for contains matching
      like_pattern <- paste0("%", pattern, "%")
    }

    # Replace common regex characters with SQL LIKE equivalents
    like_pattern <- gsub("\\*", "%", like_pattern)
    like_pattern <- gsub("\\.", "_", like_pattern)

    query <- paste0(query, " AND SampleID LIKE ?")
    params <- append(params, like_pattern)
  }

  # Add ordering
  query <- paste0(query, " ORDER BY SampleID")

  # Add limit if specified
  if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    query <- paste0(query, " LIMIT ", as.integer(limit))
  }

  # Execute query
  tryCatch({
    result <- DBI::dbGetQuery(connection, query, params = params)
    samples <- result$SampleID

    if (length(samples) == 0) {
      filter_parts <- c()
      if(data_type != "all") filter_parts <- c(filter_parts, paste0("data_type='", data_type, "'"))
      if(tumor_type != "all") filter_parts <- c(filter_parts, paste0("tumor_type='", tumor_type, "'"))
      if(data_sources != "all") filter_parts <- c(filter_parts, paste0("data_sources='", data_sources, "'"))
      if(!is.null(pattern)) filter_parts <- c(filter_parts, paste0("pattern='", pattern, "'"))

      filter_desc <- if(length(filter_parts) > 0) paste0(" with ", paste(filter_parts, collapse = " and ")) else ""

      message("No samples found for project '", project_name, "'", filter_desc)
      return(character(0))
    }

    # Print summary information
    total_query <- "SELECT COUNT(DISTINCT SampleID) as total FROM sample_anno WHERE ProjectID = ?"
    total_result <- DBI::dbGetQuery(connection, total_query, params = list(project_name))
    total_samples <- total_result$total

    filter_desc <- ""
    filters <- c()
    if (data_type != "all") filters <- c(filters, paste0("data_type='", data_type, "'"))
    if (tumor_type != "all") filters <- c(filters, paste0("tumor_type='", tumor_type, "'"))
    if (data_sources != "all") filters <- c(filters, paste0("data_sources='", data_sources, "'"))
    if (!is.null(pattern)) filters <- c(filters, paste0("pattern='", pattern, "'"))

    if (length(filters) > 0) {
      filter_desc <- paste0(" (filtered by ", paste(filters, collapse = " and "), ")")
    }

    if (!is.null(limit)) {
      message("Showing first ", length(samples), " samples out of ", total_samples,
              " total samples for project '", project_name, "'", filter_desc)
    } else {
      if (!is.null(pattern) || data_sources != "all") {
        message("Found ", length(samples), " samples out of ", total_samples,
                " total samples for project '", project_name, "'", filter_desc)
      } else {
        message("Found ", length(samples), " samples for project '", project_name, "'", filter_desc)
      }
    }

    return(samples)

  }, error = function(e) {
    stop("Error querying samples for project '", project_name, "': ", e$message)
  })
}

#' Get Annotation Data from DROMA Database
#'
#' @description Retrieves annotation data from either sample_anno or drug_anno tables
#' @param anno_type Character, type of annotation to retrieve: "sample" or "drug"
#' @param project_name Character, optional project name to filter results (default: NULL for all projects)
#' @param ids Character vector, optional specific IDs to retrieve (SampleID for samples, DrugName for drugs)
#' @param data_type Character, for sample annotations only: filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Character, for sample annotations only: filter by tumor type ("all" or specific type)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param limit Integer, maximum number of records to return (default: NULL for all records)
#' @return A data frame containing the annotation data
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Get all sample annotations
#' sample_anno <- getDROMAAnnotation("sample")
#'
#' # Get sample annotations for gCSI project only
#' gCSI_samples <- getDROMAAnnotation("sample", project_name = "gCSI")
#'
#' # Get sample annotations for specific samples
#' specific_samples <- getDROMAAnnotation("sample", ids = c("22RV1", "2313287"))
#'
#' # Get cell line samples only
#' cell_lines <- getDROMAAnnotation("sample", data_type = "CellLine")
#'
#' # Get breast cancer samples only
#' breast_samples <- getDROMAAnnotation("sample", tumor_type = "breast cancer")
#'
#' # Get all drug annotations
#' drug_anno <- getDROMAAnnotation("drug")
#'
#' # Get drug annotations for gCSI project only
#' gCSI_drugs <- getDROMAAnnotation("drug", project_name = "gCSI")
#'
#' # Get annotations for specific drugs
#' specific_drugs <- getDROMAAnnotation("drug", ids = c("Tamoxifen", "Cisplatin"))
#' }
getDROMAAnnotation <- function(anno_type, project_name = NULL, ids = NULL,
                              data_type = "all", tumor_type = "all",
                              connection = NULL, limit = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Validate anno_type
  if (!anno_type %in% c("sample", "drug")) {
    stop("anno_type must be either 'sample' or 'drug'")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Determine table name and ID column
  if (anno_type == "sample") {
    table_name <- "sample_anno"
    id_column <- "SampleID"
    project_column <- "ProjectID"
  } else {
    table_name <- "drug_anno"
    id_column <- "DrugName"
    project_column <- "ProjectID"
  }

  # Check if table exists
  all_tables <- DBI::dbListTables(connection)
  if (!table_name %in% all_tables) {
    stop("Annotation table '", table_name, "' not found in database")
  }

  # Build query
  query <- paste0("SELECT * FROM ", table_name, " WHERE 1=1")
  params <- list()

  # Add project filter
  if (!is.null(project_name)) {
    query <- paste0(query, " AND ", project_column, " = ?")
    params <- append(params, project_name)
  }

  # Add ID filter
  if (!is.null(ids) && length(ids) > 0) {
    placeholders <- paste(rep("?", length(ids)), collapse = ",")
    query <- paste0(query, " AND ", id_column, " IN (", placeholders, ")")
    params <- append(params, as.list(ids))
  }

  # Add sample-specific filters
  if (anno_type == "sample") {
    # Add data type filter
    if (data_type != "all") {
      query <- paste0(query, " AND DataType = ?")
      params <- append(params, data_type)
    }

    # Add tumor type filter
    if (tumor_type != "all") {
      query <- paste0(query, " AND TumorType = ?")
      params <- append(params, tumor_type)
    }
  }

  # Add ordering
  query <- paste0(query, " ORDER BY ", id_column)

  # Add limit if specified
  if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    query <- paste0(query, " LIMIT ", as.integer(limit))
  }

  # Execute query
  tryCatch({
    if (length(params) > 0) {
      result <- DBI::dbGetQuery(connection, query, params = params)
    } else {
      result <- DBI::dbGetQuery(connection, query)
    }

    if (nrow(result) == 0) {
      filter_desc <- ""
      filters <- c()

      if (!is.null(project_name)) {
        filters <- c(filters, paste0("project='", project_name, "'"))
      }
      if (!is.null(ids)) {
        filters <- c(filters, paste0("specific IDs (", length(ids), " requested)"))
      }
      if (anno_type == "sample") {
        if (data_type != "all") {
          filters <- c(filters, paste0("data_type='", data_type, "'"))
        }
        if (tumor_type != "all") {
          filters <- c(filters, paste0("tumor_type='", tumor_type, "'"))
        }
      }

      if (length(filters) > 0) {
        filter_desc <- paste0(" with filters: ", paste(filters, collapse = ", "))
      }

      message("No ", anno_type, " annotations found", filter_desc)
      return(data.frame())
    }

    # Print summary information
    total_query <- paste0("SELECT COUNT(*) as total FROM ", table_name)
    total_result <- DBI::dbGetQuery(connection, total_query)
    total_records <- total_result$total

    filter_desc <- ""
    if (!is.null(project_name) || !is.null(ids) ||
        (anno_type == "sample" && (data_type != "all" || tumor_type != "all"))) {
      filters <- c()

      if (!is.null(project_name)) {
        filters <- c(filters, paste0("project='", project_name, "'"))
      }
      if (!is.null(ids)) {
        filters <- c(filters, paste0("specific IDs"))
      }
      if (anno_type == "sample") {
        if (data_type != "all") {
          filters <- c(filters, paste0("data_type='", data_type, "'"))
        }
        if (tumor_type != "all") {
          filters <- c(filters, paste0("tumor_type='", tumor_type, "'"))
        }
      }

      filter_desc <- paste0(" (filtered by ", paste(filters, collapse = " and "), ")")
    }

    if (!is.null(limit)) {
      message("Retrieved first ", nrow(result), " ", anno_type, " annotations out of ",
              total_records, " total records", filter_desc)
    } else {
      message("Retrieved ", nrow(result), " ", anno_type, " annotations", filter_desc)
    }

    return(result)

  }, error = function(e) {
    stop("Error querying ", anno_type, " annotations: ", e$message)
  })
}

