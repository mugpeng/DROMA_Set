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

  # Note: Project tracking is handled separately by updateDROMAProjects()

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
#' @param connection Optional database connection object. If NULL, uses global connection
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
#' @description Provides information about omics and drug tables in the DROMA database.
#' Excludes annotation tables, projects table, and backup tables.
#' @param pattern Optional regex pattern to filter table names
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return A data frame with table information including created_date from projects table
#' @export
listDROMADatabaseTables <- function(pattern = NULL, connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Get table list
  all_tables <- DBI::dbListTables(connection)

  # Filter to only omics and drug tables, excluding system tables and backups
  tables <- all_tables[!all_tables %in% c("sample_anno", "drug_anno", "projects", "droma_metadata", "search_vectors")]
  # Remove backup tables like "_raw"
  tables <- tables[!grepl("_raw$", tables)]
  # Only keep tables that follow the pattern "project_datatype"
  tables <- tables[grepl("^[^_]+_[^_]+", tables)]

  # Apply user pattern filter if provided
  if (!is.null(pattern)) {
    tables <- grep(pattern, tables, value = TRUE)
  }

  if (length(tables) == 0) {
    message("No omics or drug tables found")
    return(data.frame())
  }

  # Get metadata for each table
  if ("droma_metadata" %in% all_tables) {
    metadata <- DBI::dbReadTable(connection, "droma_metadata")
    result <- metadata[metadata$table_name %in% tables, ]
  } else {
    # If metadata table doesn't exist, create basic info
    result <- data.frame(
      table_name = tables,
      feature_count = sapply(tables, function(t) {
        tryCatch({
          DBI::dbGetQuery(connection, paste0("SELECT COUNT(*) FROM ", t))[1,1]
        }, error = function(e) NA)
      }),
      sample_count = sapply(tables, function(t) {
        tryCatch({
          # Get column count minus feature_id column if it exists
          columns_query <- paste0("PRAGMA table_info(", t, ")")
          columns_info <- DBI::dbGetQuery(connection, columns_query)
          feature_id_cols <- sum(columns_info$name == "feature_id")
          max(0, nrow(columns_info) - feature_id_cols)
        }, error = function(e) NA)
      }),
      stringsAsFactors = FALSE
    )
  }

  # Add categorization by data type and feature type
  result$data_type <- sapply(result$table_name, function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 1) parts[1] else NA
  })

  # Handle feature type with special cases for mutation tables
  result$feature_type <- sapply(result$table_name, function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 2) {
      # Special handling for mutation tables
      if (length(parts) >= 3 && parts[2] == "mutation") {
        paste(parts[2:3], collapse = "_")  # "mutation_site" or "mutation_gene"
      } else {
        parts[2]  # Regular case
      }
    } else {
      NA
    }
  })

  # Add created_date from projects table if available
  result$created_date <- NA_character_
  result$updated_date <- NA_character_
  if ("projects" %in% all_tables) {
    projects_data <- DBI::dbReadTable(connection, "projects")
    for (i in 1:nrow(result)) {
      project_match <- which(projects_data$project_name == result$data_type[i])
      if (length(project_match) > 0) {
        result$created_date[i] <- projects_data$created_date[project_match[1]]
        result$updated_date[i] <- projects_data$updated_date[project_match[1]]
      }
    }
  }

  # Reorder columns for better readability
  column_order <- c("table_name", "data_type", "feature_type",
  "feature_count", "sample_count", "created_date", "updated_date")
  result <- result[, column_order[column_order %in% colnames(result)]]

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

  # Extract project names from table prefixes, excluding backup tables
  project_names <- unique(sapply(tables, function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 2 &&
        !t %in% c("sample_anno", "drug_anno") &&
        !grepl("_raw$", t)) {
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
      # Get all tables for this project, excluding backup tables
      project_tables <- grep(paste0("^", project_data_types, "_"), tables, value = TRUE)
      # Remove backup tables like "_raw"
      project_tables <- project_tables[!grepl("_raw$", project_tables)]
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


#' Update DROMA Projects Metadata
#'
#' @description Updates or adds project metadata to the projects table in the DROMA database.
#' This function automatically detects project information from existing tables and updates
#' the projects table accordingly. If a project already exists, it updates the metadata;
#' if it's new, it adds a new entry.
#'
#' @param project_name Character, the name of the project to update (e.g., "gCSI", "CCLE").
#'   If NULL, updates all projects found in the database.
#' @param dataset_type Character, the dataset type to assign to the project(s) (e.g., "CellLine", "PDX", "PDO").
#'   If NULL, attempts to guess from sample_anno table (default: NULL)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param create_table Logical, whether to create the projects table if it doesn't exist (default: TRUE)
#' @return Invisibly returns TRUE if successful
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Add some data tables first
#' updateDROMADatabase(expr_data, "gCSI_mRNA", overwrite = TRUE)
#' updateDROMADatabase(cnv_data, "gCSI_cnv", overwrite = TRUE)
#'
#' # Update project metadata for gCSI with specific dataset type
#' updateDROMAProjects("gCSI", dataset_type = "CellLine")
#'
#' # Update all projects in the database (auto-detect dataset type)
#' updateDROMAProjects()
#' }
updateDROMAProjects <- function(project_name = NULL, dataset_type = NULL, connection = NULL, create_table = TRUE) {
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

  # Get all tables in the database
  all_tables <- DBI::dbListTables(connection)

  # Extract project names from table names
  if (is.null(project_name)) {
    # Get all project names from tables, excluding backup tables
    project_names <- unique(sapply(all_tables, function(t) {
      parts <- strsplit(t, "_")[[1]]
      # Exclude system tables and backup tables (like "_raw")
      if (length(parts) >= 2 &&
          !t %in% c("sample_anno", "drug_anno", "search_vectors", "projects", "droma_metadata") &&
          !grepl("_raw$", t)) {
        return(parts[1])
      } else {
        return(NA)
      }
    }))
    project_names <- project_names[!is.na(project_names)]
  } else {
    # Use specified project name
    project_names <- project_name
  }

  if (length(project_names) == 0) {
    message("No projects found in database")
    return(invisible(FALSE))
  }

  # Check if projects table exists, create if needed
  if (!"projects" %in% all_tables) {
    if (create_table) {
      # Create projects table
      create_projects_query <- "
        CREATE TABLE projects (
          project_name TEXT PRIMARY KEY,
          dataset_type TEXT,
          data_types TEXT,
          sample_count INTEGER,
          drug_count INTEGER,
          created_date TEXT,
          updated_date TEXT
        )"
      DBI::dbExecute(connection, create_projects_query)
      message("Created projects table")
    } else {
      stop("Projects table does not exist. Set create_table = TRUE to create it.")
    }
  }

  # Get existing projects data
  existing_projects <- DBI::dbReadTable(connection, "projects")

  # Process each project
  updated_count <- 0
  added_count <- 0

  for (proj in project_names) {
    # Get tables for this project, excluding backup tables
    project_tables <- grep(paste0("^", proj, "_"), all_tables, value = TRUE)
    # Remove backup tables like "_raw"
    project_tables <- project_tables[!grepl("_raw$", project_tables)]

    if (length(project_tables) == 0) {
      warning("No tables found for project '", proj, "'")
      next
    }

    # Extract data types from table names
    data_types <- unique(sapply(project_tables, function(t) {
      sub(paste0("^", proj, "_"), "", t)
    }))

    # Determine dataset type - use user-provided or guess from sample_anno
    current_dataset_type <- dataset_type  # Use user-provided dataset_type
    sample_count <- 0

    if ("sample_anno" %in% all_tables) {
      # Get sample IDs for this project from all project tables
      sample_ids <- c()

      for (table in project_tables) {
        tryCatch({
          # Get column names from the table
          columns_query <- paste0("PRAGMA table_info(", table, ")")
          columns_info <- DBI::dbGetQuery(connection, columns_query)
          table_columns <- columns_info$name[columns_info$name != "feature_id"]

          if (length(table_columns) > 0) {
            sample_ids <- c(sample_ids, table_columns)
          }
        }, error = function(e) {
          # Skip tables that can't be processed
        })
      }

      sample_ids <- unique(sample_ids)
      sample_count <- length(sample_ids)

      # If dataset_type is NULL, try to guess from sample_anno
      if (is.null(current_dataset_type) && length(sample_ids) > 0) {
        tryCatch({
          # Use parameterized query to avoid SQL injection
          placeholders <- paste(rep("?", length(sample_ids)), collapse = ",")
          types_query <- paste0(
            "SELECT DISTINCT DataType FROM sample_anno WHERE SampleID IN (",
            placeholders, ") AND ProjectID = ?"
          )
          params <- c(as.list(sample_ids), proj)
          types_result <- DBI::dbGetQuery(connection, types_query, params = params)

          if (nrow(types_result) > 0) {
            current_dataset_type <- types_result$DataType[1]
          }
        }, error = function(e) {
          # If sample_anno query fails, continue without dataset type
        })
      }
    }

    # Convert to NA_character_ if still NULL
    if (is.null(current_dataset_type)) {
      current_dataset_type <- NA_character_
    }

    # Count drugs in drug table if available
    drug_count <- 0
    drug_table <- paste0(proj, "_drug")
    if (drug_table %in% all_tables) {
      tryCatch({
        drug_count_query <- paste0("SELECT COUNT(*) as count FROM ", drug_table)
        drug_count_result <- DBI::dbGetQuery(connection, drug_count_query)
        drug_count <- drug_count_result$count[1]
      }, error = function(e) {
        # If drug count query fails, keep drug_count as 0
      })
    }

    # Prepare project metadata
    current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    data_types_str <- paste(sort(data_types), collapse = ",")

    # Check if project already exists
    if (proj %in% existing_projects$project_name) {
      # Update existing project
      update_query <- "
        UPDATE projects SET
          dataset_type = ?,
          data_types = ?,
          sample_count = ?,
          drug_count = ?,
          updated_date = ?
        WHERE project_name = ?"

      DBI::dbExecute(connection, update_query, params = list(
        current_dataset_type, data_types_str, sample_count, drug_count, current_time, proj
      ))

      message("Updated project '", proj, "' with ", length(data_types), " data types (",
              data_types_str, "), ", sample_count, " samples, ", drug_count, " drugs")
      updated_count <- updated_count + 1

    } else {
      # Add new project
      insert_query <- "
        INSERT INTO projects (project_name, dataset_type, data_types, sample_count, drug_count, created_date, updated_date)
        VALUES (?, ?, ?, ?, ?, ?, ?)"

      DBI::dbExecute(connection, insert_query, params = list(
        proj, current_dataset_type, data_types_str, sample_count, drug_count, current_time, current_time
      ))

      message("Added new project '", proj, "' with ", length(data_types), " data types (",
              data_types_str, "), ", sample_count, " samples, ", drug_count, " drugs")
      added_count <- added_count + 1
    }
  }

  # Summary message
  if (added_count > 0 || updated_count > 0) {
    message("Project metadata update complete: ", added_count, " projects added, ",
            updated_count, " projects updated")
  } else {
    message("No projects were added or updated")
  }

  invisible(TRUE)
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

#' Update DROMA Annotation Tables with New Names
#'
#' @description Adds harmonized sample or drug names to the corresponding annotation tables.
#' Processes all entries from the name_mapping dataframe. For high confidence matches,
#' uses existing annotation information when available. Creates new entries with new_name
#' as the primary identifier and original_name recorded in ProjectRawName. Updates existing
#' entries if they already exist in the database.
#'
#' @param anno_type Character, type of annotation to update: "sample" or "drug"
#' @param name_mapping Data frame with name mappings containing at minimum: original_name, new_name, match_confidence columns
#' @param project_name Character, the project name to assign to new entries
#' @param data_type Character or character vector, the data type(s) to assign to new sample entries (e.g., "CellLine", "PDC", "PDX"). Can be a single value (applied to all) or a vector matching the length of name_mapping (only for sample annotations)
#' @param tumor_type Character or character vector, the tumor type(s) to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param PatientID Character or character vector, the patient ID(s) to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param Gender Character or character vector, the gender(s) to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param Age Numeric or numeric vector, the age(s) to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param FullEthnicity Character or character vector, the full ethnicity/ethnicities to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param SimpleEthnicity Character or character vector, the simple ethnicity/ethnicities to assign to new sample entries. Can be a single value (applied to all) or a vector matching the length of name_mapping (default: NA, only for sample annotations)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return Invisibly returns TRUE if successful, along with a summary of changes
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Check and harmonize sample names
#' sample_mapping <- checkDROMASampleNames(colnames(my_data))
#'
#' # Update sample annotations with single values
#' updateDROMAAnnotation("sample", sample_mapping, project_name = "MyProject",
#'                      data_type = "CellLine", tumor_type = "breast cancer",
#'                      PatientID = "Patient_001", Gender = "Female", Age = 45)
#'
#' # Update sample annotations with vectors (one value per sample)
#' updateDROMAAnnotation("sample", sample_mapping, project_name = "MyProject",
#'                      data_type = c("CellLine", "PDX", "CellLine", "PDO"),
#'                      tumor_type = c("breast cancer", "lung cancer", "breast cancer", "colon cancer"),
#'                      PatientID = c("Patient_001", "Patient_002", "Patient_003", "Patient_004"),
#'                      Gender = c("Female", "Male", "Female", "Male"),
#'                      Age = c(45, 52, 38, 41),
#'                      FullEthnicity = c("European", "Asian", "African", "Hispanic"),
#'                      SimpleEthnicity = c("Caucasian", "Asian", "African", "Hispanic"))
#'
#' # Check and harmonize drug names
#' drug_mapping <- checkDROMADrugNames(rownames(my_drug_data))
#'
#' # Update drug annotations (PatientID not used for drugs)
#' updateDROMAAnnotation("drug", drug_mapping, project_name = "MyProject")
#' }
updateDROMAAnnotation <- function(anno_type, name_mapping, project_name, data_type = NA_character_,
                                  tumor_type = NA_character_, PatientID = NA_character_, Gender = NA_character_,
                                  Age = NA_real_, FullEthnicity = NA_character_,
                                  SimpleEthnicity = NA_character_, connection = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Validate anno_type
  if (!anno_type %in% c("sample", "drug")) {
    stop("anno_type must be either 'sample' or 'drug'")
  }

  # Validate Age parameter - must be numeric when provided
  if (!all(is.na(Age)) && !is.numeric(Age)) {
    stop("Age parameter must be numeric (single value or vector) or NA")
  }

  # If Age is a vector, validate length
  if (length(Age) > 1 && length(Age) != nrow(name_mapping)) {
    stop("Age vector length (", length(Age), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
  }

  # Validate character parameters for sample annotations
  if (anno_type == "sample") {
    # Validate data_type parameter
    if (!all(is.na(data_type)) && !is.character(data_type)) {
      stop("data_type parameter must be character (single value or vector) or NA")
    }
    if (length(data_type) > 1 && length(data_type) != nrow(name_mapping)) {
      stop("data_type vector length (", length(data_type), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }

    # Validate tumor_type parameter
    if (!all(is.na(tumor_type)) && !is.character(tumor_type)) {
      stop("tumor_type parameter must be character (single value or vector) or NA")
    }
    if (length(tumor_type) > 1 && length(tumor_type) != nrow(name_mapping)) {
      stop("tumor_type vector length (", length(tumor_type), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }

    # Validate PatientID parameter
    if (!all(is.na(PatientID)) && !is.character(PatientID)) {
      stop("PatientID parameter must be character (single value or vector) or NA")
    }
    if (length(PatientID) > 1 && length(PatientID) != nrow(name_mapping)) {
      stop("PatientID vector length (", length(PatientID), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }

    # Validate Gender parameter
    if (!all(is.na(Gender)) && !is.character(Gender)) {
      stop("Gender parameter must be character (single value or vector) or NA")
    }
    if (length(Gender) > 1 && length(Gender) != nrow(name_mapping)) {
      stop("Gender vector length (", length(Gender), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }

    # Validate FullEthnicity parameter
    if (!all(is.na(FullEthnicity)) && !is.character(FullEthnicity)) {
      stop("FullEthnicity parameter must be character (single value or vector) or NA")
    }
    if (length(FullEthnicity) > 1 && length(FullEthnicity) != nrow(name_mapping)) {
      stop("FullEthnicity vector length (", length(FullEthnicity), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }

    # Validate SimpleEthnicity parameter
    if (!all(is.na(SimpleEthnicity)) && !is.character(SimpleEthnicity)) {
      stop("SimpleEthnicity parameter must be character (single value or vector) or NA")
    }
    if (length(SimpleEthnicity) > 1 && length(SimpleEthnicity) != nrow(name_mapping)) {
      stop("SimpleEthnicity vector length (", length(SimpleEthnicity), ") must match name_mapping rows (", nrow(name_mapping), ") or be a single value")
    }
  }

  # Validate name_mapping structure
  required_cols <- c("original_name", "new_name", "match_confidence")
  if (!all(required_cols %in% colnames(name_mapping))) {
    stop("name_mapping must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Determine table name and column names
  if (anno_type == "sample") {
    table_name <- "sample_anno"
    id_column <- "SampleID"
    raw_column <- "ProjectRawName"
    index_prefix <- "UM_SAMPLE_"
  } else {
    table_name <- "drug_anno"
    id_column <- "DrugName"
    raw_column <- "ProjectRawName"
    index_prefix <- "UM_DRUG_"
  }

  # Check if table exists
  all_tables <- DBI::dbListTables(connection)
  if (!table_name %in% all_tables) {
    stop("Annotation table '", table_name, "' not found in database")
  }

  # Get existing annotation data
  existing_anno <- DBI::dbReadTable(connection, table_name)

  # Find the maximum IndexID number for generating new IDs
  max_index_num <- 0
  if ("IndexID" %in% colnames(existing_anno) && nrow(existing_anno) > 0) {
    # Extract numbers from IndexID column
    index_ids <- existing_anno$IndexID[!is.na(existing_anno$IndexID)]
    if (length(index_ids) > 0) {
      # Extract numbers from IndexID strings like "UM_SAMPLE_1" or "UM_DRUG_2077"
      numbers <- as.numeric(gsub(paste0("^", index_prefix), "", index_ids))
      numbers <- numbers[!is.na(numbers)]
      if (length(numbers) > 0) {
        max_index_num <- max(numbers)
      }
    }
  }

  # Initialize counters
  added_count <- 0
  skipped_count <- 0
  current_index_num <- max_index_num

  # Process ALL entries from name_mapping - add each as a new entry
  for (i in 1:nrow(name_mapping)) {
    new_name <- name_mapping$new_name[i]
    original_name <- name_mapping$original_name[i]
    match_confidence <- name_mapping$match_confidence[i]

    # Get current parameter values (handle both single values and vectors)
    current_age <- if (length(Age) == 1) Age else Age[i]
    current_data_type <- if (length(data_type) == 1) data_type else data_type[i]
    current_tumor_type <- if (length(tumor_type) == 1) tumor_type else tumor_type[i]
    current_patient_id <- if (length(PatientID) == 1) PatientID else PatientID[i]
    current_gender <- if (length(Gender) == 1) Gender else Gender[i]
    current_full_ethnicity <- if (length(FullEthnicity) == 1) FullEthnicity else FullEthnicity[i]
    current_simple_ethnicity <- if (length(SimpleEthnicity) == 1) SimpleEthnicity else SimpleEthnicity[i]

    # Check if this SampleID/DrugName + ProjectID combination already exists
    if (anno_type == "sample") {
      existing_check <- which(existing_anno$SampleID == new_name & existing_anno$ProjectID == project_name)
    } else {
      existing_check <- which(existing_anno$DrugName == new_name & existing_anno$ProjectID == project_name)
    }

    if (length(existing_check) > 0) {
      # Skip this entry as it already exists
      skipped_count <- skipped_count + 1
      next
    }

    # Determine IndexID for this entry
    new_index_id <- NA_character_
    if (match_confidence == "high") {
      # For high confidence matches, try to use existing IndexID if available
      if (anno_type == "sample") {
        original_entry_idx <- which(existing_anno$SampleID == new_name |
                                      existing_anno$ProjectRawName == new_name)
      } else {
        original_entry_idx <- which(existing_anno$DrugName == new_name |
                                      existing_anno$ProjectRawName == new_name)
      }

      if (length(original_entry_idx) > 0) {
        original_entry <- existing_anno[original_entry_idx[1], ]
        if ("IndexID" %in% colnames(original_entry) && !is.na(original_entry$IndexID)) {
          new_index_id <- original_entry$IndexID
        }
      }
    }

    # If no existing IndexID found (or not high confidence), generate new one
    if (is.na(new_index_id)) {
      current_index_num <- current_index_num + 1
      new_index_id <- paste0(index_prefix, current_index_num)
    }

    # Add new entry - behavior depends on match_confidence
    if (match_confidence == "high") {
      # For high confidence matches, check if we can find existing info from the harmonized name
      if (anno_type == "sample") {
        original_entry_idx <- which(existing_anno$SampleID == name_mapping$harmonized_name[i] |
                                      existing_anno$ProjectRawName == name_mapping$harmonized_name[i])
      } else {
        original_entry_idx <- which(existing_anno$DrugName == name_mapping$harmonized_name[i] |
                                      existing_anno$ProjectRawName == name_mapping$harmonized_name[i])
      }

      if (length(original_entry_idx) > 0) {
        # Use existing entry info but update for this project
        # For existing entries with IndexID, use original values; user parameters only for new samples
        original_entry <- existing_anno[original_entry_idx[1], ]

        if (anno_type == "sample") {
          insert_query <- paste0("INSERT INTO ", table_name, " (",
                                 "SampleID, PatientID, ProjectID, HarmonizedIdentifier, TumorType, ",
                                 "MolecularSubtype, Gender, Age, FullEthnicity, SimpleEthnicity, ",
                                 "TNMstage, Primary_Metastasis, DataType, ProjectRawName, AlternateName, IndexID) ",
                                 "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
          DBI::dbExecute(connection, insert_query, params = list(
            new_name,
            if("PatientID" %in% colnames(original_entry)) original_entry$PatientID else current_patient_id,
            project_name,
            if("HarmonizedIdentifier" %in% colnames(original_entry)) original_entry$HarmonizedIdentifier else NA_character_,
            # Use original values from existing annotation
            if("TumorType" %in% colnames(original_entry)) original_entry$TumorType else NA_character_,
            if("MolecularSubtype" %in% colnames(original_entry)) original_entry$MolecularSubtype else NA_character_,
            if("Gender" %in% colnames(original_entry)) original_entry$Gender else NA_character_,
            if("Age" %in% colnames(original_entry)) original_entry$Age else NA_character_,
            if("FullEthnicity" %in% colnames(original_entry)) original_entry$FullEthnicity else NA_character_,
            if("SimpleEthnicity" %in% colnames(original_entry)) original_entry$SimpleEthnicity else NA_character_,
            if("TNMstage" %in% colnames(original_entry)) original_entry$TNMstage else NA_character_,
            if("Primary_Metastasis" %in% colnames(original_entry)) original_entry$Primary_Metastasis else NA_character_,
            # Use original data_type from existing annotation
            if("DataType" %in% colnames(original_entry)) original_entry$DataType else current_data_type,
            original_name,
            if("AlternateName" %in% colnames(original_entry)) original_entry$AlternateName else NA_character_,
            new_index_id
          ))
    } else {
          insert_query <- paste0("INSERT INTO ", table_name, " (",
                                 "DrugName, ProjectID, `Harmonized ID (Pubchem ID)`, `Source for Clinical Information`, ",
                                 "`Clinical Phase`, MOA, Targets, ProjectRawName, IndexID) ",
                                 "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")
      DBI::dbExecute(connection, insert_query, params = list(
            new_name,
            project_name,
            if("Harmonized ID (Pubchem ID)" %in% colnames(original_entry)) original_entry$`Harmonized ID (Pubchem ID)` else NA_character_,
            if("Source for Clinical Information" %in% colnames(original_entry)) original_entry$`Source for Clinical Information` else NA_character_,
            if("Clinical Phase" %in% colnames(original_entry)) original_entry$`Clinical Phase` else NA_character_,
            if("MOA" %in% colnames(original_entry)) original_entry$MOA else NA_character_,
            if("Targets" %in% colnames(original_entry)) original_entry$Targets else NA_character_,
            original_name,
            new_index_id
          ))
        }
      added_count <- added_count + 1
      } else {
        # Couldn't find original entry, create new one using user parameters
        if (anno_type == "sample") {
          insert_query <- paste0("INSERT INTO ", table_name, " (",
                                 "SampleID, PatientID, ProjectID, HarmonizedIdentifier, TumorType, ",
                                 "MolecularSubtype, Gender, Age, FullEthnicity, SimpleEthnicity, ",
                                 "TNMstage, Primary_Metastasis, DataType, ProjectRawName, AlternateName, IndexID) ",
                                 "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
                  DBI::dbExecute(connection, insert_query, params = list(
          new_name, current_patient_id, project_name, NA_character_, current_tumor_type,
          NA_character_, current_gender, if(is.na(current_age)) NA_character_ else as.character(current_age), current_full_ethnicity, current_simple_ethnicity,
          NA_character_, NA_character_, current_data_type, original_name, NA_character_, new_index_id
        ))
        } else {
          insert_query <- paste0("INSERT INTO ", table_name, " (",
                                 "DrugName, ProjectID, `Harmonized ID (Pubchem ID)`, `Source for Clinical Information`, ",
                                 "`Clinical Phase`, MOA, Targets, ProjectRawName, IndexID) ",
                                 "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")
          DBI::dbExecute(connection, insert_query, params = list(
            new_name, project_name, NA_character_, NA_character_,
            NA_character_, NA_character_, NA_character_, original_name, new_index_id
          ))
        }
        added_count <- added_count + 1
      }
  } else {
      # For non-high confidence matches, create new entry with new_name as identifier using user parameters
      if (anno_type == "sample") {
        insert_query <- paste0("INSERT INTO ", table_name, " (",
                               "SampleID, PatientID, ProjectID, HarmonizedIdentifier, TumorType, ",
                               "MolecularSubtype, Gender, Age, FullEthnicity, SimpleEthnicity, ",
                               "TNMstage, Primary_Metastasis, DataType, ProjectRawName, AlternateName, IndexID) ",
                               "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
        DBI::dbExecute(connection, insert_query, params = list(
          new_name, current_patient_id, project_name, NA_character_, current_tumor_type,
          NA_character_, current_gender, if(is.na(current_age)) NA_character_ else as.character(current_age), current_full_ethnicity, current_simple_ethnicity,
          NA_character_, NA_character_, current_data_type, original_name, NA_character_, new_index_id
        ))
      } else {
        insert_query <- paste0("INSERT INTO ", table_name, " (",
                               "DrugName, ProjectID, `Harmonized ID (Pubchem ID)`, `Source for Clinical Information`, ",
                               "`Clinical Phase`, MOA, Targets, ProjectRawName, IndexID) ",
                               "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")
        DBI::dbExecute(connection, insert_query, params = list(
          new_name, project_name, NA_character_, NA_character_,
          NA_character_, NA_character_, NA_character_, original_name, new_index_id
        ))
      }
      added_count <- added_count + 1
    }
  }

  # Update created_date if it's NA for this project
  if (added_count > 0) {
    all_tables <- DBI::dbListTables(connection)
    if ("projects" %in% all_tables) {
      # Check if project exists and has NA created_date
      check_query <- "SELECT created_date FROM projects WHERE project_name = ?"
      existing_project <- DBI::dbGetQuery(connection, check_query, params = list(project_name))

      if (nrow(existing_project) > 0 && is.na(existing_project$created_date[1])) {
        # Update created_date to current time
        current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        update_query <- "UPDATE projects SET created_date = ? WHERE project_name = ?"
        DBI::dbExecute(connection, update_query, params = list(current_time, project_name))
        message("  Updated created_date for project '", project_name, "' to ", current_time)
      }
    }
  }

  # Print summary
  message("Updated ", table_name, " table:")
  message("  Added: ", added_count, " new entries")
  message("  Skipped: ", skipped_count, " existing entries")
  if (current_index_num > max_index_num) {
    message("  Generated new IndexIDs from ", index_prefix, max_index_num + 1, " to ", index_prefix, current_index_num)
  }

  # Print match confidence summary for all entries
  if ("match_type" %in% colnames(name_mapping)) {
    match_summary <- table(name_mapping$match_type)
    message("  Match types for processed entries:")
    for (match_type in names(match_summary)) {
      message("    ", match_type, ": ", match_summary[match_type])
    }
  } else {
    confidence_summary <- table(name_mapping$match_confidence)
    message("  Match confidence for processed entries:")
    for (confidence in names(confidence_summary)) {
      message("    ", confidence, ": ", confidence_summary[confidence])
    }
  }

  # Additional note about sample-specific parameters
  if (anno_type == "sample" && (added_count > 0)) {
    message("  Note: For samples with existing IndexID, original annotation values were preserved.")
    message("        User-provided parameters (PatientID, tumor_type, Gender, Age, FullEthnicity, SimpleEthnicity)")
    message("        were only applied to new sample entries without existing records.")
  }

  invisible(TRUE)
}


#' Check and Harmonize Sample Names Against DROMA Database
#'
#' @description Checks sample names (column names) against the sample_anno table in the DROMA database
#' and provides harmonized mappings. Uses fuzzy matching and name cleaning approach.
#'
#' @param sample_names Character vector of sample names to check and harmonize
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param max_distance Numeric, maximum distance for fuzzy matching (default: 0.2)
#' @param min_name_length Integer, minimum name length for partial matching (default: 5)
#' @return A data frame with columns: original_name, cleaned_name, harmonized_name, match_type, match_confidence, new_name
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Check sample names from a data matrix
#' sample_names <- colnames(my_data_matrix)
#' name_mapping <- checkDROMASampleNames(sample_names)
#'
#' # View the mapping results
#' print(name_mapping)
#' }
checkDROMASampleNames <- function(sample_names, connection = NULL, max_distance = 0.2, min_name_length = 5) {
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

  # Get sample annotation data
  sample_anno <- DBI::dbReadTable(connection, "sample_anno")

  # Function to clean sample names for better matching
  clean_sample_name <- function(name) {
    # Convert to lowercase
    name <- tolower(name)

    # Remove [?]
    name <- gsub("[?]", "", name, fixed = TRUE)

    # Remove Chinese characters (if any)
    name <- gsub("[\\p{Han}]", "", name, perl = TRUE)

    # First remove [xx] format and its contents
    name <- gsub("\\[.*?\\]", "", name)

    # Handle parentheses - only remove if there's content outside them
    if (!grepl("^\\s*\\(.*\\)\\s*$", name)) {
      name <- gsub("\\s*\\([^\\)]+\\)", "", name)
    } else {
      # For names entirely in parentheses, remove the parentheses but keep content
      name <- gsub("^\\s*\\(|\\)\\s*$", "", name)
    }

    # Remove special characters and extra spaces
    name <- gsub("[^a-z0-9]", "", name)
    name <- trimws(name)

    return(name)
  }

  # Clean reference names in sample_anno
  sample_anno$clean_sampleid <- sapply(sample_anno$SampleID, clean_sample_name)
  sample_anno$clean_rawname <- sapply(sample_anno$ProjectRawName, clean_sample_name)

  # Handle AlternateName column if it exists and create alternate name mapping
  alternate_mapping <- data.frame()
  if ("AlternateName" %in% colnames(sample_anno)) {
    for (i in 1:nrow(sample_anno)) {
      harmonized_name <- sample_anno$SampleID[i]

      # Add the harmonized name itself
      alternate_mapping <- rbind(alternate_mapping, data.frame(
        raw_name = harmonized_name,
        harmonized_name = harmonized_name,
        clean_name = clean_sample_name(harmonized_name),
        stringsAsFactors = FALSE
      ))

      # Process alternate names if they exist
      if (!is.na(sample_anno$AlternateName[i]) && sample_anno$AlternateName[i] != "") {
        alt_names <- strsplit(sample_anno$AlternateName[i], ":|:")[[1]]
        alt_names <- alt_names[!alt_names %in% "|" & alt_names != ""]
        # Add each alternate name
        for (alt_name in alt_names) {
          alternate_mapping <- rbind(alternate_mapping, data.frame(
            raw_name = alt_name,
            harmonized_name = harmonized_name,
            clean_name = clean_sample_name(alt_name),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Clean input sample names
  sample_names_clean <- sapply(sample_names, clean_sample_name)

  # Create result data frame
  result <- data.frame(
    original_name = sample_names,
    cleaned_name = sample_names_clean,
    harmonized_name = NA_character_,
    match_type = NA_character_,
    match_confidence = NA_character_,
    stringsAsFactors = FALSE
  )

  # Match each sample name
  for (i in 1:length(sample_names)) {
    clean_name <- sample_names_clean[i]
    original_name <- sample_names[i]

    # For very long names, keep original
    if (nchar(original_name) > 30) {
      result$harmonized_name[i] <- original_name
      result$match_type[i] <- "keep_original_long"
      result$match_confidence[i] <- "medium"
      next
    }

    # Try exact match with SampleID
    exact_sampleid <- which(sample_anno$clean_sampleid == clean_name)
    if (length(exact_sampleid) > 0) {
      result$harmonized_name[i] <- sample_anno$SampleID[exact_sampleid[1]]
      result$match_type[i] <- "exact_sampleid"
      result$match_confidence[i] <- "high"
      result$new_name[i] <- sample_anno$SampleID[exact_sampleid[1]]
      next
    }

    # Try exact match with ProjectRawName
    exact_rawname <- which(sample_anno$clean_rawname == clean_name)
    if (length(exact_rawname) > 0) {
      result$harmonized_name[i] <- sample_anno$SampleID[exact_rawname[1]]
      result$match_type[i] <- "exact_rawname"
      result$match_confidence[i] <- "high"
      result$new_name[i] <- sample_anno$SampleID[exact_rawname[1]]
      next
    }

    # Try exact match with AlternateName if available
    if (nrow(alternate_mapping) > 0) {
      exact_alternate <- which(alternate_mapping$clean_name == clean_name)
      if (length(exact_alternate) > 0) {
        result$harmonized_name[i] <- alternate_mapping$harmonized_name[exact_alternate[1]]
        result$match_type[i] <- "exact_alternate"
        result$match_confidence[i] <- "high"
        result$new_name[i] <- alternate_mapping$harmonized_name[exact_alternate[1]]
        next
      }
    }

    # Try fuzzy match with SampleID
    if (nchar(clean_name) >= 3) {
      fuzzy_sampleid <- agrep(clean_name, sample_anno$clean_sampleid, max.distance = max_distance, value = FALSE)
      if (length(fuzzy_sampleid) > 0) {
        result$harmonized_name[i] <- sample_anno$SampleID[fuzzy_sampleid[1]]
        result$match_type[i] <- "fuzzy_sampleid"
        result$match_confidence[i] <- "medium"
        result$new_name[i] <- sample_anno$SampleID[fuzzy_sampleid[1]]
        next
      }
    }

    # Try fuzzy match with ProjectRawName
    if (nchar(clean_name) >= 3) {
      fuzzy_rawname <- agrep(clean_name, sample_anno$clean_rawname, max.distance = max_distance, value = FALSE)
      if (length(fuzzy_rawname) > 0) {
        result$harmonized_name[i] <- sample_anno$SampleID[fuzzy_rawname[1]]
        result$match_type[i] <- "fuzzy_rawname"
        result$match_confidence[i] <- "medium"
        result$new_name[i] <- sample_anno$SampleID[fuzzy_rawname[1]]
        next
      }
    }

    # Try fuzzy match with AlternateName if available
    if (nrow(alternate_mapping) > 0 && nchar(clean_name) >= 3) {
      fuzzy_alternate <- agrep(clean_name, alternate_mapping$clean_name, max.distance = max_distance, value = FALSE)
      if (length(fuzzy_alternate) > 0) {
        result$harmonized_name[i] <- alternate_mapping$harmonized_name[fuzzy_alternate[1]]
        result$match_type[i] <- "fuzzy_alternate"
        result$match_confidence[i] <- "medium"
        result$new_name[i] <- alternate_mapping$harmonized_name[fuzzy_alternate[1]]
        next
      }
    }

    # Try partial match (sample name is contained in annotation)
    if (nchar(clean_name) >= min_name_length) {
      # Check if sample name is contained in SampleID names
      partial_sampleid <- which(sapply(sample_anno$clean_sampleid, function(x) grepl(clean_name, x, fixed = TRUE)))
      if (length(partial_sampleid) > 0) {
        result$harmonized_name[i] <- sample_anno$SampleID[partial_sampleid[1]]
        result$match_type[i] <- "partial_sampleid"
        result$match_confidence[i] <- "low"
        result$new_name[i] <- sample_anno$SampleID[partial_sampleid[1]]
        next
      }

      # Check if sample name is contained in ProjectRawName names
      partial_rawname <- which(sapply(sample_anno$clean_rawname, function(x) grepl(clean_name, x, fixed = TRUE)))
      if (length(partial_rawname) > 0) {
        result$harmonized_name[i] <- sample_anno$SampleID[partial_rawname[1]]
        result$match_type[i] <- "partial_rawname"
        result$match_confidence[i] <- "low"
        result$new_name[i] <- sample_anno$SampleID[partial_rawname[1]]
        next
      }
    }

    # No match found - use cleaned name
    result$harmonized_name[i] <- clean_name
    result$match_type[i] <- "no_match"
    result$match_confidence[i] <- "none"
  }

  # Set new_name based on match_confidence
  result$new_name <- ifelse(result$match_confidence == "high",
                           result$harmonized_name,
                           result$original_name)

  # Print summary
  match_summary <- table(result$match_type)
  message("Sample name matching summary:")
  for (match_type in names(match_summary)) {
    message("    ", match_type, ": ", match_summary[match_type])
  }

  # Warn about low confidence matches
  low_confidence <- result[result$match_confidence %in% c("medium", "low"), ]
  if (nrow(low_confidence) > 0) {
    message("\nWarning: ", nrow(low_confidence), " samples have medium/low confidence matches.")
    message("Consider manual review of these matches:")
    for (i in 1:min(5, nrow(low_confidence))) {
      message("    ", low_confidence$original_name[i], " -> ", low_confidence$harmonized_name[i],
              " (", low_confidence$match_type[i], ")")
    }
    if (nrow(low_confidence) > 5) {
      message("    ... and ", nrow(low_confidence) - 5, " more")
    }
  }

  return(result)
}

#' Check and Harmonize Drug Names Against DROMA Database
#'
#' @description Checks drug names (row names) against the drug_anno table in the DROMA database
#' and provides harmonized mappings. Uses fuzzy matching and name cleaning approach.
#'
#' @param drug_names Character vector of drug names to check and harmonize
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param max_distance Numeric, maximum distance for fuzzy matching (default: 0.2)
#' @param min_name_length Integer, minimum name length for partial matching (default: 5)
#' @param keep_long_names_threshold Integer, names longer than this will be kept as original (default: 17)
#' @return A data frame with columns: original_name, cleaned_name, harmonized_name, match_type, match_confidence, new_name
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Check drug names from a drug response matrix
#' drug_names <- rownames(my_drug_matrix)
#' name_mapping <- checkDROMADrugNames(drug_names)
#'
#' # View the mapping results
#' print(name_mapping)
#' }
checkDROMADrugNames <- function(drug_names, connection = NULL, max_distance = 0.2,
                               min_name_length = 5, keep_long_names_threshold = 17) {
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

  # Check if drug_anno table exists
  all_tables <- DBI::dbListTables(connection)
  if (!"drug_anno" %in% all_tables) {
    stop("Drug annotation table 'drug_anno' not found in database")
  }

  # Get drug annotation data
  drug_anno <- DBI::dbReadTable(connection, "drug_anno")

  # Function to clean drug names for better matching
  clean_drug_name <- function(name) {
    # Convert to lowercase
    name <- tolower(name)

    # Remove [?]
    name <- gsub("[?]", "", name, fixed = TRUE)

    # Remove Chinese characters (if any)
    name <- gsub("[\\p{Han}]", "", name, perl = TRUE)

    # Handle parentheses - only remove if there's content outside them
    if (!grepl("^\\s*\\(.*\\)\\s*$", name)) {
      name <- gsub("\\s*\\([^\\)]+\\)", "", name)
    } else {
      # For names entirely in parentheses, remove the parentheses but keep content
      name <- gsub("^\\s*\\(|\\)\\s*$", "", name)
    }

    # Remove special characters and extra spaces
    name <- gsub("[^a-z0-9]", " ", name)
    name <- gsub("\\s+", " ", name)
    name <- trimws(name)

    return(name)
  }

  # Clean reference names in drug_anno
  drug_anno$clean_drugname <- sapply(drug_anno$DrugName, clean_drug_name)
  drug_anno$clean_rawname <- sapply(drug_anno$ProjectRawName, clean_drug_name)

  # Clean input drug names
  drug_names_clean <- sapply(drug_names, clean_drug_name)

  # Create result data frame
  result <- data.frame(
    original_name = drug_names,
    cleaned_name = drug_names_clean,
    harmonized_name = NA_character_,
    match_type = NA_character_,
    match_confidence = NA_character_,
    stringsAsFactors = FALSE
  )

  # Match each drug name
  for (i in 1:length(drug_names)) {
    clean_name <- drug_names_clean[i]
    original_name <- drug_names[i]

    # For very long drug names, keep original
    if (nchar(original_name) > keep_long_names_threshold) {
      result$harmonized_name[i] <- original_name
      result$match_type[i] <- "keep_original_long"
      result$match_confidence[i] <- "medium"
      next
    }

    # Try exact match with DrugName
    exact_drugname <- which(drug_anno$clean_drugname == clean_name)
    if (length(exact_drugname) > 0) {
      result$harmonized_name[i] <- drug_anno$DrugName[exact_drugname[1]]
      result$match_type[i] <- "exact_drugname"
      result$match_confidence[i] <- "high"
      next
    }

    # Try exact match with ProjectRawName
    exact_rawname <- which(drug_anno$clean_rawname == clean_name)
    if (length(exact_rawname) > 0) {
      result$harmonized_name[i] <- drug_anno$DrugName[exact_rawname[1]]
      result$match_type[i] <- "exact_rawname"
      result$match_confidence[i] <- "high"
      next
    }

    # Try fuzzy match with DrugName
    if (nchar(clean_name) >= 3) {
      fuzzy_drugname <- agrep(clean_name, drug_anno$clean_drugname, max.distance = max_distance, value = FALSE)
      if (length(fuzzy_drugname) > 0) {
        result$harmonized_name[i] <- drug_anno$DrugName[fuzzy_drugname[1]]
        result$match_type[i] <- "fuzzy_drugname"
        result$match_confidence[i] <- "medium"
        next
      }
    }

    # Try fuzzy match with ProjectRawName
    if (nchar(clean_name) >= 3) {
      fuzzy_rawname <- agrep(clean_name, drug_anno$clean_rawname, max.distance = max_distance, value = FALSE)
      if (length(fuzzy_rawname) > 0) {
        result$harmonized_name[i] <- drug_anno$DrugName[fuzzy_rawname[1]]
        result$match_type[i] <- "fuzzy_rawname"
        result$match_confidence[i] <- "medium"
        next
      }
    }

    # Try partial match (drug name is contained in annotation)
    if (nchar(clean_name) >= min_name_length) {
      # Check if drug name is contained in DrugName names
      partial_drugname <- which(sapply(drug_anno$clean_drugname, function(x) grepl(clean_name, x, fixed = TRUE)))
      if (length(partial_drugname) > 0) {
        result$harmonized_name[i] <- drug_anno$DrugName[partial_drugname[1]]
        result$match_type[i] <- "partial_drugname"
        result$match_confidence[i] <- "low"
        next
      }

      # Check if drug name is contained in ProjectRawName names
      partial_rawname <- which(sapply(drug_anno$clean_rawname, function(x) grepl(clean_name, x, fixed = TRUE)))
      if (length(partial_rawname) > 0) {
        result$harmonized_name[i] <- drug_anno$DrugName[partial_rawname[1]]
        result$match_type[i] <- "partial_rawname"
        result$match_confidence[i] <- "low"
        next
      }
    }

    # No match found - use cleaned name
    result$harmonized_name[i] <- clean_name
    result$match_type[i] <- "no_match"
    result$match_confidence[i] <- "none"
  }

  # Set new_name based on match_confidence
  result$new_name <- ifelse(result$match_confidence == "high",
                           result$harmonized_name,
                           result$original_name)

  # Print summary
  match_summary <- table(result$match_type)
  message("Drug name matching summary:")
  for (match_type in names(match_summary)) {
    message("    ", match_type, ": ", match_summary[match_type])
  }

  # Warn about low confidence matches
  low_confidence <- result[result$match_confidence %in% c("medium", "low"), ]
  if (nrow(low_confidence) > 0) {
    message("\nWarning: ", nrow(low_confidence), " drugs have medium/low confidence matches.")
    message("Consider manual review of these matches:")
    for (i in 1:min(5, nrow(low_confidence))) {
      message("    ", low_confidence$original_name[i], " -> ", low_confidence$harmonized_name[i],
              " (", low_confidence$match_type[i], ")")
    }
    if (nrow(low_confidence) > 5) {
      message("    ... and ", nrow(low_confidence) - 5, " more")
    }
  }

  return(result)
}
