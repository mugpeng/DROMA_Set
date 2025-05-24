#' Create a MultiDromaSet from Database
#'
#' @description Creates a MultiDromaSet object containing multiple projects from a SQLite database
#' @param project_names Character vector, the names of the projects/datasets to include (e.g., c("gCSI", "CCLE"))
#' @param db_path Path to the SQLite database
#' @param db_groups Optional character vector of group names in the database, if different from project_names
#' @param load_metadata Logical, whether to load sample and treatment metadata (default: TRUE)
#' @param dataset_types Optional character vector of dataset types for each project (e.g., c("CellLine", "PDX"))
#' @param auto_load Logical, whether to automatically load treatment response and molecular profiles (default: FALSE)
#' @return A MultiDromaSet object linked to the database
#' @export
#' @examples
#' \dontrun{
#' # Create a MultiDromaSet for gCSI and CCLE data from database
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "~/droma.sqlite")
#'
#' # Create a MultiDromaSet and automatically load all data
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "~/droma.sqlite", auto_load = TRUE)
#'
#' # Create with specific dataset types
#' multi_set <- createMultiDromaSetFromDatabase(
#'   project_names = c("gCSI", "PDX_data"),
#'   db_path = "~/droma.sqlite",
#'   dataset_types = c("CellLine", "PDX")
#' )
#' }
createMultiDromaSetFromDatabase <- function(project_names,
                                          db_path = file.path(path.expand("~"), "droma.sqlite"),
                                          db_groups = NULL,
                                          load_metadata = TRUE,
                                          dataset_types = NULL,
                                          auto_load = FALSE) {

  if (!file.exists(db_path)) {
    stop("Database file not found: ", db_path)
  }

  if (length(project_names) == 0) {
    stop("At least one project name must be provided")
  }

  # Set db_groups to project_names if not specified
  if (is.null(db_groups)) {
    db_groups <- project_names
  }

  if (length(db_groups) != length(project_names)) {
    stop("Length of db_groups must match length of project_names")
  }

  # Set dataset_types if not provided
  if (is.null(dataset_types)) {
    dataset_types <- rep(NA_character_, length(project_names))
  }

  if (length(dataset_types) != length(project_names)) {
    stop("Length of dataset_types must match length of project_names")
  }

  # Connect to database to validate projects exist
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  all_tables <- DBI::dbListTables(con)

  # Validate that each project has tables in the database
  for (i in seq_along(project_names)) {
    project_tables <- grep(paste0("^", db_groups[i], "_"), all_tables, value = TRUE)
    if (length(project_tables) == 0) {
      stop("No tables found for project '", project_names[i], "' with group prefix '", db_groups[i], "'")
    }
  }

  # Create individual DromaSet objects for each project
  droma_sets <- list()

  for (i in seq_along(project_names)) {
    message("Creating DromaSet for project: ", project_names[i])

    tryCatch({
      ds <- createDromaSetFromDatabase(
        project_name = project_names[i],
        db_path = db_path,
        db_group = db_groups[i],
        load_metadata = load_metadata,
        dataset_type = dataset_types[i],
        auto_load = auto_load
      )
      droma_sets[[project_names[i]]] <- ds
    }, error = function(e) {
      warning("Problem with creating DromaSet for project '", project_names[i], "': ", e$message)
    })
  }

  if (length(droma_sets) == 0) {
    stop("No DromaSet objects could be created from the specified projects")
  }

  # Create the MultiDromaSet object
  multi_set <- MultiDromaSet(
    name = names(droma_sets),
    DromaSets = droma_sets,
    db_info = list(db_path = db_path)
  )

  message("Created MultiDromaSet with ", length(droma_sets), " projects: ",
          paste(names(droma_sets), collapse = ", "))

  return(multi_set)
}

#' Create a MultiDromaSet from Existing DromaSet Objects
#'
#' @description Creates a MultiDromaSet object from a list of existing DromaSet objects
#' @param droma_sets List of DromaSet objects or individual DromaSet objects passed as separate arguments
#' @param project_names Optional character vector of project names (will be extracted from DromaSet objects if not provided)
#' @param merge_metadata Logical, whether to merge sample and treatment metadata from all DromaSets (default: TRUE)
#' @return A MultiDromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Create individual DromaSet objects
#' gCSI <- createDromaSetFromDatabase("gCSI", "~/droma.sqlite")
#' CCLE <- createDromaSetFromDatabase("CCLE", "~/droma.sqlite")
#'
#' # Method 1: Pass as a list
#' multi_set <- createMultiDromaSetFromObjects(list(gCSI = gCSI, CCLE = CCLE))
#'
#' # Method 2: Pass as separate arguments
#' multi_set <- createMultiDromaSetFromObjects(gCSI, CCLE)
#'
#' # Method 3: Specify custom project names
#' multi_set <- createMultiDromaSetFromObjects(
#'   list(gCSI, CCLE),
#'   project_names = c("Genomics_CSI", "Cancer_Cell_Lines")
#' )
#' }
createMultiDromaSetFromObjects <- function(..., project_names = NULL, merge_metadata = TRUE) {

  # Handle different input formats
  args <- list(...)

  if (length(args) == 1 && is.list(args[[1]]) && !inherits(args[[1]], "DromaSet")) {
    # Input is a single list of DromaSet objects
    droma_sets <- args[[1]]
  } else {
    # Input is individual DromaSet objects
    droma_sets <- args
  }

  if (length(droma_sets) == 0) {
    stop("At least one DromaSet object must be provided")
  }

  # Validate that all inputs are DromaSet objects
  if (!all(sapply(droma_sets, function(x) inherits(x, "DromaSet")))) {
    stop("All inputs must be DromaSet objects")
  }

  # Extract project names if not provided
  if (is.null(project_names)) {
    if (is.null(names(droma_sets))) {
      # Extract names from DromaSet objects
      project_names <- sapply(droma_sets, function(x) x@name)
    } else {
      # Use list names
      project_names <- names(droma_sets)
    }
  }

  # Validate project names
  if (length(project_names) != length(droma_sets)) {
    stop("Length of project_names must match number of DromaSet objects")
  }

  # Set names for the list
  names(droma_sets) <- project_names

  # Check for duplicate project names
  if (any(duplicated(project_names))) {
    stop("Duplicate project names found: ", paste(project_names[duplicated(project_names)], collapse = ", "))
  }

  # Create the MultiDromaSet object
  if (merge_metadata) {
    multi_set <- MultiDromaSet(
      name = project_names,
      DromaSets = droma_sets
    )
  } else {
    # Don't merge metadata, create empty metadata
    multi_set <- MultiDromaSet(
      name = project_names,
      DromaSets = droma_sets,
      sampleMetadata = data.frame(),
      treatmentMetadata = data.frame()
    )
  }

  message("Created MultiDromaSet with ", length(droma_sets), " projects: ",
          paste(project_names, collapse = ", "))

  return(multi_set)
}

#' Add DromaSet to Existing MultiDromaSet
#'
#' @description Adds a new DromaSet object to an existing MultiDromaSet
#' @param multi_set A MultiDromaSet object
#' @param droma_set A DromaSet object to add
#' @param project_name Optional character, name for the new project (will use DromaSet name if not provided)
#' @param update_metadata Logical, whether to update merged metadata (default: TRUE)
#' @return Updated MultiDromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Create initial MultiDromaSet
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "~/droma.sqlite")
#'
#' # Create a new DromaSet
#' new_set <- createDromaSetFromDatabase("GDSC", "~/droma.sqlite")
#'
#' # Add to MultiDromaSet
#' multi_set <- addDromaSetToMulti(multi_set, new_set)
#' }
addDromaSetToMulti <- function(multi_set, droma_set, project_name = NULL, update_metadata = TRUE) {

  if (!inherits(multi_set, "MultiDromaSet")) {
    stop("multi_set must be a MultiDromaSet object")
  }

  if (!inherits(droma_set, "DromaSet")) {
    stop("droma_set must be a DromaSet object")
  }

  # Get project name
  if (is.null(project_name)) {
    project_name <- droma_set@name
  }

  # Check if project already exists
  if (project_name %in% multi_set@name) {
    stop("Project '", project_name, "' already exists in MultiDromaSet. Use a different project_name or remove the existing project first.")
  }

  # Add the new DromaSet
  multi_set@DromaSets[[project_name]] <- droma_set
  multi_set@name <- c(multi_set@name, project_name)

  # Update metadata if requested
  if (update_metadata) {
    multi_set@sampleMetadata <- mergeSampleMetadata(multi_set@DromaSets)
    multi_set@treatmentMetadata <- mergeTreatmentMetadata(multi_set@DromaSets)

    # Update dataset types
    new_dataset_types <- unique(sapply(multi_set@DromaSets, function(x) x@datasetType))
    new_dataset_types <- new_dataset_types[!is.na(new_dataset_types)]
    multi_set@datasetType <- new_dataset_types
  }

  message("Added project '", project_name, "' to MultiDromaSet")

  return(multi_set)
}

#' Remove DromaSet from MultiDromaSet
#'
#' @description Removes a DromaSet object from an existing MultiDromaSet
#' @param multi_set A MultiDromaSet object
#' @param project_name Character, name of the project to remove
#' @param update_metadata Logical, whether to update merged metadata (default: TRUE)
#' @return Updated MultiDromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Remove a project from MultiDromaSet
#' multi_set <- removeDromaSetFromMulti(multi_set, "CCLE")
#' }
removeDromaSetFromMulti <- function(multi_set, project_name, update_metadata = TRUE) {

  if (!inherits(multi_set, "MultiDromaSet")) {
    stop("multi_set must be a MultiDromaSet object")
  }

  if (!project_name %in% multi_set@name) {
    stop("Project '", project_name, "' not found in MultiDromaSet. Available projects: ",
         paste(multi_set@name, collapse = ", "))
  }

  if (length(multi_set@name) == 1) {
    stop("Cannot remove the last project from MultiDromaSet")
  }

  # Remove the DromaSet
  multi_set@DromaSets[[project_name]] <- NULL
  multi_set@name <- multi_set@name[multi_set@name != project_name]

  # Update metadata if requested
  if (update_metadata) {
    multi_set@sampleMetadata <- mergeSampleMetadata(multi_set@DromaSets)
    multi_set@treatmentMetadata <- mergeTreatmentMetadata(multi_set@DromaSets)

    # Update dataset types
    new_dataset_types <- unique(sapply(multi_set@DromaSets, function(x) x@datasetType))
    new_dataset_types <- new_dataset_types[!is.na(new_dataset_types)]
    multi_set@datasetType <- new_dataset_types
  }

  message("Removed project '", project_name, "' from MultiDromaSet")

  return(multi_set)
}

#' Create MultiDromaSet from All Available Projects in Database
#'
#' @description Creates a MultiDromaSet object containing all available projects found in a database
#' @param db_path Path to the SQLite database
#' @param exclude_projects Character vector of project names to exclude (optional)
#' @param include_projects Character vector of project names to include (if specified, only these will be included)
#' @param load_metadata Logical, whether to load sample and treatment metadata (default: TRUE)
#' @param auto_load Logical, whether to automatically load treatment response and molecular profiles (default: FALSE)
#' @return A MultiDromaSet object containing all available projects
#' @export
#' @examples
#' \dontrun{
#' # Create MultiDromaSet with all projects in database
#' multi_set <- createMultiDromaSetFromAllProjects("~/droma.sqlite")
#'
#' # Create MultiDromaSet excluding specific projects
#' multi_set <- createMultiDromaSetFromAllProjects("~/droma.sqlite",
#'                                                exclude_projects = c("test_project"))
#'
#' # Create MultiDromaSet with only specific projects
#' multi_set <- createMultiDromaSetFromAllProjects("~/droma.sqlite",
#'                                                include_projects = c("gCSI", "CCLE"))
#' }
createMultiDromaSetFromAllProjects <- function(db_path = file.path(path.expand("~"), "droma.sqlite"),
                                             exclude_projects = NULL,
                                             include_projects = NULL,
                                             load_metadata = TRUE,
                                             auto_load = FALSE) {

  if (!file.exists(db_path)) {
    stop("Database file not found: ", db_path)
  }

  # Get all available projects
  con <- connectDROMADatabase(db_path)
  on.exit(closeDROMADatabase(con), add = TRUE)

  all_projects <- listDROMAProjects(con, show_names_only = TRUE)

  if (length(all_projects) == 0) {
    stop("No projects found in database")
  }

  # Filter projects based on include/exclude criteria
  if (!is.null(include_projects)) {
    # Only include specified projects
    missing_projects <- setdiff(include_projects, all_projects)
    if (length(missing_projects) > 0) {
      warning("The following projects were not found in database: ",
              paste(missing_projects, collapse = ", "))
    }
    selected_projects <- intersect(include_projects, all_projects)
  } else {
    # Include all projects except excluded ones
    if (!is.null(exclude_projects)) {
      selected_projects <- setdiff(all_projects, exclude_projects)
    } else {
      selected_projects <- all_projects
    }
  }

  if (length(selected_projects) == 0) {
    stop("No projects selected after applying include/exclude criteria")
  }

  message("Creating MultiDromaSet with projects: ", paste(selected_projects, collapse = ", "))

  # Create MultiDromaSet
  multi_set <- createMultiDromaSetFromDatabase(
    project_names = selected_projects,
    db_path = db_path,
    load_metadata = load_metadata,
    auto_load = auto_load
  )

  return(multi_set)
}
