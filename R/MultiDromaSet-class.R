#' MultiDromaSet Class
#'
#' @description A class to represent multiple DROMA projects with drug response and omics data.
#' Each MultiDromaSet contains multiple DromaSet objects, allowing for cross-project analyses
#' and handling of overlapping samples between different datasets.
#'
#' @slot name Character vector, the names of the projects included in this MultiDromaSet
#' @slot DromaSets List containing DromaSet objects for each project
#' @slot sampleMetadata Data frame with merged sample annotations from all projects
#' @slot treatmentMetadata Data frame with merged drug annotations from all projects
#' @slot datasetType Character vector, the types of datasets included (e.g., "CellLine", "PDX", "PDO")
#' @slot db_info List containing database connection information
#' @export
setClass("MultiDromaSet",
  slots = c(
    name = "character",
    DromaSets = "list",
    sampleMetadata = "data.frame",
    treatmentMetadata = "data.frame",
    datasetType = "character",
    db_info = "list"
  ),
  prototype = list(
    name = character(),
    DromaSets = list(),
    sampleMetadata = data.frame(),
    treatmentMetadata = data.frame(),
    datasetType = character(),
    db_info = list()
  )
)

#' Create a MultiDromaSet Object
#'
#' @description Creates a MultiDromaSet object to store multiple DromaSet objects for cross-project analysis
#' @param name Character vector, the names of the projects included
#' @param DromaSets List containing DromaSet objects for each project
#' @param sampleMetadata Data frame with merged sample annotations (optional, will be merged from DromaSets if not provided)
#' @param treatmentMetadata Data frame with merged drug annotations (optional, will be merged from DromaSets if not provided)
#' @param datasetType Character vector, the types of datasets included (optional, will be inferred from DromaSets if not provided)
#' @param db_info List containing database connection information
#' @return A MultiDromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Create individual DromaSet objects
#' gCSI <- createDromaSetFromDatabase("gCSI", "~/droma.sqlite")
#' CCLE <- createDromaSetFromDatabase("CCLE", "~/droma.sqlite")
#'
#' # Create a MultiDromaSet object
#' multi_set <- MultiDromaSet(
#'   name = c("gCSI", "CCLE"),
#'   DromaSets = list(gCSI = gCSI, CCLE = CCLE)
#' )
#' }
MultiDromaSet <- function(name,
                         DromaSets = list(),
                         sampleMetadata = data.frame(),
                         treatmentMetadata = data.frame(),
                         datasetType = character(),
                         db_info = list()) {

  # Validate inputs
  if (length(DromaSets) > 0) {
    # Check that all elements are DromaSet objects
    if (!all(sapply(DromaSets, function(x) inherits(x, "DromaSet")))) {
      stop("All elements in DromaSets must be DromaSet objects")
    }

    # If name not provided, extract from DromaSets
    if (missing(name) || length(name) == 0) {
      name <- sapply(DromaSets, function(x) x@name)
    }

    # Ensure names match
    if (length(name) != length(DromaSets)) {
      stop("Length of name must match length of DromaSets")
    }

    # Set names for DromaSets list
    if (is.null(names(DromaSets))) {
      names(DromaSets) <- name
    }

    # Merge metadata if not provided
    if (nrow(sampleMetadata) == 0) {
      sampleMetadata <- mergeSampleMetadata(DromaSets)
    }

    if (nrow(treatmentMetadata) == 0) {
      treatmentMetadata <- mergeTreatmentMetadata(DromaSets)
    }

    # Infer datasetType if not provided
    if (length(datasetType) == 0) {
      datasetType <- unique(sapply(DromaSets, function(x) x@datasetType))
      datasetType <- datasetType[!is.na(datasetType)]
    }

    # Merge db_info if not provided
    if (length(db_info) == 0) {
      db_paths <- unique(sapply(DromaSets, function(x) x@db_info$db_path))
      db_paths <- db_paths[!is.null(db_paths)]
      if (length(db_paths) > 0) {
        db_info <- list(db_path = db_paths[1])  # Use first database path
      }
    }
  }

  # Create new MultiDromaSet object
  object <- new("MultiDromaSet",
              name = name,
              DromaSets = DromaSets,
              sampleMetadata = sampleMetadata,
              treatmentMetadata = treatmentMetadata,
              datasetType = datasetType,
              db_info = db_info)

  # Return the object
  return(object)
}

#' Helper function to merge sample metadata from multiple DromaSets
#' @param DromaSets List of DromaSet objects
#' @return Merged data frame of sample metadata
mergeSampleMetadata <- function(DromaSets) {
  if (length(DromaSets) == 0) {
    return(data.frame())
  }

  # Extract sample metadata from each DromaSet
  sample_list <- lapply(names(DromaSets), function(proj_name) {
    ds <- DromaSets[[proj_name]]
    if (nrow(ds@sampleMetadata) > 0) {
      # Add ProjectID column to track source
      metadata <- ds@sampleMetadata
      metadata$ProjectID <- proj_name
      return(metadata)
    } else {
      return(NULL)
    }
  })

  # Remove NULL entries
  sample_list <- sample_list[!sapply(sample_list, is.null)]

  if (length(sample_list) == 0) {
    return(data.frame())
  }

  # Merge all sample metadata (keep all entries, don't remove duplicates)
  merged_samples <- do.call(rbind, sample_list)

  return(merged_samples)
}

#' Helper function to merge treatment metadata from multiple DromaSets
#' @param DromaSets List of DromaSet objects
#' @return Merged data frame of treatment metadata
mergeTreatmentMetadata <- function(DromaSets) {
  if (length(DromaSets) == 0) {
    return(data.frame())
  }

  # Extract treatment metadata from each DromaSet
  treatment_list <- lapply(names(DromaSets), function(proj_name) {
    ds <- DromaSets[[proj_name]]
    if (nrow(ds@treatmentMetadata) > 0) {
      # Add ProjectID column to track source
      metadata <- ds@treatmentMetadata
      metadata$ProjectID <- proj_name
      return(metadata)
    } else {
      return(NULL)
    }
  })

  # Remove NULL entries
  treatment_list <- treatment_list[!sapply(treatment_list, is.null)]

  if (length(treatment_list) == 0) {
    return(data.frame())
  }

  # Merge all treatment metadata (keep all entries, don't remove duplicates)
  merged_treatments <- do.call(rbind, treatment_list)

  return(merged_treatments)
}

#' Show Method for MultiDromaSet objects
#'
#' @description Displays information about a MultiDromaSet object
#' @param object A MultiDromaSet object
#' @return NULL, prints information to console
#' @export
setMethod("show", "MultiDromaSet", function(object) {
  cat("MultiDromaSet Object\n")
  cat("Projects:", paste(object@name, collapse = ", "), "\n")
  cat("Dataset Types:", paste(unique(object@datasetType), collapse = ", "), "\n\n")

  # Show information for each DromaSet
  cat("Individual DromaSets:\n")
  for (i in seq_along(object@DromaSets)) {
    ds <- object@DromaSets[[i]]
    cat("  ", names(object@DromaSets)[i], ":\n")

    # Treatment Response
    if (length(ds@treatmentResponse) > 0) {
      tr_info <- lapply(ds@treatmentResponse, function(x) {
        if(is.matrix(x) || is.data.frame(x)) {
          paste0("(", nrow(x), " drugs x ", ncol(x), " samples)")
        } else {
          "(data in database)"
        }
      })
      cat("    Treatment Response:", paste(names(ds@treatmentResponse), tr_info, sep = ": ", collapse = ", "), "\n")
    } else {
      cat("    Treatment Response: None loaded\n")
    }

    # Molecular Profiles
    if (length(ds@molecularProfiles) > 0) {
      mp_info <- lapply(ds@molecularProfiles, function(x) {
        if(is.matrix(x) || is.data.frame(x)) {
          paste0("(", nrow(x), " features x ", ncol(x), " samples)")
        } else {
          "(data in database)"
        }
      })
      cat("    Molecular Profiles:", paste(names(ds@molecularProfiles), mp_info, sep = ": ", collapse = ", "), "\n")
    } else {
      cat("    Molecular Profiles: None loaded\n")
    }
  }

  # Sample overlap information
  cat("\nSample Information:\n")
  if (nrow(object@sampleMetadata) > 0) {
    # Get unique samples across all projects
    unique_samples <- length(unique(object@sampleMetadata$SampleID))
    total_samples <- nrow(object@sampleMetadata)
    cat("  Total unique samples:", unique_samples, "\n")

    # Show sample overlap between projects
    if (length(object@DromaSets) > 1) {
      sample_lists <- lapply(object@DromaSets, function(ds) {
        if (nrow(ds@sampleMetadata) > 0) {
          return(ds@sampleMetadata$SampleID)
        } else {
          return(character(0))
        }
      })

      cat("  Sample counts per project:\n")
      for (i in seq_along(sample_lists)) {
        cat("    ", names(sample_lists)[i], ":", length(sample_lists[[i]]), "samples\n")
      }

      # Show pairwise overlaps
      if (length(sample_lists) == 2) {
        overlap_count <- length(intersect(sample_lists[[1]], sample_lists[[2]]))
        cat("  Overlapping samples between projects:", overlap_count, "\n")
      }
    }
  } else {
    cat("  No sample metadata available\n")
  }

  # Treatment overlap information
  cat("\nTreatment Information:\n")
  if (nrow(object@treatmentMetadata) > 0) {
    # Get unique drugs across all projects
    unique_drugs <- length(unique(object@treatmentMetadata$DrugName))
    total_drugs <- nrow(object@treatmentMetadata)
    cat("  Total unique drugs:", unique_drugs, "\n")

    # Show drug overlap between projects
    if (length(object@DromaSets) > 1) {
      drug_lists <- lapply(object@DromaSets, function(ds) {
        if (nrow(ds@treatmentMetadata) > 0) {
          return(ds@treatmentMetadata$DrugName)
        } else {
          return(character(0))
        }
      })

      cat("  Drug counts per project:\n")
      for (i in seq_along(drug_lists)) {
        cat("    ", names(drug_lists)[i], ":", length(drug_lists[[i]]), "drugs\n")
      }

      # Show pairwise overlaps for drugs
      if (length(drug_lists) == 2) {
        overlap_count <- length(intersect(drug_lists[[1]], drug_lists[[2]]))
        cat("  Overlapping drugs between projects:", overlap_count, "\n")
      }
    }
  } else {
    cat("  No treatment metadata available\n")
  }

  # Database information
  cat("\nDatabase Connection Information:\n")
  if (length(object@db_info) > 0) {
    cat("  Path:", ifelse(is.null(object@db_info$db_path), "Not specified", object@db_info$db_path), "\n")
  } else {
    cat("  No database information available\n")
  }

  invisible(NULL)
})

#' Get Available Projects in MultiDromaSet
#'
#' @description Returns the names of projects available in a MultiDromaSet
#' @param object A MultiDromaSet object
#' @return Character vector of project names
#' @export
setGeneric("availableProjects", function(object) standardGeneric("availableProjects"))

#' @rdname availableProjects
#' @export
setMethod("availableProjects", "MultiDromaSet", function(object) {
  return(object@name)
})

#' Get DromaSet by Project Name
#'
#' @description Retrieves a specific DromaSet from a MultiDromaSet by project name
#' @param object A MultiDromaSet object
#' @param project_name Character, the name of the project to retrieve
#' @return A DromaSet object
#' @export
setGeneric("getDromaSet", function(object, project_name) standardGeneric("getDromaSet"))

#' @rdname getDromaSet
#' @export
setMethod("getDromaSet", "MultiDromaSet", function(object, project_name) {
  if (!project_name %in% object@name) {
    stop("Project '", project_name, "' not found. Available projects: ", paste(object@name, collapse = ", "))
  }

  return(object@DromaSets[[project_name]])
})

#' Get Overlapping Samples Between Projects
#'
#' @description Finds samples that are present in multiple projects within a MultiDromaSet
#' @param object A MultiDromaSet object
#' @param projects Character vector, specific projects to check for overlap (default: all projects)
#' @return A list containing overlapping sample information
#' @export
setGeneric("getOverlappingSamples", function(object, projects = NULL) standardGeneric("getOverlappingSamples"))

#' @rdname getOverlappingSamples
#' @export
setMethod("getOverlappingSamples", "MultiDromaSet", function(object, projects = NULL) {
  if (is.null(projects)) {
    projects <- object@name
  }

  # Validate project names
  invalid_projects <- setdiff(projects, object@name)
  if (length(invalid_projects) > 0) {
    stop("Invalid project names: ", paste(invalid_projects, collapse = ", "))
  }

  if (length(projects) < 2) {
    stop("At least 2 projects are required to find overlaps")
  }

  # Get sample lists for each project
  sample_lists <- lapply(projects, function(proj) {
    ds <- object@DromaSets[[proj]]
    if (nrow(ds@sampleMetadata) > 0) {
      return(ds@sampleMetadata$SampleID)
    } else {
      return(character(0))
    }
  })
  names(sample_lists) <- projects

  # Find overlapping samples
  if (length(projects) == 2) {
    # For two projects, simple intersection
    overlap <- intersect(sample_lists[[1]], sample_lists[[2]])
    result <- list(
      projects = projects,
      overlapping_samples = overlap,
      overlap_count = length(overlap),
      project_sample_counts = sapply(sample_lists, length)
    )
  } else {
    # For multiple projects, find samples present in all projects
    overlap <- Reduce(intersect, sample_lists)
    result <- list(
      projects = projects,
      overlapping_samples = overlap,
      overlap_count = length(overlap),
      project_sample_counts = sapply(sample_lists, length),
      pairwise_overlaps = combn(projects, 2, function(pair) {
        list(
          projects = pair,
          overlap = intersect(sample_lists[[pair[1]]], sample_lists[[pair[2]]]),
          count = length(intersect(sample_lists[[pair[1]]], sample_lists[[pair[2]]]))
        )
      }, simplify = FALSE)
    )
  }

  return(result)
})

#' Load Molecular Profiles Across Multiple Projects
#'
#' @description Loads molecular profiles from multiple projects,
#' returning only data for overlapping samples
#' @param object A MultiDromaSet object
#' @param molecular_type Character, the type to load (e.g., "mRNA", "cnv") or "all" to load all available types
#' @param features Character vector, specific features to load (optional)
#' @param projects Character vector, specific projects to load from (default: all projects)
#' @param overlap_only Logical, whether to return only overlapping samples (default: FALSE)
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")
#' @return A list containing molecular profile matrices from each project
#' @export
setGeneric("loadMultiProjectMolecularProfiles", function(object, molecular_type, features = NULL, projects = NULL, overlap_only = FALSE, data_type = "all", tumor_type = "all")
  standardGeneric("loadMultiProjectMolecularProfiles"))

#' @rdname loadMultiProjectMolecularProfiles
#' @export
setMethod("loadMultiProjectMolecularProfiles", "MultiDromaSet", function(object, molecular_type, features = NULL, projects = NULL, overlap_only = FALSE, data_type = "all", tumor_type = "all") {
  if (is.null(projects)) {
    projects <- object@name
  }

  # Validate project names
  invalid_projects <- setdiff(projects, object@name)
  if (length(invalid_projects) > 0) {
    stop("Invalid project names: ", paste(invalid_projects, collapse = ", "))
  }

  # Handle "all" molecular_type option
  if (molecular_type == "all") {
    # Get all available molecular types across all projects
    all_mol_types <- character()
    for (proj in projects) {
      ds <- object@DromaSets[[proj]]
      proj_types <- availableMolecularProfiles(ds, include_db = TRUE)
      all_mol_types <- unique(c(all_mol_types, proj_types))
    }

    if (length(all_mol_types) == 0) {
      stop("No molecular profile types found across the specified projects")
    }

    # Load each molecular type separately
    all_results <- list()

    for (mol_type in all_mol_types) {
      message("Loading molecular profile type: ", mol_type)

      tryCatch({
        mol_results <- loadMultiProjectMolecularProfiles(
          object = object,
          molecular_type = mol_type,
          features = features,
          projects = projects,
          overlap_only = overlap_only,
          data_type = data_type,
          tumor_type = tumor_type
        )

        # Store results with molecular type as top-level key
        all_results[[mol_type]] <- mol_results

      }, error = function(e) {
        warning("Failed to load molecular profile '", mol_type, "': ", e$message)
      })
    }

    message("Loaded ", length(all_results), " molecular profile types across ",
            length(projects), " projects")
    return(all_results)
  }

  # Original single molecular_type loading logic
  # Load data from each project
  data_list <- list()

  for (proj in projects) {
    ds <- object@DromaSets[[proj]]

    tryCatch({
      data_list[[proj]] <- loadMolecularProfiles(ds, molecular_type = molecular_type,
                                                features = features, return_data = TRUE,
                                                data_type = data_type, tumor_type = tumor_type)
    }, error = function(e) {
      warning("Problem with loading molecular profiles from project '", proj, "': ", e$message)
      data_list[[proj]] <- NULL
    })
  }

  # Remove NULL entries
  data_list <- data_list[!sapply(data_list, is.null)]

  if (length(data_list) == 0) {
    stop("No molecular profile data could be loaded from any project")
  }

  # If overlap_only is TRUE, filter to overlapping samples
  if (overlap_only && length(data_list) > 1) {
    # Get sample names from each dataset
    sample_lists <- lapply(data_list, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        return(colnames(x))
      } else {
        return(character(0))
      }
    })

    # Find overlapping samples
    overlapping_samples <- Reduce(intersect, sample_lists)

    if (length(overlapping_samples) == 0) {
      warning("No overlapping samples found between projects. Return all data.")
      return(data_list)
    }

    # Filter each dataset to overlapping samples
    data_list <- lapply(data_list, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        common_samples <- intersect(colnames(x), overlapping_samples)
        if (length(common_samples) > 0) {
          return(x[, common_samples, drop = FALSE])
        }
      }
      return(x)
    })

    message("Filtered molecular profile data to ", length(overlapping_samples), " overlapping samples")
  }

  return(data_list)
})

#' Load Treatment Response Across Multiple Projects
#'
#' @description Loads treatment response data from multiple projects,
#' returning only data for overlapping samples
#' @param object A MultiDromaSet object
#' @param drugs Character vector, specific drugs to load (optional)
#' @param projects Character vector, specific projects to load from (default: all projects)
#' @param overlap_only Logical, whether to return only overlapping samples (default: FALSE)
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")
#' @return A list containing treatment response matrices from each project
#' @export
setGeneric("loadMultiProjectTreatmentResponse", function(object, drugs = NULL, projects = NULL, overlap_only = FALSE, data_type = "all", tumor_type = "all")
  standardGeneric("loadMultiProjectTreatmentResponse"))

#' @rdname loadMultiProjectTreatmentResponse
#' @export
setMethod("loadMultiProjectTreatmentResponse", "MultiDromaSet", function(object, drugs = NULL, projects = NULL, overlap_only = FALSE, data_type = "all", tumor_type = "all") {
  if (is.null(projects)) {
    projects <- object@name
  }

  # Validate project names
  invalid_projects <- setdiff(projects, object@name)
  if (length(invalid_projects) > 0) {
    stop("Invalid project names: ", paste(invalid_projects, collapse = ", "))
  }

  # Load data from each project
  data_list <- list()

  for (proj in projects) {
    ds <- object@DromaSets[[proj]]

    tryCatch({
      data_list[[proj]] <- loadTreatmentResponse(ds, drugs = drugs, return_data = TRUE,
                                                data_type = data_type, tumor_type = tumor_type)
    }, error = function(e) {
      warning("Problem with loading treatment response from project '", proj, "': ", e$message)
      data_list[[proj]] <- NULL
    })
  }

  # Remove NULL entries
  data_list <- data_list[!sapply(data_list, is.null)]

  if (length(data_list) == 0) {
    stop("No treatment response data could be loaded from any project")
  }

  # If overlap_only is TRUE, filter to overlapping samples
  if (overlap_only && length(data_list) > 1) {
    # Get sample names from each dataset
    sample_lists <- lapply(data_list, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        return(colnames(x))
      } else {
        return(character(0))
      }
    })

    # Find overlapping samples
    overlapping_samples <- Reduce(intersect, sample_lists)

    if (length(overlapping_samples) == 0) {
      warning("No overlapping samples found between projects. Return all data.")
      return(data_list)
    }

    # Filter each dataset to overlapping samples
    data_list <- lapply(data_list, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        common_samples <- intersect(colnames(x), overlapping_samples)
        if (length(common_samples) > 0) {
          return(x[, common_samples, drop = FALSE])
        }
      }
      return(x)
    })

    message("Filtered treatment response data to ", length(overlapping_samples), " overlapping samples")
  }

  return(data_list)
})

#' Subset MultiDromaSet Object
#'
#' @description Creates a subset of a MultiDromaSet object containing only specified projects
#' @param x A MultiDromaSet object
#' @param projects Character vector, project names to include in the subset
#' @param drop Logical, whether to drop projects that don't exist (default: FALSE)
#' @param ... Additional arguments (currently unused)
#' @return A new MultiDromaSet object containing only the specified projects
#' @export
#' @examples
#' \dontrun{
#' # Create a subset with specific projects
#' multi_subset <- subset(multi_set, projects = c("gCSI", "CCLE"))
#'
#' # Create a subset with non-existent projects (drop = FALSE will raise error)
#' multi_subset <- subset(multi_set, projects = c("gCSI", "UNKNOWN"), drop = FALSE)
#'
#' # Create a subset with non-existent projects (drop = TRUE will ignore them)
#' multi_subset <- subset(multi_set, projects = c("gCSI", "UNKNOWN"), drop = TRUE)
#' }
setMethod("subset", "MultiDromaSet", function(x, projects, drop = FALSE, ...) {

  # Validate input
  if (!is.character(projects)) {
    stop("projects must be a character vector")
  }

  if (length(projects) == 0) {
    stop("At least one project must be specified")
  }

  # Check which projects exist
  existing_projects <- projects %in% x@name

  if (!any(existing_projects)) {
    stop("None of the specified projects exist in MultiDromaSet")
  }

  # Handle non-existent projects based on drop parameter
  if (!drop && !all(existing_projects)) {
    missing_projects <- projects[!existing_projects]
    stop("The following projects do not exist: ", paste(missing_projects, collapse = ", "))
  }

  # Get only existing projects
  valid_projects <- projects[existing_projects]

  # Extract DromaSets for valid projects
  subset_droma_sets <- x@DromaSets[valid_projects]

  # Filter sample metadata
  if (nrow(x@sampleMetadata) > 0) {
    subset_sample_metadata <- x@sampleMetadata[x@sampleMetadata$ProjectID %in% valid_projects, , drop = FALSE]
  } else {
    subset_sample_metadata <- data.frame()
  }

  # Filter treatment metadata
  if (nrow(x@treatmentMetadata) > 0) {
    subset_treatment_metadata <- x@treatmentMetadata[x@treatmentMetadata$ProjectID %in% valid_projects, , drop = FALSE]
  } else {
    subset_treatment_metadata <- data.frame()
  }

  # Update dataset types based on remaining projects
  subset_dataset_types <- unique(sapply(subset_droma_sets, function(ds) ds@datasetType))
  subset_dataset_types <- subset_dataset_types[!is.na(subset_dataset_types)]

  # Create new MultiDromaSet object
  subset_object <- new("MultiDromaSet",
                     name = valid_projects,
                     DromaSets = subset_droma_sets,
                     sampleMetadata = subset_sample_metadata,
                     treatmentMetadata = subset_treatment_metadata,
                     datasetType = subset_dataset_types,
                     db_info = x@db_info)

  return(subset_object)
})
