# Helper function for dealing with columns
process_column <- function(data = metadata, 
                           column, 
                           new_name = NULL, 
                           type = "character", 
                           factor_levels = NULL, 
                           factor_order = NULL, 
                           replace_na = NULL, 
                           rename_factors = NULL, 
                           boolean_map = NULL, 
                           domain = NULL, 
                           print_table = TRUE, 
                           log_changes = TRUE) {
  
  # Ensure the column exists
  if (!column %in% names(data)) {
    stop(paste("Column", column, "does not exist in the data frame."))
  }
  
  # Store the original column name for tracking
  original_column <- column
  
  # Rename the column if a new name is provided
  if (!is.null(new_name)) {
    names(data)[names(data) == column] <- new_name
    column <- new_name
  }
  
  # Set the type of the column
  if (!is.null(type)) {
    if (type == "factor") {
      
      # Convert to factor and optionally set levels
      if (!is.null(factor_levels)) {
        data[[column]] <- factor(data[[column]], levels = factor_levels)
      } else {
        data[[column]] <- as.factor(data[[column]])
      }
      
      # Rename factor levels if provided
      if (!is.null(rename_factors)) {
        levels(data[[column]]) <- plyr::mapvalues(levels(data[[column]]), 
                                                  from = names(rename_factors), 
                                                  to = rename_factors)
      }
      
      # Reorder factor levels if provided (after renaming)
      if (!is.null(factor_order)) {
        data[[column]] <- factor(data[[column]], levels = factor_order, ordered = TRUE)
      }
      
    } else if (type == "character") {
      data[[column]] <- as.character(data[[column]])
    } else if (type == "numeric") {
      data[[column]] <- as.numeric(data[[column]])
    } else if (type == "double") {
      data[[column]] <- as.double(data[[column]])
    } else if (type == "integer") {
      data[[column]] <- as.integer(data[[column]])
    } else if (type == "boolean") {
      # Handle boolean type with custom mapping
      if (is.null(boolean_map)) {
        stop("For boolean type, a boolean_map must be provided.")
      }
      data[[column]] <- ifelse(
        tolower(data[[column]]) %in% tolower(boolean_map$true), 
        TRUE, 
        ifelse(
          tolower(data[[column]]) %in% tolower(boolean_map$false), 
          FALSE, 
          NA
        )
      )
    } else if (type == "date") {
      # Convert to Date type, default format is "%Y-%m-%d"
      data[[column]] <- as.Date(data[[column]], format = "%Y-%m-%d")
    } else {
      stop("Unsupported type specified.")
    }
  }
  
  # Replace specific strings with NA for character columns
  if (!is.null(replace_na) && is.character(data[[column]])) {
    data[[column]][data[[column]] %in% replace_na] <- NA
  }
  
  # Print a table of the column's values to the console
  if (print_table) {
    cat("\nTable of", column, ":\n")
    print(table(data[[column]], useNA = "ifany"))
  }
  
  # Compile parameters used into a list
  parameters_used <- list(
    original_column = original_column,
    new_name = new_name,
    type = type,
    factor_levels = factor_levels,
    factor_order = factor_order,
    replace_na = replace_na,
    rename_factors = rename_factors,
    boolean_map = boolean_map,
    domain = domain
  )
  
  # Log the changes if enabled
  if (log_changes) {
    global_log_entry <- list(
      changes = parameters_used,
      timestamp = Sys.time()
    )
    # Use the original column name as the name for the log entry
    column_change_log[[original_column]] <<- global_log_entry
  }
  
  # Return the updated data
  return(data)
}