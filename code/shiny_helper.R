# conditionalPanel that only shows field is not 0
# using array bracket notation to be compatible with namespaces
notZeroConditionalPanel <- function(fieldname, ...) {
  condition_str <- paste0("input['", fieldname, "'] != 0")
  conditionalPanel(condition = condition_str, ...)  
}

# toggleConditionPanel <- function(fieldname, ...) {
#   condition_str <- paste0("input['", fieldname, "'] %% 2 != 0")
#   conditionalPanel(condition = condition_str, ...)
# }

# dynamicSwitch <- function(fieldname, ...) {
#   prettySwitch(
#     inputId = fieldname,
#     value = FALSE,
#     label = "Dynamic",
#     fill = TRUE, 
#     status = "primary")
# }

# SPINNERS ----------
withSpinnerColor <- function(ui_element, plot_type, ...) { #plot type = "gene", "cell", "protein", "compound"
  withSpinner(ui_element,
              type = 3, 
              color = ddh_colors(plot_type), #use ddh_colors fun to give back hex code
              color.background = "#FFFFFF")
}
