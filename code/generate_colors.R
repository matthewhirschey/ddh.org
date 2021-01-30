## A function to generate color palettes
generate_colors <- function(hex) {
  #' @param hex Hex code for main color as string
  #' @return A set of three colors
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .6), .2),
    hex,
    colorspace::darken(hex, .6, space = "HLS")
  )
  return(color_set)
}

## MAIN COLORS -----------------------------------------------------------------
##2EC09C  ## cyan
##BE34EF  ## violet
##E06B12  ## orange
##004AAB  ## blue
##F0CE44  ## yellow

## GENES -----------------------------------------------------------------------

## Color sets  for genes
color_set_gene <- generate_colors("#2EC09C")  ## cyan

## Palette function for genes
pal_gene <- grDevices::colorRampPalette(color_set_gene)
assign("color_pal_gene", pal_gene, envir = .GlobalEnv)


## CELLS -----------------------------------------------------------------------

## Color sets for cells
color_set_cell <- generate_colors("#BE34EF")  ## violet

## Palette function for cells
pal_cell <- grDevices::colorRampPalette(color_set_cell)
assign("color_pal_cell", pal_cell, envir = .GlobalEnv)


## DRUGS -----------------------------------------------------------------------

## Color sets for drugs
color_set_drug <- generate_colors("#E06B12")  ## orange

## Palette function for drugs
pal_drug <- grDevices::colorRampPalette(color_set_drug)
assign("color_pal_drug", pal_drug, envir = .GlobalEnv)
