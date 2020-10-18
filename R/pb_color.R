#' Color Palette for PLSDA-batch
#'
#' The function outputs a vector of colors.
#'
#' @param num.vector An integer vector specifying which color to use in the palette
#' (there are only 25 colors available).
#'
#' @return
#' A vector of colors (25 colors max.)
#'
#' @examples
#' my.colors = pb_color(1:5)
#'
#' @export
pb_color <- function(num.vector){
  hex_codes1 <- scales::hue_pal(l = 65, c = 100)(10)
  hex_codes2 <- scales::hue_pal(l = 40, c = 60)(5)

  colorlist <- c(hex_codes1[c(1,4,7,2,8,3,9,5,10,6)], mixOmics::color.mixo(c(1,2,6,3,9,4,5,7,10)),
                 hex_codes2[1:3], mixOmics::color.mixo(8), hex_codes2[4:5])

  mycolor <- colorlist[num.vector]
  return(mycolor)
}
