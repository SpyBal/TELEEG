altcoord_curvedpolar =function (theta = "x", start = 0, direction = 1, clip = "on", 
          halign = c("center")) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")   "y" else "x"
  ggproto(NULL, CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction), clip = clip, halign = halign, 
          render_fg = function(self, panel_params, theme) {
            if (is.null(panel_params$theta.major)) {
              return(element_render(theme, "panel.border"))
            }
            txt_el <- calc_element("axis.text.x", theme)
            if (inherits(txt_el, "element_blank")) {
              out <- grobTree(zeroGrob(), element_render(theme, 
                                                         "panel.border"))
              return(out)
            }
            theta <- ggplot2:::theta_rescale(self, panel_params$theta.major, 
                                   panel_params)
            labels <- panel_params$theta.labels
            theta <- theta[!is.na(theta)]
            ends_apart <- (theta[length(theta)] - theta[1])%%(2 * 
                                                                pi)
            if (length(theta) > 0 && ends_apart < 0.05) {
              n <- length(labels)
              if (is.expression(labels)) {
                combined <- substitute(paste(a, " ", b), list(a = labels[[1]], 
                                                              b = labels[[n]]))
              }
              else {
                combined <- paste(labels[1], labels[n], sep = " ")
              }
              labels[[n]] <- combined
              labels <- labels[-1]
              theta <- theta[-1]
            }
            element_gp <- gpar(fontsize = rep(txt_el$size, length(labels)), 
                               col = rep(txt_el$colour, length(labels)), fontfamily = txt_el$family, 
                               fontface = txt_el$face, lineheight = txt_el$lineheight)
            wid <- mean(diff(theta))
            path_t <- seq(-wid/2, wid/2, len = 1000)
            id <- rep(seq_along(labels), each = length(path_t))
            theta <- as.vector(t(outer(theta, path_t, "+")))
            x <- 0.45 * sin(theta) + 0.5
            y <- 0.45 * cos(theta) + 0.5
            grobTree(if (length(labels) > 0) 
              textpathGrob(labels, x = x, y = y, id = id, 
                           hjust = txt_el$hjust, vjust = txt_el$vjust, 
                           halign = halign, gp_text = element_gp, gp_path = gpar(linetype = 0, 
                                                                                 lty = 0), upright = TRUE, default.units = "native"), 
              element_render(theme, "panel.border"))
          })
}
