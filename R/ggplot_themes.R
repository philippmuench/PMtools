#' Custom ggplto2 theme
#' @export
themePM <- function(base.size = 11, axis.family = "", base.family = "")
{
  half_line <- base.size/2
  line.size = 0.2
  ggplot2::theme(
    line = ggplot2::element_line(colour = "black", size = line.size,
                        linetype = 1, lineend = "butt"),
    rect = ggplot2::element_rect(fill = "white", colour = "black",
                        size = 0.5, linetype = 1),
    text = ggplot2::element_text(family = base.family, face = "plain",
                        colour = "black", size = base.size,
                        lineheight = 0.9,  hjust = 0.5,
                        vjust = 0.5, angle = 0,
                        margin = ggplot2::margin(), debug = FALSE),

    axis.line = ggplot2::element_line(colour = "black", size = line.size, linetype = 1),
    axis.text = ggplot2::element_text(family = axis.family, size = ggplot2::rel(0.8), colour = "black"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.8*half_line/2),
                               vjust = 1),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.8*half_line/2),
                               hjust = 1),
    axis.ticks = ggplot2::element_line(colour = "black", size = line.size),
    axis.ticks.length = ggplot2::unit(half_line/2, "pt"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.8 * half_line,
                                                b = 0.8 * half_line/2)),
    axis.title.y = ggplot2::element_text(angle = 90,
                                margin = ggplot2::margin(r = 0.8 * half_line,
                                                l = 0.8 * half_line/2)),

    legend.background = ggplot2::element_rect(colour = NA),
    legend.spacing = ggplot2::unit(0.2, "cm"),
    legend.key = ggplot2::element_rect(fill = "white", colour = "white"),
    legend.key.size = ggplot2::unit(1, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
    legend.text.align = NULL,
    legend.title = ggplot2::element_text(hjust = 0, size = base.size),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,

    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(half_line, "pt"), panel.margin.x = NULL,
    panel.margin.y = NULL, panel.ontop = FALSE,

    strip.background = ggplot2::element_rect(colour = "#f0f0f0",fill = "#f0f0f0"),
    strip.text = ggplot2::element_text(family = axis.family, colour = "black", size = ggplot2::rel(0.8)),
    strip.text.x = ggplot2::element_text(size = base.size, family = axis.family, margin = ggplot2::margin(t = half_line,
                                                b = half_line)),
    strip.text.y = ggplot2::element_text(size = base.size, family = axis.family, angle = -90,
                                margin = ggplot2::margin(l = half_line,
                                                r = half_line)),
    strip.switch.pad.grid = ggplot2::unit(0.1, "cm"),
    strip.switch.pad.wrap = ggplot2::unit(0.1, "cm"),

    plot.background = ggplot2::element_rect(colour = "white"),
    plot.title = ggplot2::element_text(size = ggplot2::rel(1.2),
                              margin = ggplot2::margin(b = half_line * 1.2)),
    plot.margin = ggplot2::margin(half_line, half_line, half_line, half_line),
    complete = TRUE)
}

#' Custom ggplto2 break function (exponent)
#' @export
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

#' Custom pseudolog transformation
pseudolog_trans <- function(base = exp(1), from=0)
{
  trans <- function(x) log(x, base) - from
  inv <- function(x) base^(x + from)
  scales::trans_new("mylog", trans, inv, scales::log_breaks(base = base), domain = c(base^from, Inf))
}

#' Signed Pseudo Logarithm
pseudoLog10 <- function(x) { asinh(x/2)/log(10) }
