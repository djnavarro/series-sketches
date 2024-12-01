
# set up ------------------------------------------------------------------

name <- "sketches"
version <- 8

# define common helper functions & core tools
source(here::here("source", "common.R"), echo = FALSE)
source(here::here("source", "sketches_base.R"), echo = FALSE)

# make sure we haven't accidentally messed up the versions
# assert_version_consistency(name)

# define the sketches system ----------------------------------------------

bernstein <- function(beta, t = seq(0, 1, .01)) {
  n <- length(beta) - 1
  w <- choose(n, 0:n)
  b <- rep(0, length(t))
  for(v in 0:n) {
    b = b + beta[v + 1] * w[v + 1] * t^v * (1 - t)^(n-v)
  }
  b
}

bezier <- new_class(
  name = "bezier",
  parent = S7_object,
  properties = list(
    x = class_numeric,
    y = class_numeric,
    n = new_property(class = class_numeric, default = 100L),
    curve = new_property(
      class = class_data.frame,
      getter = function(self) {
        t <- seq(0, 1, length.out = self@n)
        data.frame(
          x = bernstein(self@x, t),
          y = bernstein(self@y, t)
        )
      }
    )
  ),
  validator = function(self) {
    if (length(self@x) != length(self@y)) return("x and y must have same length")
    if (length(self@x) < 2) return("at least two control points are required")
    if (length(self@n) != 1) return("n must be length 1")
    if (self@n <= 0) return("n must be a non-negative number")
  })

bezier_ribbon <- new_class(
  name = "bezier_ribbon",
  parent = drawable,
  properties = list(
    x          = class_numeric,
    y          = class_numeric,
    xend       = class_numeric,
    yend       = class_numeric,
    xctr_1     = class_numeric,
    yctr_1     = class_numeric,
    xctr_2     = class_numeric,
    yctr_2     = class_numeric,
    width      = class_numeric,
    smooth     = class_numeric,
    n          = class_integer,
    frequency  = class_numeric,
    octaves    = class_integer,
    seed       = class_integer,
    bezier = new_property(
      class = bezier,
      getter = function(self) {
        bezier(
          x = c(self@x, self@xctr_1, self@xctr_2, self@xend),
          y = c(self@y, self@yctr_1, self@yctr_2, self@yend),
          n = self@n
        )
      }
    ),
    points = new_property(
      class = points,
      getter = function(self) {
        x <- self@bezier@curve$x
        y <- self@bezier@curve$y
        displacement <- ambient::fracture(
          noise = ambient::gen_simplex,
          fractal = ambient::fbm,
          x = x,
          y = y,
          frequency = self@frequency,
          seed = self@seed,
          octaves = self@octaves
        ) |>
          ambient::normalize(to = c(0, 1))
        taper <- sqrt(
          seq(0, 1, length.out = self@n) * seq(1, 0, length.out = self@n)
        )
        width <- displacement * taper * self@width
        dx <- self@xend - self@x
        dy <- self@yend - self@y
        points(
          x = c(x - width * dy, x[self@n:1L] + width[self@n:1L] * dy),
          y = c(y + width * dx, y[self@n:1L] - width[self@n:1L] * dx)
        )
      }
    )
  ),
  constructor = function(x = 0,
                         y = 0,
                         xend = 1,
                         yend = 1,
                         xctr_1 = .5,
                         yctr_1 = .5,
                         xctr_2 = 0,
                         yctr_2 = 0,
                         width = 0.2,
                         smooth = 3L,
                         n = 100L,
                         frequency = 1,
                         octaves = 2L,
                         seed = 1L,
                         ...) {
    new_object(
      drawable(),
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      xctr_1 = xctr_1,
      yctr_1 = yctr_1,
      xctr_2 = xctr_2,
      yctr_2 = yctr_2,
      width = width,
      smooth = smooth,
      n = n,
      frequency = frequency,
      octaves = octaves,
      seed = seed,
      style = style(...)
    )
  },
  validator = function(self) {
    if (length(self@x) != 1) return("x must be length 1")
    if (length(self@y) != 1) return("y must be length 1")
    if (length(self@xend) != 1) return("xend must be length 1")
    if (length(self@yend) != 1) return("yend must be length 1")
    if (length(self@width) != 1) return("width must be length 1")
    if (length(self@n) != 1) return("n must be length 1")
    if (length(self@frequency) != 1) return("frequency must be length 1")
    if (length(self@octaves) != 1) return("octaves must be length 1")
    if (length(self@seed) != 1) return("seed must be length 1")
    if (self@width < 0) return("width must be a non-negative number")
    if (self@frequency < 0) return("frequency must be a non-negative number")
    if (self@n < 1L) return("n must be a positive integer")
    if (self@octaves < 1L) return("octaves must be a positive integer")
  }
)


make_sketch <- function(seed, name, version) {

  # specify the output path and message the user
  output <- output_path(name, version, seed, format = "png")
  message("generating art at ", output)

  # data frame with one row per twist
  set.seed(seed)
  palettes <- list(
    readr::read_csv(
      here::here("source", "palette_02.csv"),
      show_col_types = FALSE
    ),
    readr::read_csv(
      here::here("source", "palette_01.csv"),
      col_select = tidyselect::starts_with("col"),
      show_col_types = FALSE
    )
  ) |>
    dplyr::bind_rows()

  row <- sample.int(nrow(palettes), 1)
  palette <- unlist(palettes[row, ])
  palette <- sample(palette)
  n_ribbons <- 100L
  angle <- runif(1, min = 0, max = 2*pi)
  x_focus_1 <- runif(1, min = -2, max = 2)
  y_focus_1 <- runif(1, min = -2, max = 2)
  x_focus_2 <- runif(1, min = -2, max = 2)
  y_focus_2 <- runif(1, min = -2, max = 2)
  pull_1 <- .8
  pull_2 <- .4

  values <- tibble::tibble(
    x = rnorm(n_ribbons, sd = 2),
    y = rnorm(n_ribbons, sd = 2),
    xend = x + cos(angle) * 1,
    yend = y + sin(angle) * 1,
    xctr_1 = (1 - pull_1) * (x + 2 * xend) / 3 + pull_1 * x_focus_1,
    yctr_1 = (1 - pull_1) * (y + 2 * yend) / 3 + pull_1 * y_focus_1,
    xctr_2 = (1 - pull_2) * (x + xend) / 2 + pull_2 * x_focus_2,
    yctr_2 = (1 - pull_2) * (y + yend) / 2 + pull_2 * y_focus_2,
    width = runif(n_ribbons, min = .01, max = 1),
    smooth = 6L,
    n = 100L,
    fill = sample(palette, n_ribbons, replace = TRUE),
    color = fill
  )

  # list of things to draw
  drawables <- purrr::pmap(values, bezier_ribbon)

  # create a sketch from the ribbons and then draw the sketch
  seed_str <- stringr::str_pad(seed, width = 4, pad = "0")
  png(
    filename = output,
    width = 2000,
    height = 2000,
    units = "px",
    bg = palette[1]
  )
  drawables |>
    sketch() |>
    draw(xlim = c(-2, 2), ylim = c(-2, 2))
  dev.off()

}

for(s in 801:900) make_sketch(s, name, version)

