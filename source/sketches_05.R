
# set up ------------------------------------------------------------------

name <- "sketches"
version <- 5

# define common helper functions & core tools
source(here::here("source", "common.R"), echo = FALSE)
source(here::here("source", "sketches_base.R"), echo = FALSE)

# make sure we haven't accidentally messed up the versions
# assert_version_consistency(name)

# define the sketches system ----------------------------------------------

# helper function
smooth_bridge <- function(n, scale = .1, smooth = 0, seed = 1L) {
  withr::with_seed(
    seed = seed,
    code = {b <- c(0, e1071::rbridge(1, n - 1))}
  )
  b <- b * scale
  if (smooth > 0) {
    for(i in 1:smooth) {
      b <- (b + c(b[-1], 0)/2 + c(0, b[-n])/2)/2
    }
  }
  b
}

rbezier <- function(n, x, y, xend, yend, xctr, yctr) {
  df <- bezier::bezier(
    t = seq(0, 1, length.out = n),
    p = c(xctr, yctr),
    start = c(x, y),
    end = c(xend, yend)
  ) |>
    matrix(ncol = 2) |>
    as.data.frame()
  names(df) <- c("x", "y")
  tibble::as_tibble(df)
}

bezier_ribbon <- new_class(
  name = "bezier_ribbon",
  parent = drawable,
  properties = list(
    x          = class_numeric,
    y          = class_numeric,
    xend       = class_numeric,
    yend       = class_numeric,
    xctr       = class_numeric,
    yctr       = class_numeric,
    width      = class_numeric,
    smooth     = class_numeric,
    n          = class_integer,
    frequency  = class_numeric,
    octaves    = class_integer,
    seed       = class_integer,
    path = new_property(
      class = points,
      getter = function(self) {
        bezier <- rbezier(
          n = self@n,
          x = self@x,
          y = self@y,
          xend = self@xend,
          yend = self@yend,
          xctr = self@xctr,
          yctr = self@yctr
        )
        x_disp <- smooth_bridge(
          n = self@n,
          smooth = self@smooth,
          scale = 0.1 * self@width,
          seed = self@seed
        )
        y_disp <- smooth_bridge(
          n = self@n,
          smooth = self@smooth,
          scale = 0.1 * self@width,
          seed = self@seed + 1
        )
        points(
          x = bezier$x + x_disp,
          y = bezier$y + y_disp
        )
      }
    ),
    points = new_property(
      class = points,
      getter = function(self) {
        x <- self@path@x
        y <- self@path@y
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
                         xctr = .5,
                         yctr = .5,
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
      xctr = xctr,
      yctr = yctr,
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
    if (length(self@xctr) != 1) return("xctr must be length 1")
    if (length(self@yctr) != 1) return("yctr must be length 1")
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
  palettes <- readr::read_csv(
    here::here("source", "palette_02.csv"),
    show_col_types = FALSE
  )
  row <- sample.int(nrow(palettes), 1)
  palette <- unlist(palettes[row, ])
  n_ribbons <- 400L
  values <- tibble::tibble(
    x = rnorm(n_ribbons, sd = 2),
    y = rnorm(n_ribbons, sd = 2),
    xend = x + 1,
    yend = y + rnorm(n_ribbons, sd = 1),
    xctr = rnorm(n_ribbons),
    yctr = rnorm(n_ribbons),
    width = runif(n_ribbons, min = .1, max = .3),
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

for(s in 101:200) make_sketch(s, name, version)

