suppressPackageStartupMessages({
  library(rmarkdown)
})

files <- grep("Rmd$", list.files("./notebooks/"), value = TRUE)
for (file in files) {
  tryCatch(render(paste0("./notebooks/", file)), error = function(e) print(e))
  rm(list = ls())
  gc()
}
