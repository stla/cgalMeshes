myinstall <- function() {
  try(pkgload::unload("cgalMeshes"))
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function() {
  try(pkgload::unload("cgalMeshes"))
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)" 
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
  source(here::here("inst", "essais", "clean-Rd3.R"))
}
