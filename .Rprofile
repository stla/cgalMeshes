dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "cgalMeshes.dll", package = "cgalMeshes")
  )
}

myinstall <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    try(dllunload())
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)" 
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
  source(here::here("inst", "essais", "clean-Rd3.R"))
}
