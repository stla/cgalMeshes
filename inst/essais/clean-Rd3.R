#' get i-th character of a string
ithCharacter <- function(string, i) {
  substr(string, i, i)
}

#' get closing brace in a string
#' @param index the index of an opening brace
getClosingBrace <- function(string, index) {
  if(ithCharacter(string, index) != "{") {
    stop("Not an opening brace at the given index.")
  }  
  stack <- logical(0L)
  for(i in index:nchar(string)) {
    char <- ithCharacter(string, i)
    if(char == "{") {
      stack <- c(stack, TRUE)
    } else if(char == "}") {
      stack <- stack[-1L]
      if(length(stack) == 0L) {
        return(i)
      }
    }
  }
  stop("No closing brace found.")
}

#' clean Rd file containing R6 examples
#' @param RdFile Rd file to be cleaned
cleanRd <- function(RdFile){
  
  # read Rd file in single string
  rd <- paste0(readLines(RdFile), collapse = "\n")
  
  # split the string at the closing brace of the \examples section 
  examples_openingBrace <- stringr::str_locate(rd, "examples")[1L, "end"] + 1L
  examples_closingBrace <- getClosingBrace(rd, examples_openingBrace)
  Rd_before <- substr(rd, 1L, examples_closingBrace)
  Rd_after  <- substr(rd, examples_closingBrace + 1L, nchar(rd))
  
  # locate all opening braces of \donttest and \dontrun and get closing braces
  from <- to <- integer(0L)
  if(stringr::str_detect(Rd_after, "donttest")){
    # locate all opening braces of the \donttest and get the closing braces
    openingBraces <- 
      1L + stringr::str_locate_all(Rd_after, "donttest")[[1L]][, "end"]
    closingBraces <- vapply(openingBraces, function(index) {
      getClosingBrace(Rd_after, index)
    }, integer(1L))
    from <- c(openingBraces - 9L, closingBraces)
    to   <- c(openingBraces, closingBraces)
  }
  if(stringr::str_detect(Rd_after, "dontrun")){
    # locate all opening braces of the \dontrun and get the closing braces
    openingBraces <- 
      1L + stringr::str_locate_all(Rd_after, "dontrun")[[1L]][, "end"]
    closingBraces <- vapply(openingBraces, function(index) {
      getClosingBrace(Rd_after, index)
    }, integer(1L))
    from <- c(from, openingBraces - 8L, closingBraces)
    to   <- c(to, openingBraces, closingBraces)
  }
  
  # delete '\donttest{', '\dontrun{' and '}'
  o <- order(from)
  Rd_after <- stringi::stri_sub_replace_all(
    Rd_after, from = from[o], to = to[o], replacement = ""
  )
  
  # write the Rd file
  Rd_new <- paste0(Rd_before, Rd_after)
  writeLines(Rd_new, RdFile)
}

cleanRd("man/cgalMesh.Rd")


