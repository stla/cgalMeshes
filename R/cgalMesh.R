cgalMesh <- R6Class(
  "cgalMesh",
  
  lock_class = TRUE,
  
  cloneable = FALSE,
  
  private = list(
    .meshXPtr = NULL
  ),
  
  public = list(
    
    initialize = function(mesh, vertices, faces){
      if(!missing(mesh)) {
        if(inherits(mesh, "mesh3d")) {
          mesh <- getVF(mesh)
        }
        VF <- checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
      } else {
        VF <- checkMesh(vertices, faces, aslist = TRUE)
      }
      private[[".meshXPtr"]] <- CGALmesh$new(VF[["vertices"]], VF[["faces"]])
      invisible(self)
    }
    
  )
  
)