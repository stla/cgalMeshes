# cgalMeshes 2.4.0

- New method `$smoothShape`, to smooth a mesh.

- Now the `torusMesh` function can return a mesh of a torus constructed from a conformal parameterization.

- New function `conformalCyclideMesh`, to get a mesh of a Dupin cyclide constructed from a conformal parameterization.


# cgalMeshes 2.3.0

- The function `parametricMesh` was slow with the argument `clean=TRUE`. This has been improved.

- Implementation of the 3D alpha wrapping.

- Optimal bounding box of a mesh.

- New method `$sampleInMesh` to uniformly sample in the volume bounded by a mesh. The `$sample` method has been renamed to `$sampleOnMesh`.


# cgalMeshes 2.2.0

- Implementation of the Minkowski sum of two meshes.

- The two methods `$getEdges` and `$getFacesInfo` return more information.

- The function `parametricMesh` gains three new arguments:
  * `clean`, to merge duplicated vertices;
  * `fcolor`, a function to assign a color at each vertex from its `(u,v)`-coordinates;
  * `fnormal`, a function to assign a normal at each vertex from its `(u,v)`-coordinates.


# cgalMeshes 2.1.0

- The `$clip` method with `volume=TRUE` didn't correctly preserve the face properties (e.g. the face colors). Now it partially works: the face properties of the clipped mesh are correctly preserved, but not the face properties of the clipping mesh (help wanted).

- It is now possible to write a binary PLY file representing a mesh.

- Hausdorff distance between two meshes.

- Random sampling on a mesh.


# cgalMeshes 2.0.0

- It is now possible to assign colors and scalars to the faces and the vertices of a mesh. These properties are preserved by certain operations. For example, the computation of the connected components preserves all properties, the union preserves the face properties. The normals are also assigned to the mesh as such a property, and they are preserved when computing the connected components.

- Writing a mesh to a file includes the colors and the normals in this file if they are present, and reading a mesh file including colors or normals includes these properties in the mesh. 

- The function `torusMesh` now allows to construct a mesh of the torus with a given minor radius and whose centerline passes through three given points.

- New function `parametricMesh`, to get a mesh of a parametric surface.

- New function `revolutionMesh`, to get a mesh of a surface of revolution.

- New function `sphericalTriangle`, to get a mesh of a spherical triangle.

- New function `gyroTriangle`, to get a mesh of a hyperbolic triangle.

- The function `AFSreconstruction` gains a new argument `jetSmoothing` allowing to smooth the points cloud before the surface reconstruction.

- New function `SSSreconstruction` to perform scale-space surface reconstruction.

- New methods for `cgalMesh` objects:
  * `clipToPlane`, to clip the mesh to a plane;
  * `clipToIsoCuboid`, to clip the mesh to an iso-oriented cuboid;
  * `getBorders`, to extract the boundary cycles from the mesh;
  * `fillBoundaryHole`, to fill the holes of the mesh;
  * `isotropicRemeshing`, to remesh the mesh;
  * subdivision methods: `CatmullClark`, `LoopSubdivision`, `DooSabin`, `Sqrt3Subdivision`;
  * `whereIs`, to check whether a point is inside or outside the mesh, when the mesh is closed.


# cgalMeshes 1.0.0

First release.
