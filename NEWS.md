# cgalMeshes 1.1.0

- It is now possible to assign colors and scalars to the faces and the vertices of a mesh. These properties are preserved by certain operations. For example, the computation of the connected components preserves all properties, the union preserves the face properties. The normals are also assigned to the mesh as a property, and they are preserved when computing the connected components.

- Writing a mesh to a file includes the colors and the normals in this file if they are present, and reading a mesh file including colors or normals includes these properties in the mesh. 

- New method `clipToPlane` for `cgalMesh` objects, to clip the mesh to a plane.

- New method `getBorders` for `cgalMesh` objects, to extract the boundary cycles from the mesh.

- New method `fillBoundaryHole` for `cgalMesh` objects, to fill the holes of the mesh.


# cgalMeshes 1.0.0

First release.
