# cgalMeshes

R6 based utilities for 3D meshes.

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/cgalMeshes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/cgalMeshes/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/cgalMeshes/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/cgalMeshes/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->


### Geodesic distance 

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/trefoilKnot.gif)

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/knot-2-5.gif)

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/dragon.png)


### Clipping

![](https://raw.githubusercontent.com/stla/MeshesTools/main/inst/screenshots/Togliatti.gif)

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/clippedCylinders.gif)


### Fairing

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/HopfTorus.gif)


### Decomposition into convex parts

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/pentagrammicPrism.png)


### Convex hull

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/oloid.gif)


### Subdivision methods

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/Hopf_LoopSubdivision.png)


### Hole filling

![](https://raw.githubusercontent.com/stla/cgalMeshes/main/inst/screenshots/holeFilling.png)


### Boolean operations

#### Intersection

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Intersection.png)

![](https://laustep.github.io/stlahblog/posts/figures/tetrahedraCompoundIntersection.gif)

#### Difference

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Difference.png)

#### Union

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Union.png)

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/tetrahedraCompound.gif)


### Advanced front surface reconstruction

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/Bunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/StanfordDragon.png)

*Dummy head:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/DummyHead.png)

*Skull:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/AFSexamples/Skull.png)


### Other tools

Volume, area, centroid, distance between a point and a mesh, connected 
components.


## More features

There are more features in the **github** branch, to install with:

```r
remotes::install_github("stla/cgalMeshes@github")
```

### Poisson reconstruction

*Spider cage:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/SpiderCage.png)

Here is a series of three images which show the effect of this `spacing` 
parameter (0.05, 0.02, 0.005):

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/SolidMobiusStrip_spacings.png)

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/StanfordBunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/SurfaceReconstruction/master/inst/PoissonExamples/StanfordDragon.png)


### Minkowski addition

*Octahedron + sphere:*

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/OctahedronPlusSphere.gif)

*Tetrahedron + truncated icosahedron:*

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/TetrahedronPlusTruncatedIcosahedron.gif)

*Septuaginta + great stellated dodecahedron:*

![](https://raw.githubusercontent.com/stla/MinkowskiSum/main/inst/screenshots/septuaginta_gsdodecahedron.gif)

*Stanford bunny + sphere:*

![](https://raw.githubusercontent.com/stla/MinkowskiSum/main/inst/screenshots/bunny.png)


### Shape smoothing

*Hopf torus:*

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/HopfTorusSmoothed.gif)

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/StanfordBunnySmoothed.gif)



## Blog posts

- ['CGAL' meets 'R6': the 'cgalMeshes' package](https://laustep.github.io/stlahblog/posts/cgalMeshes.html)

- [Update of 'cgalMeshes'](https://laustep.github.io/stlahblog/posts/cgalMeshes2.html)



## License

This package is provided under the GPL-3 license but it uses the C++ library 
CGAL. If you wish to use CGAL for commercial purposes, you must obtain a 
license from the [GeometryFactory](https://geometryfactory.com).
