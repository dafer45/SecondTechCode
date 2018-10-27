This code accompanies the blog post [Using TBTK to calculate the density of states (DOS) of a 1D, 2D, and 3D square lattice](http://second-tech.com/wordpress/index.php/2018/10/23/using-tbtk-to-calculate-the-density-of-states-dos-of-a-1d-2d-and-3d-square-lattice/) on [second-tech.com](https://www.second-tech.com/wordpress).

To build and run this project, first download the full Second Tech code package by typing
```bash
git clone https://www.github.com/dafer45/SecondTechCode
```
[TBTK](https://github.com/dafer45/TBTK) and [BLAS](http://www.netlib.org/blas/), [LAPACK](http://www.netlib.org/lapack/), and [OpenCV](https://opencv.org/) also needs to be installed. This code has been compiled and tested with TBTK v1.0.2.

Next enter the folder 2018_10_23 and type
```bash
cmake .
make
./build/Application
```

The resulting output can be found in the figures folder.

<b>Contact:</b> kristofer.bjornson@gmail.com
