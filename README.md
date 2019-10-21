# Computing Minimum Persistent 2-Cycles

This project implements the algorithms described in the following paper:

[Computing Minimal Persistent Cycles: Polynomial and Hard Cases](https://arxiv.org/pdf/1907.04889.pdf)

by Tamal K. Dey, Tao Hou, and Sayan Mandal which appears on the conference SODA20'.

## Group Information

This project is developed by [Sayan Mandal](http://web.cse.ohio-state.edu/%7Emandal.25/) and [Tao Hou](https://taohou01.github.io) under the **Jyamiti group** at the Ohio State University. For more information about our group, please visit Professor [Tamal K. Dey](http://web.cse.ohio-state.edu/~dey.8/)'s page.

## About the Implementation

The implemented algorithms compute minimum persistent 2-cycles for designated (finite/infinite) intervals from a volume data, where the volume data is given in the '[perseus](http://people.maths.ox.ac.uk/nanda/perseus/index.html)' format. In a nutshell, a .perseus file specifies the scalar value of each voxel in the volume, where each voxel is treated as a *3-dimensional* cell in the cubical complex formed by the volume. For more information on the perseus format, please visit [The Perseus Project](http://people.maths.ox.ac.uk/nanda/perseus/index.html) or [GUDHI](https://gudhi.inria.fr/doc/latest/fileformats.html#FileFormatsPerseus).

Note that a cubical complex from a volume data does not produce infinite intervals. Our implemention of the infinite interval algorithm takes a finite interval as input and computes *the minimum 2-cycle born at the birth index*.

## Building

This project contains two independent executables, i.e., **pers2cyc_fin** and **pers2cyc_inf**, for finite and infinite intervals respectively. While having different building dependencies, both softwares are built by [CMake](https://cmake.org) with the following commands:

```
cd pers2cyc_fin [or] pers2cyc_inf
mkdir build
cd build
cmake ../src
make
```
Both softwares are developed and tested under MacOS. Since the C++ codes and CMake files are somewhat cross-platform, the codes should be able to compile and run on Linux (or possibly Windows) if the dependencies are installed correctly. Most building issues can be resolved by inspecting (or changing) the 'CMakeLists.txt' files. For problems about the softwares, please contact the author of each one listed at the end.

### Building dependencies

**pers2cyc_fin** has the following dependencies:

* [Boost](https://www.boost.org)

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

* [GUDHI](https://gudhi.inria.fr)

**pers2cyc_inf** has the following dependencies:

* [Boost](https://www.boost.org)

* [CGAL](https://www.cgal.org)

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

* [GUDHI](https://gudhi.inria.fr)

where Eigen and GUDHI are header files only and Boost and CGAL need the prebuilt binaries.

## Usage

We provide a sample perseus file 'sphere_50_50_50.perseus' where only one interval will be produced. The minimal 2-cycle produced for the interval is a uniformly shaped 2-shpere.

### pers2cyc_fin

To produce cycles for the top 10 longest finite intervals:
```
./pers2cyc_fin -i ../../sphere_50_50_50.perseus -n 10
```

This outputs a folder named "sphere_50_50_50.perseusloops/" which stores 2-cycles of the top 10 intervals in .off format.

To print help information:

```
./pers2cyc_fin -h
```


### pers2cyc_inf

To produce cycles for the top 10 longest (=finite) intervals:

```
./pers2cyc_inf -i 10 ../../sphere_50_50_50.perseus
```

To produce cycles for the 5th longest (finite) interval and the 9 intervals following it:

```
./pers2cyc_inf -i 5,10 ../../sphere_50_50_50.perseus
```

All cycles (.off files) will be written to the current folder of execution and cycle files will be named by the input .perseus file, the sorting criterion of the intervals, and the index of the interval. By default, the intervals are sorted by the difference of the scalar values on the corresponding paired simplices. To sort the intervals by the difference of the filtration indices of the paired simplices, use the option '--idx'.

To print help information:

```
./pers2cyc_inf -h
```

## Raw to Perseus Script

There is a script '*vol2perseus.py*' written in python3 which can convert a raw volume file to a perseus file. The input can be edited in the script.

## Authors

* **pers2cyc_fin**: Sayan Mandal, sayan.mndl90 at gmail.com
* **pers2cyc_inf**: Tao Hou, taohou01 at gmail.com

## License

THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND. NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. 

This software was developed (and is copyrighted by) the Jyamiti group at The Ohio State University. Please do not redistribute this software. This program is for academic research use only. This software uses the Boost library, Gudhi library, CGAL library and Eigen library which are covered under their own licenses.


