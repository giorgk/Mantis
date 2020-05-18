# Overview
MantisServer is a computer program for rapid calculation of concentrations via convoluting loading functions with precomputed Unit Response Functions. 



# Download and compiling
MantisServer depends on [boost](https://www.boost.org/) and [CGAL](https://www.cgal.org/).

Most likely the MantisServer can be built on top of the standard version of the libraries that come any linux distribution. For example one can obtain the libraries via sudo apt-get.

First clone the repositories
```
git clone https://github.com/giorgk/Mantis.git
cd CPP/MantisServer/MantisServer/
```

Once these two dependencies boost and cgal are met MantisServer can be configured using cmake as follows:
```
cmake .
```
If the libraries cannot be found under the standard path then one should define the paths as:
```
cmake -DCGAL_DIR:PATH=/path/to/cgal -DBOOST_DIR=/path/to/boost .
```

Sometimes it is difficult to set the paths following the default installation.

## Using Spack
An alternative solution to this is via [spack](https://spack.readthedocs.io/en/latest/)

The installation of spack involves two steps ([see more](https://spack.readthedocs.io/en/latest/getting_started.html)):
```
git clone https://github.com/spack/spack.git
```
and then run `$SPACK_ROOT/share/spack/setup-env.sh` which allow bash to be aware of the spack commands. This is better to run in the .bashrc startup file

Next we create an environment and make it active
```
spack env create mantis
spack env activate mantis
```

Install the nessecary libraries:
```
spack install boost@1.69.0 target=x86_64
spack install cgal@4.13 ^boost@1.69.0 target=x86_64
```
Sometimes cgal failed to get build because the mprf failed to apply a patch. One way to bypass this is by using an older version of mpfr
```
spack install cgal@4.13 ^boost@1.69.0 ^mpfr@4.0.0
```

After the installation of the packages, as long as the mantis enviroment is the active one, you can verify using
```
spack env list
```
The active if any will be highlighted.

then from the folder `CPP/MantisServer/MantisServer/` type
```
cmake .
```

## Using Vcpkg
Coming soon...



