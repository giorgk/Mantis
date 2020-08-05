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
spack install cmake@3.9.4 target=x86_64
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
Using vcpkg seems a bit easier.
First you build vcpkg as explained [here](https://github.com/microsoft/vcpkg#quick-start-unix)

Because of some recent updates of vcpkg, prior to install boost, make sure
`
autoconf autopoint libtool libtool-bin
`
are installed. If not:

```
sudo apt-get install autoconf autopoint libtool libtool-bin
```


Next install boost
```
./vcpkg install boost
```

Sometimes boost failed to get built because there was an error with python. IN such cases, since we dont need python add the flag ```./vcpkg install boost --keep-going```

and cgal
```
./vcpkg install cgal
```
Prior to installing cgal you may have to install yasm if its not already installed.
```
sudo apt-get install yasm
```

During the build if an error log appears like this:
```
configure: error: No usable m4 in $PATH or /usr/5bin (see config.log for reasons).
```
Do install the m4 and re-run install
```
sudo apt-get install m4
```

### Building MantisServer
If cmake hasn't been installed
```
sudo apt install cmake
```

If both Boost and Cgal are successfully built make a build directory under ```
Mantis/CPP/MantisServer/MantisServer ```
and run cmake

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake ..
```
It's very important to buid the release version.

If there are no errors then compile it as
```
make
```

### Running MantisServer
```
./MantisServer -v
```
will show the current version

```
./MantisServer -h
```
will return a help message and a list of options that one can configure at runtime

```
./MantisServer -c configfile.dat
```
will start the server using the configuration described in configfile.dat. An example configuration file can be found [here](https://github.com/giorgk/Mantis/blob/master/CPP/MantisServer/config.dat)




