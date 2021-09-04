# Mantis
Mantis is a server application that simulates within seconds multicentury goundwater contaminant transport from diffuse pollution sources. 

A detailed discussion about the concepts of Mantis and our Non Point Source Assessment Toolbox (NPSAT) can be found on the paper

[Kourakos, G., F.Klein, A.Cortis, and T.Harter (2012), A groundwater nonpoint source pollution modeling framework to evaluate long-term dynamics of pollutant exceedance probabilities in wells and other discharge locations, Water Resour. Res., 48, W00L13, doi:10.1029/2011WR010813](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR010813)



# How to use Mantis
To use Mantis one need to prepare a large number of data. The [wiki](https://github.com/giorgk/Mantis/wiki) is the best source to obtain information about that.

Once the initialization data are properly setup one can launch the server with the following command

```
MantisServer.exe -c config_file
```
where `config_file` is the configuration`initialization file of the server. Details about the file can be found [here](https://github.com/giorgk/Mantis/wiki/Initialization-Data).

Other flags the the server supports are:
```
MantisServer.exe -h
```
To obtain a list of options execute the following
and
```
MantisServer.exe -v
```
to get the current version number.

The Server does nothing until it receives an [input message](https://github.com/giorgk/Mantis/wiki/Runtime-Data).

The input message triggers the simulation and the server replies with an [output message](https://github.com/giorgk/Mantis/wiki#output-data).

The easiest way to run the simulations is via the web interface. However it may be usefull if one for example wants to run in batch mode multiple simulations to use a local server. To make it easier to communicate with a local server we provide a client program that can send/receive message to/from server.

## Download Mantis
A windows executable is available under [BIN](https://github.com/giorgk/Mantis/tree/master/CPP/BIN). To run this download the *.exe and all dlls in the same folder and make sure antivirus is not blocking it.

For linux and windows one can build it by using the cmake for configuring with the following options:
```
cmake \
-DCMAKE_TOOLCHAIN_FILE=path\to\vcpkg\scripts\buildsystems\vcpkg.cmake \
-DUSEHF=True
```
The flag `USEHF` enable or disable the HDF5 loading functionality.



 # Mantis Client

 The MantisClient utility program can be found under [BIN](https://github.com/giorgk/Mantis/tree/master/CPP/BIN). It is used to send  messages to the server and receive the replies.
 If the messages are valid it prints them into a file so that can be read by other programs for further analysis/validation.

 ## Build client.
 The TestClient program requires only boost library. Under the directory ```CPP/TestClient/TestClient``` there is the CMakeLIsts.txt file.
 
 ### Spack
 First, if you use spack load the correct environment. (TEST is the name of the environment where I built all the libraries for Mantis. It would make more sence to name this environment MANTIS for example)

 e.g
 ```
spack env activate TEST
 ```
Make sure you are under the ```CPP/TestClient/TestClient``` directory. Then do
 ```
mkdir build
cd build
cmake ..
 ```
 If there is a system default installation of boost but cmake picks the default instead of the one that comes under the spack environment then do the following
 ```
cmake -DBoost_NO_BOOST_CMAKE=TRUE \
-DBoost_NO_SYSTEM_PATHS=TRUE \
-DBOOST_ROOT:PATHNAME=/path/to/spack/var/spack/environments/myenv/.spack-env/view/ ..
 ```


### Vcpkg
The workflow is the same with vcpkg the difference is the cmake configuration
```
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake ..
```
If everything is succesfull then run
```
make
```
to build the program

### Run Client
1. Default
    ```
    ./MantisClient
    ```
    Running the client without inputs will send a default scenario to the server and receive the output. The outputs are also printed in the file testClientResults.dat

2. Using incoming message
    ```
    ./MantisClient incomingMsg_v1.dat outputfile.dat
    ```
    where the ```incomingMsg.dat``` is a file with the incoming message. The file can be split into multiple lines for readability and the MantisClient will read and reshape the message and send it to the server as one line. 

3. Stop server
    Simply send
    ```
    ./MantisClient quit
    ```





