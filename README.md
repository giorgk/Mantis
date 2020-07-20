# Mantis
This is a playground for the development of the forward implementation 
phase of the NPSAT. At the moment there are three implementations.

1. [Matlab Desktop](https://github.com/giorgk/Mantis#matlab-desktop)
2. [Matlab server](https://github.com/giorgk/Mantis#matlab-server)
3. [C++ Server](https://github.com/giorgk/Mantis#c-server)

## Matlab Desktop
Matlab Desktop is essentially a prototype for the web application.  The following line loads the GUI
```
Mantis
```
The first step is to load the data.
To keep the size of this repo small I'm not uploading to github the 
required data but the 
[links where to download the data](https://github.com/giorgk/Mantis/blob/master/Local/Readme.md).

The second step is to adjust the loading reduction for each crop using
the sliders.

The last step is to press run and wait...

### Dependencies
This code depends on some functions of the 
[mSim toolbox](http://subsurface.gr/software/msim/). Simply 
[download](http://subsurface.gr/software/msim/msim-download/) and follow
few simple steps to install it. 
The `msim_compile` command is quite important as it improves the performance almost 8 times, but the code works even without.


## Matlab Server
This is going to be a temporary implementation. This is a script that runs only once and does nothing until it receives a trigger file.

The matlab server version runs with the following command
```
MantisTestServer
```
#### Setup
The program requires alot of data and loading can take up to a minute or more depending the system.

When the program starts, reads the option file **MantisServer.opt** which includes the following list of options in the form of KEYWORD;value

|KEYWORD| Description|
| --- | ----------- | 
|NGWS| The full path for the folder where the Ngw_####.tif files are|
|LUS|The full path for the folder where the model_input_LU####.tif files are|
|URFS| The full path including the filename of the URFs mat file|
|WELLS| The full path including the filename of the Wells mat file|
|MAPS| The full path including the filename of the MAPS mat file|
|WAIT| The waiting time in seconds before checking for input files|

**NOTE**: Don't forget to include the / at the end for the folders
An example configuration file can be find in the repository


#### Trigger computation
To trigger any computation the program expects to find an input file with the name **MantisServer.inp**.

The format of the file is as follows:

The first line is **MapID regionIDs**

where:

| Map ID | Region | Valid Region IDS |
| --- | ----------- | ----------- |
| 1 | Central Valley | not needed |
| 2 | Basins | 1 - 3 |
| 3 | Counties| 1 - 58 *|
| 4 | B118 | 1 - 45 |
| 5 | CVHM farms | 1 - 21 |

*Only the counties within the central valley would return BTCs
 
The matlab file **MAPS** contains the definition for the region IDs.


To trigger the computation for the Sacramento Valley the first line should be 
```
2 2
``` 

while for running the 19th and 21st CVHM farm the first line becomes 
```
5 21 1 
```
because the 19th farm is in the 21 record in the CVHM farm gis layer and the 21st farm is the 1st record in the same file.

The second line is the **number of crops** to change the percentage.

Next repeat **number of crops** times the following pattern:

**CropID percentage**

An example of the input file which defines loading reductions for 12 crops for two of the farms is:
```
5 19 21
12
301 0.5
302 0.5
303 0.4
400 0.7
401 0.7
402 0.6
605 0.3
1451 0.75
1452 0.253
1460 0.5501
212910 0.102
10003 0.458
```

When the program finds an input file creates an empty **MantisServer.lock** file.
This file will be deleted by program once the computations are finished.
It is not safe to attempt to read the output of the program **MantisServer.out** while the **MantisServer.lock** is present. 

Therefore the safest conditions to read the results is that the output file is present and lock file is not. In addition before sending a new input file its best to have either deleted or moved the output file.

The format of the **MantisServer.out** is the following:

The first line is the **number of wells** and the **number of time steps**.

Next repeat **number of wells** times the well breakthrough curve which consist of **number of time steps** numbers separated by space.

For example
```
855 156
0.685492 0.685492 0.685494 0.685496 0.685499 0.685503 ...
0.000000 0.161214 0.322429 0.483646 0.644866 0.806089 ...
0.076698 0.277653 0.478834 0.680057 0.881323 1.082631 ...
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ...
0.000000 0.000004 0.000017 0.000048 0.000102 0.000187 ...
.
.
.
```


For both files input and output the separation character is just a space.

#### Stop MantisServer
To terminate the program one way is to copy a file with name **MantisServer.quit** 

The file can be empty and it will be deleted before the program terminates.

## C++ Server
This is a c++ server implementation of the Mantis code.
Detailed documentation should be available [here](https://codedocs.xyz/giorgk/Mantis/) once I have figure out how to properly set up my doxygen documentation via the codedocs.

In the mean time here are some usefull info about it.

To obtain a list of options execute the following
```
MantisServer.exe -h
```

To run the server in test mode use the following
```
MantisServer.exe -c config_file -t
```
In test mode it doesn't load the NGW and LU maps and therefore does not execute the simulation. However it reads the incoming message and if there are no errors will return a message similar to the one in the actuall mode populated with random numbers.

An example of configuration file can be found [here](https://github.com/giorgk/Mantis/blob/master/CPP/MantisServer/config.dat)

Last, ommit the *t* flag to run the actual server
```
MantisServer.exe -c config_file
```
The loading of the NGW and LU can take several minutes

 Under [CPP](https://github.com/giorgk/Mantis/tree/master/CPP) there is a test client program that sends an input message and receives the result

 ### Format of input message
* __Nyears__ [int]: The number of years to simulate (For the time being this should not be less than 100 and no more than 500 although we can remove easily those limits if desired).
* __ReductionYear__ [int as YYYY]: The year to start the reduction.
* __WaterContent__ [float]: This is the unsaturated zone mobile water content.
* __Scenario Name__ [string]: This is a name that would correspond to the user selected steady state model/period. For example it can be CVHM_70_03 or C2VSIM_99_09
* __MapID__ [int]: This is the id of the user selected background map (see table above for the background map ids) 
* __Nregions__ [int]: The number of selected regions of the MapID, 
* __Region ids__ [int]: _Nregions_ numbers that correspond to the ids of the selected regions
* __Number of categories for reduction__ [int]: The total number of crops that the loading will be changed 
* __Crop id__ [int]: This is the id of the crop (_We should add a file with the list of the crops ids_)
* __Reduction__ [float]: This is the percentage of how much nitrate we want to keep. 1-> no reduction, 0-> 100% reduction

Repeat the last  two lines _Number of categories for reduction_ times.

* __ENDofMSG\n__ append at the end the keyword ```ENDofMSG``` along with the endline character. 
Use space to separate the different values. Therefore the _Scenario Name_ should not have spaces.

 ### Format of output message
 * __STATUS__ 1 for success and 0 for failure. </br>
 If the status is 0 then a message with some info follows: </br>
 If the status is 1 then:

 * __Number of wells in the selected regions__ [int]
 * __BreakThrough Curve values__ [float]: Repeat this for _Nwells_ x _Nyears_. The first _Nyears_ values correspond to the BTC of the first well, the next _Nyears_ values correspond to the second well and so on so forth.

 ## TestClient
 TestClient is a small utility program which sends message to the server and receives the replies.
 If the messages are valid it prints them into a file so that can be read by other programs for further analysis/validation.

 ### Build client.
 The TestClient program requires only boost library. Under the directory ```CPP/TestClient/TestClient``` there is the CMakeLIsts.txt file.
 
 #### Spack
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
 If there is a system default installation of boost and cmake pick the default instead of the one that comes under the spack environment then do the following
 ```
cmake -DBoost_NO_BOOST_CMAKE=TRUE -DBoost_NO_SYSTEM_PATHS=TRUE ..
 ```
#### Vcpkg
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
    ./TestClient
    ```
    Running the client without inputs will send a default scenario to the server and receive the output. The outputs are also printed in the file testClientResults.dat

2. Using incoming message
    ```
    ./TestClient incomingMsg.dat
    ```
    where the ```incomingMsg.dat``` is a file with the incoming message. The file can be split into multiple lines for readability and the TestClient will read reshape the message and send it to the server as one line. The TestClient will append the ENDofMSG keyword.
3. Stop server
    Simply send
    ```
    ./TestClient quit
    ```



