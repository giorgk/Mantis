# Mantis
Is a server application that executes the forward phase of the NPSAT. 

During the development we developed a couple of different versions. However the focus is on the c++ server implementation.

1. [C++ Server](https://github.com/giorgk/Mantis#c-server)
2. [Matlab Desktop](https://github.com/giorgk/Mantis#matlab-desktop) (obsolete)
3. [Matlab server](https://github.com/giorgk/Mantis#matlab-server) (obsolete)

## C++ Server
The C++ version is the one that we actually support. 

Eventually a detailed documentation of the code will be available [here](https://codedocs.xyz/giorgk/Mantis/) once I have figure out how to properly set up my doxygen documentation via the codedocs.

In the mean time here are some usefull info about it.

To obtain a list of options execute the following
```
MantisServer.exe -h
```

To run the server in test mode use 
```
MantisServer.exe -c config_file -t
```
In test mode it doesn't load the N loading data and therefore does not execute the simulation. However it reads the incoming message and if there are no errors will return a message similar to the one in the actuall mode populated with random numbers.

To run the program in simulation mode just ommit the *t* flag.
```
MantisServer.exe -c config_file
```


The list of inputs of the configuration file is shown with the help option `MantisServer.exe -h`. </br> 
An example of configuration file can also be found [here](https://github.com/giorgk/Mantis/blob/master/CPP/MantisServer/mantisConfig.ini)

### Configuration options
The configuration options are divided into sections:

#### CV_Raster
The CV Raster section contains information about hte active area in Central valley and 
the discretization. Contains the following options:

* __Ncells__ : is the total number of active pixels of the Central Valley.
* __Nrows__ : is the number of rows of the raster
* __Ncols__ : is the number of the columns of the raster
* __Raster__ : is the the name of the file that containts the active area. </br>
The raster is a 2D array that contains -1 to the pixels that correspond to the area outside the domain. Every pixel in the active area contains a unique id from 0 to __Ncells-1__ 

#### Data
The Data section provides the input files of the main data
* __MAPS__ : This is a file which containts all the background maps.
* __NO3__ a filename that contains a list of loading scenarios.</br>
Each line in the file describes an N loading scenario as follows </br>
[TYPE] [NAME] [filename] </br>
where </br>[TYPE] is one of the `GNLM` or `SWAT` flags </br>
[NAME] is a descriptive name of the loading scenario that should not have any space or weird symbols. This name it is used in the incoming message in the `loadScen` variable. </br>
[filename] is the path  where the data live. e.g </br>
GNLM GNLM MantisData/GNLM_LU_NGW.dat </br>
SWAT SWAT1 MantisData/SWAT_LOADING_SCEN_1.dat </br>
SWAT SWAT2 MantisData/SWAT_LOADING_SCEN_2.dat </br>
If the line starts with ```#``` then it is ignore. This is a way to comment lines.

* __WELLS__ A file that lists the files (one in each row) that contain the wells for each flow scenario.
* __URFS__ A file that lists the files (one in each row) that contain the URFS for each flow scenario.
* __UNSAT__ A file that contains the values of `Depth/Recharge`. 

* __Path__ If the path is not empty then all the input files will be relative to this path

#### ServerOptions
Options specific to server
* __PORT__ The port number
* __NTHREADS__ The number of threads.
* __Prefix__ This is used as a prefix if ```DebugID``` is provided in the input scenario. It will use this prefix for the various files that are printed

_The most recent input files are in the MantisData folder under the Mantis Google Drive main folder_.

### Input message 
The server does nothing until it receives an input message. </br>
The input message has the KEYWORD VALUE format using spaces as separation character. The input message has to be one line. The end line character `\n`  indicates the end of the message.
#### Input keywords values
The keywords have to use the exact lower Capital case as it appears on the list
* __endSimYear__ [integer YYYY] The simulation always starts at 1945 and continues up to this year. The value has to be greater than 1945 of course. Currently there is an upper hardcoded limit of 2500. Yet if  the end year is less than 1990 or greater than 2500 it gets reset to 2100.
* __startRed__ [integer YYYY] The year to start the reduction applications. This should be always between 1945 and `endSimYear`. If not it gets reset to 2020
* __endRed__ [integer ####] The year to fully implement the reduction application rates. This should always be between `startRed+1` and `endSimYear`. If not it gets reset to `startRed` + 1. 
* __flowScen__ [string] This is a keyword from the following list: 

|Flow scenarios| Description |
|--|---|
|C2VsimRun01Ref6 | Simulation based on C2Vsim average flow conditions. Each basin has beem averaged on a different period |


* __loadScen__ [string] This is a keyword from the following list:

|N Loading scenarios | Description |
|--|---|
|GNLM | The N loading is based on [GNLM](https://ucd-cws.github.io/nitrates/maps/) historic and future predictions. It covers a period between 1945 - 2050 with 15 years increments |
|SWAT1 | Concentrations history (1990 - 2015) based on _Baseline_|
|SWAT2 | Concentrations history (1990 - 2015) based on _High Fertilization_|
|SWAT3 | Concentrations history (1990 - 2015) based on _High Irrigation_|
|SWAT4 | Concentrations history (1990 - 2015) based on _High Irrigation and High Fertilization_|
Eventually the SWAT scenarios will use a mixed of GNLM and SWAT loading. At the moment the SWAT period 1990-2015 is repeated during the simulation.

* __unsatScen__[string] The unsaturated scenario name. Valid options are

| Unsaturated scenarios | Description |
|--|--|
|C2VSIM_SPRING_2000 | Depth and recharge for spring 2000. This corresponds to slow travel times |
|C2VSIM_SPRING_2015 | Depth and recharge for spring 2015. This corresponds to fast travel times |

For any other non valid option the Unsaturated travel is disabled

* __unsatWC__ [float] This is the unsaturated mobile water content coefficient. THis is a multiplier coefficient to the Depth/Recharge value. 

* __bMap__ [string] The name of the background map. This should be one of the following values:

|Background map keys | Description |
|--|---|
|CentralValley | This is the entire Central Valley |
|Basins | The CV is divided into 3 Subbasins|
|Counties | The CV is divided into 58 counties|
|B118 | The CV is divided into 45 groundwater basins|
|Townships | The CV is divided into 703 townships|
|Subregions | The CV is divided into 21 subregions named as _Subregions_|

* __Nregions__ [integer string1 string2,...,stringN] This is the number of subregions to consider during the simulation. This number is followed by `Nregions` names of the regions. Therefore the format should look like the following:</br>
1 CentralValley (If the user has selected _CentralValley_ as background map) </br>
2 SanJoaquinValley TulareLakeBasin (if the user has selected _Basins_ as background map)</br>
4 Subregion2 Subregion17 Subregion12 Subregion8 (if the user has selected _Subregions_ as background map)

#### Codes for Regions
1. _CentralValley_: has only one option which is  `CentralValley`

2. _Basins_: is divided into `SacramentoValley`, `SanJoaquinValley` and `TulareLakeBasin`

3. _Counties_: The list of counties can be found in the shapefile _counties_simple_ under the field _name_. The names in the field name containts spaces, which have to be stripped.

4. _B118_: The list of B118 can be found in the shapefile _B118_simple_ under the field _Basin_Subb_. The names have a format similar to 5-22.13, 2-31 etc. The dashes and dots have to be replaced by `_` For example the above codes will be converted to 5_22_13, 2_31 etc

5. _Townships_: The list of Townships can be found in the shapefile _CVHM_Townships_3310_simplified_ under the field _CO_MTR_.

6. _Subregions_: The C2Vsim subregions take their names by appending to the word `Subregion` the _IRGE_ field of the shapefile _C2Vsim_Subregions_3310_.

* __Ncrops__ [integer, integer1 float1, integer2 float2,...,integerN floatN] The first integer is the number of crops to select for loading reduction. Then `Ncrops` pair of [int float] values which correspond to crop ids reduction percent. </br>
__VERY IMPORTANT NOTE:__ </br>
 >So far the percentage was interpreted as the amount of loading to keep. This was a bit confusing and could not address the option to increase the loading. So in the new version the percentage corresponds to reduction. If the loading of the base case is 30 mg/l and the reduction is 0.6, the final loading will be 30*0.6 e.g keep 40% of the base case.

 #### Codes for crops
 For the GNLM the list of crops is identical to the one used in the LanduseTable_2017_0515 file which can be found [here](https://github.com/thharter/GNLM/tree/master/Input_Data). The column DWR/CAMLCode is what the server expects to find.

 For the SWAT the codes the crop IDs are in the SWAT_LULC.csv under the Local folder of the repository. Valid ID codes are 1-57

 At some point the crops codes will be identical for all loading scenarios.

* __minRch__ [optional] The default value is 0.000027 which corresponds to 10 mm/year. </br>
This has effect on GNLM loading only. During the conversion from kg/ha to mg/l if the recharge is less than this value the concentration is set to zero


* __DebugID__ [optional] If this is present in the input message then the simulation prints 4 files with the following names 
    - DebugID_urf.dat with the all the urfs expanded
    - DebugID_lf_dat with all the loadinf functions
    - DebugID_btc.dat with the streamline breaktrhough curves
    - DebugID_well_btc with the averaged well btcs.
    
    Printing all the files is time consuming while the urf, lf and btc files tend to be very large for areas with many wells. </br>
    __This should be used for debug purposes only__


* __ENDofMSG__ This is a keyword that when found indicates that the message has been read in a correct way.




### Output message 
* __STATUS__ 1 for success and 0 for failure. </br>
 If the status is 0 then a message will follow. </br>
 If the status is 1 then the following info is listed:

 * __Number of wells in the selected regions__ [int]

 * __Number of years__ [int]

 * __BreakThrough Curve values__ [float]: Repeat this for _Nwells_ x _Nyears_. The first _Nyears_ values correspond to the BTC of the first well, the next _Nyears_ values correspond to the second well and so on so forth.

 * __ENDofMSG\n__ is appended at the end of the message along with the endline character. 
 

 ## Client Test program

 A TestClient utility program can be found under [CPP](https://github.com/giorgk/Mantis/tree/master/CPP). It is used to send  messages to the server and receive the replies.
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
 If there is a system default installation of boost but cmake picks the default instead of the one that comes under the spack environment then do the following
 ```
cmake -DBoost_NO_BOOST_CMAKE=TRUE \
-DBoost_NO_SYSTEM_PATHS=TRUE \
-DBOOST_ROOT:PATHNAME=/path/to/spack/var/spack/environments/myenv/.spack-env/view/ ..
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
    ./TestClient incomingMsg_v1.dat
    ```
    where the ```incomingMsg.dat``` is a file with the incoming message. The file can be split into multiple lines for readability and the TestClient will read and reshape the message and send it to the server as one line. The TestClient will append the ENDofMSG keyword.
3. Stop server
    Simply send
    ```
    ./TestClient quit
    ```

----
__Anything below this line is obsolete__
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





