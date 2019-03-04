# Mantis
This is a playground for the development of the forward implementation 
phase of the NPSAT. At the moment there are two implementations.

1. Matlab 
2. Python

## Matlab implementation
Matlab requires just one line code that loads the GUI
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


## Python implementation
Not alot to say here as this is in primitive stage.
For the time being one need to 
[download](https://drive.google.com/drive/u/2/folders/1OH0R6OH5piws8l9tvqBvuX3fu_K-IawE) 
the _URFdata.mat_ and _data4python.mat_ files and add them inside the *Local* 
folder.

## Matlab Server Implementation
This is going to be a temporary implementation (maybe). This is a script that runs only once and does nothing until it receives a trigger file.

The matlab server version runs with the following command
```
MantisTestServer
```
#### Setup
The program requires alot of data and loading can take up to a minute or more depending the system.

When the program starts, reads the option file **MantisServer.opt** which includes the following list of options in the form of KEYWORD:value

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
| 2 | SubBasins | 1 - 3 |
| 3 | CVHM farms | 1 - 21 |
| 4 | B118 | 1 - 45 |
 
The matlab file **MAPS** contains the definition for the region IDs.

To trigger the computation for the Sacramento Valley the first line should be 
```
2 2
``` 

while for running the 19th and 21st CVHM farm the first line becomes 
```
3 21 1 
```
because the 19th farm is in the 21 record in the CVHM farm gis layer and the 21st farm is the 1st record in the same file.

The second line is the **number of crops** to change the percentage.

Next repeat **number of crops** times the following pattern:

**CropID percentage**

An example of the input file which defines loading reductions for 12 crops for two of the farms is:
```
3 19 21
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

Next repeat **number of wells** times the well breakthrough curve which is consist of **number of time steps** numbers separated by space.

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
