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
This is going to be a temporary implementation. This is a script that runs only once and does nothing until it receives a trigger file.
The matlab server version runs with the following command
```
MantisTestServer
```
#### Setup
First is going to load the require data and this can take up to a minute or more depending the system.

The following is a list of data needed and their location relative to the folder this script is running

1. **Unit response functions** under *Local/CVHM/CVHM_ALLURFS_TA*
2. **N loading maps**. These are 8 images under *Local/Ngw_XXXX.tif*
3. **Land use maps** These are 5 images under *Local/model_input_LUXXXX.tif*
4. **Well locations** under *Local/CVHM/CVHMWells*
5. **Central Valley Maps** under *MAPS*

Note that the last name in the path is the actually file. 

#### Trigger computation
To trigger any computation the program expects to find an input file with the name **MantisServer.inp**.
The format of the file is as follows
In one line provide 
MapID {regionIDs}

| Map ID | Region | Valid Region IDS |
| --- | ----------- | ----------- |
| 1 | Central Valley | not needed |
| 2 | SubBasins | 1 - 3 |
| 3 | CVHM farms | 1 - 21 |
| 4 | B118 | 1 - 45 |
 
The matlab file **MAPS** contains the definition for the region IDs.
To trigger the computation for the Sacramento Valley the first line should be 
1 2, while for running the 19th CVHM far the first line becomes 
2 21 because the 19th farm is in the 21 record in the CVHM farm gis layer.

The second line is the *number of crops* to change the percentage.
Next repeat *number of crops* times the following pattern:
CropID percentage

When the program fins the input files creates an empty **MantisServer.lock** file.
This file will be deleted by program once the computations are finished.
It is not safe to attempt to read the output of the program **MantisServer.out** while the **MantisServer.lock** is present.

The format of the **MantisServer.out** is the following:
The first line is the *number of wells* and the *number of time steps*.
Next repeat *number of wells* times the well breakthrough curve which is consist of *number of time steps* numbers separated by space.

For both files input and outpus the separation character is the space.

To terminate the program one way is to copy a file with name **MantisServer.quit**. The file can be empty
