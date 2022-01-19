# All the R functions live under the gwtools (https://github.com/UCD-GW-Nitrate/gwtools) package.
# Either install and load the entire package or download only the mantis.R file and source it.
# Here I'll load the package
library(gwtools)

# Setup server and path -------------------------------------
# First launch the server.
# Open a terminal and run MantisServer.exe -c config.ini
# The following is valid only from my PC
# From the MantisData/v1_8 directory run:
#..\..\CPP\MantisServer\MantisServer\cmake-build-release\MantisServer.exe -c .\mantisConfig.ini

# Configure the scenario options ----------------------------
# The following command creates a list with all available options to configure a scenario
scenario <- mantis.ScenarioOptions()

# First setup the path for the Mantis client program
scenario$client <- "../CPP/TestClient/TestClient/cmake-build-release/MantisClient.exe"

# Quit Server -------------------------
res <- mantis.Quit(scenario)
# If you did that you must restart the server

# Configure the input file name. This is the file where all the options will be written.
# This file name will be used as input to the Mantis client program.
scenario$infile <- 'inputMessage.dat'

# We need to setup a file name for the output file where the Mantis client will print the results
scenario$outfile <- 'BTCresults.dat'

# Test the connection ----------------
res <- mantis.Ping(scenario)
# if the result looks like this [1] "0 0"    "0 pong" its all good


# Optionally we can add an one line short description of the scenario.
# This is printed in the input file and it is not used.
scenario$descr <- 'This is just a test'

# The remaining options are described in the wiki page https://github.com/giorgk/Mantis/wiki/Runtime-Data
# The simulation always start the year 1945, while the end of the simulation can be configured with the option
scenario$endSimYear <- 2150

# Selecting the simulation area -----------------------
# Mantis provides a number of options to choose the study area.
# First we should choose a background map which defines how the Central Valley is divided.
# The available options are given in the table https://github.com/giorgk/Mantis/wiki/Runtime-Data#background-maps
# Let's select a county level
scenario$bMap <- 'Counties'

# Then we can choose which counties we want to include in the simulation. Here we will choose only one
scenario$Regions <- c('Madera')

# Select Nitrate Loading -------------------------------
# The table provides the available options https://github.com/giorgk/Mantis/wiki/Runtime-Data#nitrate-loading-scenarios
# Let's select the GNLM
scenario$loadScen <- 'GNLM'

# Using the GNLM and SWAT loading type we can modify the Nitrate loading by selecting crop types and reduction rates
# Crop variable is a nx2 matrix where the first column is the crop id and the second a multiplier to loading
# Here we will run the base loading by setting to all crops (id = -9) a reduction of 1
scenario$Crops <- matrix(c(-9, 1),nrow = 1, ncol = 2)

# Select well type ------------
# There are three options. See this table  https://github.com/giorgk/Mantis/wiki/Runtime-Data#general-options
# We select here the Virtual Irrigarion/public supply wells
scenario$wellType <- 'VI'


# Run the simulation ---------------
# The simulation involves the following steps:
# 1 write the input file with the scenario options
# 2 Run Mantis client
# 3 The Client will send the options to the Server
# 4 The server will run the simulation and send the results to client
# 5 The client will read the results and return them to the R workspace
# All this sequence is executed with one go as follows:
btc.VI <- mantis.Run(scenario)

# Switch to domestic wells and repeat the simulation
scenario$wellType <- 'VD'
btc.VD <- mantis.Run(scenario)

# Calculate percentiles -----------------
# To calculate percentiles we can use the package matrixStats
library(matrixStats)
btc_prc.VI <- colQuantiles(as.matrix(btc.VI), probs = c(0.5, 0.70, 0.80, 0.90))
btc_prc.VD <- colQuantiles(as.matrix(btc.VD), probs = c(0.5, 0.70, 0.80, 0.90))

# Plot ---------------------
# For the plot my preference is the plotly library which makes nice interactive plots
library(plotly)
{ # Plot percentiles
  tm <- seq.Date(from = as.Date(paste0(1946,"/",1,"/1")),to = as.Date(paste0(2150,"/",1,"/1")),by = "year")
  plotDF <- data.frame("Time" = tm)
  plotDF$VIp50 <- btc_prc.VI[,1]
  plotDF$VIp70 <- btc_prc.VI[,2]
  plotDF$VIp80 <- btc_prc.VI[,3]
  plotDF$VIp90 <- btc_prc.VI[,4]
  plotDF$VDp50 <- btc_prc.VD[,1]
  plotDF$VDp70 <- btc_prc.VD[,2]
  plotDF$VDp80 <- btc_prc.VD[,3]
  plotDF$VDp90 <- btc_prc.VD[,4]
  p <- plot_ly(plotDF)
  clr <- c('#e41a1c','#377eb8')
  cc <- clr[1]
  for (i in 2:(dim(plotDF)[2])) {
    if (i > 5){
      cc <- clr[2]
    }
    p <- add_trace(p, x = ~Time, y = plotDF[,i],  type = 'scatter', mode = 'lines',
                   name = names(plotDF)[i], line = list(color = cc, width = 4, dash = 'solid'))
  }
  p %>%
    layout(title = "Percentiles BTC" , #and 1% from Delta"
           xaxis = list(title = "Time"),
           yaxis = list(title = "[mg/l]"))
}

# Simulate for swat SWAT1 --------------
scenario$loadScen <- 'SWAT1'
scenario$wellType <- 'VI'
btc.VI <- mantis.Run(scenario)
scenario$wellType <- 'VD'
btc.VD <- mantis.Run(scenario)
btc_prc.VI <- colQuantiles(as.matrix(btc.VI), probs = c(.5, .7, 0.8, 0.9))
btc_prc.VD <- colQuantiles(as.matrix(btc.VD), probs = c(.5, .7, 0.8, 0.9))

# RASTER LOADING type --------------------------
# For the raster loading we need to know the background raster
library(hdf5r)

# Load the CVraster file which containts the correspondance between (row, col) and cell linear indices
cvraster <- H5File$new("../MantisData/v1_8/CVraster.h5",mode = "r")
CVR <- cvraster[["Raster"]]
CVR_D <- CVR[,]
IJ <- which(CVR_D != -1, arr.ind = T)
cvraster$close_all()
# the nth row of the raster loading corresponds to the IJ[n,] cell

# For Raster loading we have to specify the load scenario name
scenario$loadScen <- 'TEST01'
# and the subscenario name
scenario$loadSubScen <- 'test50'
# Similarly set the well type and run
scenario$wellType <- 'VI'
btc.VI <- mantis.Run(scenario)
scenario$wellType <- 'VD'
btc.VD <- mantis.Run(scenario)

# Set RASTER loading at Runtime ------------------

# Create new Raster Loading
user_load <- matrix(runif(dim(IJ)[1], min = 0, max = 200), nrow = dim(IJ)[1], ncol = 1)
userld.h5 <- H5File$new("userLoading.h5", mode = "w")
userld.h5[["Data"]] <- user_load
userld.h5[["Names"]] <- 'userload'
userld.h5$close_all()

scenario$modifierName <-'../../RmantisAPI/userLoading.h5'
scenario$modifierType <- 'Replace'
scenario$modifierUnit <- 'CONC'
scenario$wellType <- 'VI'

btc.VI <- mantis.Run(scenario)
scenario$wellType <- 'VD'
btc.VD <- mantis.Run(scenario)

# Create a multiplier Raster
mult_load <- matrix(runif(dim(IJ)[1], min = 0, max = 1), nrow = dim(IJ)[1], ncol = 1)
mult.h5 <- H5File$new("userMult.h5", mode = "w")
mult.h5$create_dataset(mult_load, name = "Data", dtype = h5types$H5T_NATIVE_FLOAT)
mult.h5[["Names"]] <- 'multload'
mult.h5$close_all()

scenario$modifierName <-'../../RmantisAPI/userMult.h5'
scenario$modifierType <- 'Multiply'
scenario$modifierUnit <- 'CONC'
scenario$wellType <- 'VI'
btc.VI <- mantis.Run(scenario)

scenario$wellType <- 'VD'
btc.VD <- mantis.Run(scenario)

# Test raster printed as integers
mult_load_int <- matrix(runif(dim(IJ)[1], min = 0, max = 100), nrow = dim(IJ)[1], ncol = 1)
mult_load_int <- round(mult_load_int)
multInt.h5 <- H5File$new("userMultInt.h5", mode = "w")
multInt.h5$create_dataset(mult_load_int, name = "Data", dtype = h5types$H5T_NATIVE_INT32)
multInt.h5[["Names"]] <- 'multloadInt'
multInt.h5$close_all()