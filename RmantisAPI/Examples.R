# All the R functions live under the gwtools (https://github.com/UCD-GW-Nitrate/gwtools) package.
# Either install and load the entire package or download only the mantis.R file and source it.
# Here I'll load the package
library(gwtools)

#------------- Setup server and path -------------------------------------
# First launch the server.
# Open a terminal and run MantisServer.exe -c config.ini

# ----------Configure the scenario options ----------------------------
# The following command creates a list with all available options to configure a scenario
scenario <- mantis.ScenarioOptions()

# First setup the path for the Mantis client program
scenario$client <- "../CPP/TestClient/TestClient/cmake-build-release/MantisClient.exe"

# -----------Quit Server ------------------
res <- mantis.Quit(scenario)

# Configure the input file name. This is the file where all the options will be written.
# This file name will be used as input to the Mantis client program.
scenario$infile <- 'inputMessage.dat'

# We need to setup a file name for the output file where the Mantis client will print the results
scenario$outfile <- 'BTCresults.dat'

# -------- Test the connection ----------------
res <- mantis.Ping(scenario)



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
scenario$Regions <- c('Merced')

# Select Nitrate Loading -------------------------------
# The table provides the available options https://github.com/giorgk/Mantis/wiki/Runtime-Data#nitrate-loading-scenarios
# Let's select the GNLM
scenario$loadScen <- 'GNLM'

# Using the GNLM and SWAT loading type we can modify the Nitrate loading by selecting crop types and reduction rates
# However here we will run the base loading by setting to all crops (id = -9) a reduction of 1
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
# 4 The server will run the simulation ans send the results to client
# 5 The client will read the results and return them to the R workspace
# All this sequence is executed with one go as follows
btc <- mantis.Run(scenario)

# Switch to domestic wells and repeat the simulation
scenario$wellType <- 'VD'
btc1 <- mantis.Run(scenario)

# To calculate percentiles we an use the package matrixStats
library(matrixStats)
btc_prc <- colQuantiles(as.matrix(btc), probs = c(.25, .5, .75, 0.85, 0.95))
btc_prc1 <- colQuantiles(as.matrix(btc1), probs = c(.25, .5, .75, 0.85, 0.95))

# For the plot I use plotly which makes nice interactive plots
library(plotly)



