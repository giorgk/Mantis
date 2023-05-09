//
// Created by giorg on 2/24/2023.
//

#ifndef MANTISSERVER_REGION_H
#define MANTISSERVER_REGION_H

#include <boost/program_options.hpp>

#include "BRaster.h"
#include "BMaps.h"
#include "wells.h"
#include "Nload.h"
#include "URFs.h"

namespace po = boost::program_options;

namespace mantisServer{
    class Region{
    public:
        Region(){}
        bool readRegionData(std::string inputfile);
        bool validateScenario(std::string &outmsg, Scenario &scenario);
        bool runSimulation(int threadid, Scenario &scenario);

    private:
        std::string path;
        BackgroundRaster raster;
        BMapCollection Bmaps;
        LinearData Unsat;
        RechargeScenarioList Rch;
        FlowWellCollection FWC;
        NLoadList NLL;
        int nThreads;
    };

    bool Region::readRegionData(std::string inputfile) {
        po::options_description regionOptions("regionOptions");
        po::variables_map vm_ro;
        regionOptions.add_options()
            // Raster Options
            ("Raster.Ncells", po::value<int>()->default_value(0), "Number of cells with index")
            ("Raster.Nrows", po::value<int>()->default_value(0), "Number of Raster rows")
            ("Raster.Ncols", po::value<int>()->default_value(0), "Number of Raster Columns")
            ("Raster.Xorig", po::value<double>()->default_value(0), "X coordinate of left lower raster corner")
            ("Raster.Yorig", po::value<double>()->default_value(0), "Y coordinate of left lower raster corner")
            ("Raster.CellSize", po::value<double>()->default_value(0), "Cell size of raster")
            ("Raster.File", po::value<std::string>(), "Name of the file with the raster values")

            // Data Options
            ("Data.BMAPS", po::value<std::string>(), "The name with the background maps")
            ("Data.NO3", po::value<std::string>(), "The main input file for Nitrate loading")
            ("Data.UNSAT", po::value<std::string>(), "The file with the unsaturated data")
            ("Data.WELLS", po::value<std::string>(), "The main input file for Wells")
            ("Data.URFS", po::value<std::string>(), "The main input file for URFs")
            ("Data.RCH", po::value<std::string>(), "The main input file for Recharge")
            ("Data.RFURF", po::value<std::string>(), "The main input file for RFURF")
            ("Data.Path", po::value<std::string>(), "The Relative path for the data")
        ;

        po::store(po::parse_config_file<char>(inputfile.c_str(), regionOptions,true), vm_ro);

       path =  vm_ro["Data.Path"].as<std::string>();

        {// Read Raster
            int Ncells = vm_ro["Raster.Ncells"].as<int>();
            int Nr = vm_ro["Raster.Nrows"].as<int>();
            int Nc = vm_ro["Raster.Ncols"].as<int>();
            double Xorig = vm_ro["Raster.Xorig"].as<double>();
            double Yorig = vm_ro["Raster.Yorig"].as<double>();
            double cs = vm_ro["Raster.CellSize"].as<double>();
            std::string rasterfile =  path + vm_ro["Raster.File"].as<std::string>();
            std::cout << "-------- Raster --------" << std::endl;
            raster.setGridLocation(Xorig, Yorig, cs);
            bool tf = raster.readData(rasterfile, Nr, Nc, Ncells);
            if (!tf)
                return false;
        }

        {// Read Background Maps
            std::cout << "-------- Background map --------" << std::endl;
            std::string bmapsfile =  path + vm_ro["Data.BMAPS"].as<std::string>();
            bool tf = Bmaps.readData(bmapsfile);
            if (!tf)
                return false;
        }

        {//Read Unsaturated data
            std::cout << "-------- Unsaturated travel time --------" << std::endl;
            std::string unsatfile =  path + vm_ro["Data.UNSAT"].as<std::string>();
            Unsat.setNoDataValue(0.0);
            bool tf = Unsat.readData(unsatfile,raster.Ncell());
            if (!tf)
                return false;
        }

        {//Read Recharge
            std::cout << "-------- Recharge --------" << std::endl;
            std::string rchfile =  vm_ro["Data.RCH"].as<std::string>();
            bool tf = Rch.readData(path,rchfile, raster.Ncell());
            if (!tf)
                return false;
        }


        {// Read the wells
            std::cout << "-------- Wells --------" << std::endl;
            std::string wellfile =  vm_ro["Data.WELLS"].as<std::string>();
            bool tf = FWC.readMainfile(path, wellfile, true, Bmaps);
            if (!tf)
                return false;
        }

        {// Read the Urfs
            std::cout << "-------- URFs --------" << std::endl;
            std::string urfsfile =  vm_ro["Data.URFS"].as<std::string>();
            bool tf = FWC.readMainfile(path, urfsfile, false, Bmaps);
            if (!tf)
                return false;
            //FWC.FlowScenarios["Modesto"].Wells[1].streamlines[0].print();
        }

        {// Run postprocess routines for the wells
            std::cout << "-------- Well sources --------" << std::endl;
            FWC.calcWellWeights();
            FWC.calcWellSourceArea(raster, Rch);
        }


        {// Read the loading maps
            std::cout << "-------- NO3 Loading --------" << std::endl;
            std::string nLoadfile =  vm_ro["Data.NO3"].as<std::string>();
            bool tf = NLL.readData(path, nLoadfile);
            if (!tf)
                return false;
        }

        return true;
    }

    bool Region::validateScenario(std::string &outmsg, Scenario &scenario) {
        bool tf;

        // Test for background map
        tf = Bmaps.hasMap(scenario.mapID);
        if (!tf){
            outmsg += "0 ERROR: The Region [" + scenario.region + "] ";
            outmsg += "does not have the map [" + scenario.mapID + "]";
            return false;
        }

        scenario.flowWellScen = scenario.flowScen + '_' + scenario.wellType;
        tf = Bmaps.validateMapsWells(scenario.mapID,scenario.regionIDs,scenario.flowWellScen, outmsg);
        if (!tf){
            outmsg = "0 ERROR: " + outmsg;
            return false;
        }

        //Test for the wells under the selected flow scenario
        tf = FWC.hasFlowScenario(scenario.flowScen);
        if (!tf){
            outmsg += "0 ERROR: The Region [" + scenario.region + "] ";
            outmsg += "does not have the flow Scenario [" + scenario.flowScen + "]";
            return false;
        }

        //Test for recharge map
        tf = Rch.hasRechargeMaps(scenario.rchName);
        if (!tf){
            outmsg += "0 ERROR: The recharge [" + scenario.rchName + "]  could not be found";
            return false;
        }
        tf = Unsat.hasScenario((scenario.unsatScenario));
        if (!tf){
            outmsg += "0 ERROR: The Unsat [" + scenario.unsatScenario + "]  could not be found";
            return false;
        }

        return true;
    }

    bool Region::runSimulation(int threadid, Scenario &scenario){
        std::vector<int> wellids;
        int NsimulationYears = scenario.endSimulationYear - 1945;
        Bmaps.getWells(scenario.mapID, scenario.regionIDs, scenario.flowWellScen, wellids);
        int startWell = 0;
        int endWell = 0;
        double rch,clprc;
        getStartEndIndices(threadid, nThreads,wellids.size(),startWell, endWell);

        // Get an iterator to the flowCollection
        std::map<std::string ,WellList>::iterator fwcit;
        std::map<int, Well>::iterator wellit;
        std::map<int, Streamline>::iterator strmlit;
        std::map<std::string, NLoad>::iterator mainLoadit, preLoadit;
        mainLoadit = NLL.NLoadMaps.find(scenario.loadScen);
        if (scenario.buseLoadTransition){
            preLoadit = NLL.NLoadMaps.find(scenario.LoadTransitionName);
        }
        std::map<std::string, RechargeScenario>::iterator rchit = Rch.RechargeList.find(scenario.rchName);




        fwcit = FWC.FlowScenarios.find(scenario.flowWellScen);
        for (int iw = startWell; iw < endWell; ++iw){
            wellit = fwcit->second.Wells.find(wellids[iw]);
            std::vector<double> WellBTC(NsimulationYears, 0);

            for (strmlit = wellit->second.streamlines.begin(); strmlit != wellit->second.streamlines.end(); ++ strmlit){
                std::vector<int> cell_lin_ind;
                std::vector<double> rch_val;
                std::vector<double> cln_rch;

                for (std::map<int,cell>::iterator cellit = strmlit->second.SourceArea.begin(); cellit != strmlit->second.SourceArea.end(); ++cellit){
                    cell_lin_ind.push_back(cellit->first);
                    rchit->second.getValues(cellit->first, rch, clprc);
                    rch_val.push_back(rch);
                    cln_rch.push_back(clprc);
                }


            }


        }

    }

}


#endif //MANTISSERVER_REGION_H
