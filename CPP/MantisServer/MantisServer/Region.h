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
        bool runSimulation(int threadid, int nThreads,
                           Scenario &scenario,
                           std::vector<std::string> &replymsg, int &nWellBTC);

    private:
        std::string path;
        BackgroundRaster raster;
        BMapCollection Bmaps;
        //LinearData Unsat;
        LinearData UnsDepth;
        RechargeScenarioList Rch;
        FlowWellCollection FWC;
        NLoadList NLL;
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
            //("Data.UNSAT", po::value<std::string>(), "The file with the unsaturated data")
            ("Data.UNSDEPTH", po::value<std::string>(), "The file with the Depth of unsaturated zone")
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
            //std::cout << "-------- Unsaturated travel time --------" << std::endl;
            //std::string unsatfile =  path + vm_ro["Data.UNSAT"].as<std::string>();
            //Unsat.setNoDataValue(0.0);
            //bool tf = Unsat.readData(unsatfile,raster.Ncell());
            //if (!tf)
            //    return false;

            std::cout << "-------- Unsaturated Depth --------" << std::endl;
            std::string unsDepthfile =  path + vm_ro["Data.UNSDEPTH"].as<std::string>();
            UnsDepth.setNoDataValue(0.0);
            bool tf = UnsDepth.readData(unsDepthfile,raster.Ncell());
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
            std::cout << "-------- URFs --------" << std::endl;std::string urfsfile =  vm_ro["Data.URFS"].as<std::string>();
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
            outmsg += "0 ERROR: The Region [" + scenario.modelArea + "] ";
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
        tf = FWC.hasFlowScenario(scenario.flowWellScen, scenario.porosity, scenario.porosityIndex);
        if (!tf){
            outmsg += "0 ERROR: The Region [" + scenario.modelArea + "] ";
            outmsg += "does not have the flow Scenario [" + scenario.flowScen + "]";
            return false;
        }

        //Test for recharge map
        if (scenario.bUseFlowRch){
            scenario.rchName = scenario.flowScen; //.substr(0,scenario.flowScen.size()-3);
        }
        tf = Rch.hasRechargeMaps(scenario.rchName);
        if (!tf){
            outmsg += "0 ERROR: The recharge [" + scenario.rchName + "]  could not be found";
            return false;
        }
        // Test for the Unsaturated option
        tf = UnsDepth.hasScenario((scenario.unsatScenario));
        if (!tf){
            outmsg += "0 ERROR: The Unsat [" + scenario.unsatScenario + "]  could not be found";
            return false;
        }
        else{
            scenario.unsatScenarioID = UnsDepth.ScenarioIndex(scenario.unsatScenario);
        }

        {// Loading validation
            tf = NLL.hasLoading(scenario.loadScen);
            if (!tf){
                outmsg += "0 ERROR: The Region [" + scenario.modelArea + "] ";
                outmsg += "does not have the Loading Scenario [" + scenario.loadScen + "]";
                return false;
            }
            if (scenario.loadScen.compare(scenario.LoadTransitionName) == 0){
                scenario.buseLoadTransition = false;
            }else{
                if (scenario.LoadTransitionName.compare("NONE") != 0){
                    tf = NLL.hasLoading(scenario.LoadTransitionName);
                    if (!tf){
                        outmsg += "0 ERROR: The Region [" + scenario.modelArea + "] ";
                        outmsg += "does not have the Loading Scenario [" + scenario.LoadTransitionName + "]";
                        return false;
                    }
                }
            }
        }

        return true;
    }

    bool Region::runSimulation(int threadid, int nThreads,
                               Scenario &scenario,
                               std::vector<std::string> &replymsg, int &nWellBTC){
        double dLTs = static_cast<double>(scenario.LoadTransitionStart);
        double dLTe = static_cast<double>(scenario.LoadTransitionEnd);
        double dSY = static_cast<double>(scenario.startSimulationYear);
        std::ofstream urf_file;
        std::ofstream lf_file;
        std::ofstream btc_file;
        std::ofstream well_btc_file;
        //Initialize debugging output files
        if (scenario.printAdditionalInfo){
            std::string root_name;
            root_name = scenario.debugPath + "_" + scenario.debugID + "_" + num2Padstr(threadid, 2);
            if (scenario.printLF){
                std::string lf_file_name = root_name + "_lf.dat";
                lf_file.open(lf_file_name.c_str());
            }
            if (scenario.printURF){
                std::string urf_file_name = root_name + "_urf.dat";
                urf_file.open(urf_file_name.c_str());
            }
            if (scenario.printBTC){
                std::string btc_file_name = root_name + "_btc.dat";
                btc_file.open(btc_file_name.c_str());
            }
            if (scenario.printWellBTC){
                std::string well_btc_file_name = root_name + "_well_btc.dat";
                well_btc_file.open(well_btc_file_name.c_str());
            }
        }

        std::vector<int> wellids;
        int NsimulationYears = scenario.endSimulationYear - scenario.startSimulationYear;
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
        nWellBTC = 0;
        for (int iw = startWell; iw < endWell; ++iw){
            wellit = fwcit->second.Wells.find(wellids[iw]);

            if (!wellit->second.bsimulateThis(scenario))
                continue;

            std::vector<double> WellBTC(NsimulationYears, 0);

            for (strmlit = wellit->second.streamlines.begin(); strmlit != wellit->second.streamlines.end(); ++strmlit){

                // If the streamline originates from stream it has zero loading
                if (strmlit->second.inRiver){
                    continue;
                }
                // If the travel time of the streamline is higher than 400 years the assume zero loading
                if (strmlit->second.age[scenario.porosityIndex] > 400){
                    continue;
                }

                std::vector<int> cell_lin_ind;
                std::vector<double> rch_val;
                std::vector<double> cln_rch;
                std::vector<double> depth_val;

                // Make a list of this streamline source area
                for (std::vector<cell>::iterator cellit = strmlit->second.SourceArea.begin(); cellit != strmlit->second.SourceArea.end(); ++cellit){
                    cell_lin_ind.push_back(cellit->lin_ind);
                    rchit->second.getValues(cellit->lin_ind, rch, clprc);
                    if (rch < scenario.minRecharge){
                        rch = scenario.minRecharge;
                    }
                    rch_val.push_back(rch);
                    cln_rch.push_back(clprc);

                    double depth = UnsDepth.getValue(scenario.unsatScenarioID,cellit->lin_ind);
                    //tau = tau * scenario.unsatZoneMobileWaterContent;
                    //tau = std::ceil(tau);
                    if (depth < scenario.unsatMinDepth){
                        depth = scenario.unsatMinDepth;
                    }
                    depth_val.push_back(depth);

                    if (rch_val.size() >= scenario.maxSourceCells){
                        break;
                    }
                }

                //Build the main load
                std::vector<double> mainload(NsimulationYears, 0);
                mainLoadit->second.buildLoadingFunction(cell_lin_ind, rch_val, cln_rch, depth_val, mainload, scenario);
                if (scenario.buseLoadTransition){
                    std::vector<double> preload(NsimulationYears, 0);
                    preLoadit->second.buildLoadingFunction(cell_lin_ind, rch_val, cln_rch, depth_val, preload, scenario);
                    // blend the two
                    double dYear = 0;
                    for (int i = 0; i < NsimulationYears; ++i){
                        if (i + scenario.startSimulationYear < scenario.LoadTransitionStart ){
                            mainload[i] = preload[i];
                        }
                        else if (i + scenario.startSimulationYear > scenario.LoadTransitionEnd){
                            break;
                        }
                        else{
                            double u = (dYear + dSY - dLTs)/(dLTe - dLTs);
                            mainload[i] = (1-u)*preload[i] + u*mainload[i];
                        }
                        dYear = dYear + 1.0;
                    }
                }

                if (scenario.printLF){
                    lf_file << wellit->first << " " << strmlit->first << " ";
                    for (int ii = 0; ii < NsimulationYears; ++ii)
                        lf_file << std::scientific << std::setprecision(10) << mainload[ii] << " ";
                    lf_file << std::endl;
                }

                //Build the URF
                URF urf;
                if (scenario.urfType == URFTYPE::LGNRM){
                    urf.init(NsimulationYears, strmlit->second.mu[scenario.porosityIndex],
                             strmlit->second.std[scenario.porosityIndex],URFTYPE::LGNRM);

                }
                else if (scenario.urfType == URFTYPE::ADE){
                    ADEoptions ade_opt;
                    ade_opt.lambda = scenario.adeLambda;
                    ade_opt.R = scenario.adeR;
                    urf.init(NsimulationYears, strmlit->second.len,
                             strmlit->second.age[scenario.porosityIndex],URFTYPE::ADE, ade_opt);
                }

                if (scenario.printURF){
                    urf_file << wellit->first << " " << strmlit->first << " ";
                    urf.print_urf(urf_file);
                }

                std::vector<double> BTC(NsimulationYears, 0);
                urf.convolute(mainload, BTC);

                if (scenario.printBTC){
                    btc_file << wellit->first << " " << strmlit->first << " "
                             << std::scientific << std::setprecision(10) << strmlit->second.w << " ";
                    for (int ii = 0; ii < NsimulationYears; ++ii)
                        btc_file << std::scientific << std::setprecision(10) << BTC[ii] << " ";
                    btc_file << std::endl;
                }

                for (unsigned int i = 0; i < BTC.size(); ++i){
                    WellBTC[i] = WellBTC[i] + BTC[i] * strmlit->second.w;
                }
            }
            if (scenario.printWellBTC){
                well_btc_file << wellit->first << " ";
                for (int ii = 0; ii < NsimulationYears; ++ii)
                    well_btc_file << std::scientific << std::setprecision(10) << WellBTC[ii] << " ";
                well_btc_file << std::endl;
            }

            if (scenario.printWellIds){
                replymsg[threadid] += std::to_string(wellit->first);
                replymsg[threadid] += " ";
            }

            for (unsigned int i = 0; i < WellBTC.size(); ++i){
                replymsg[threadid] += std::to_string(static_cast<float>(WellBTC[i]));
                replymsg[threadid] += " ";
            }
            nWellBTC = nWellBTC + 1;
        }

        if (scenario.printAdditionalInfo){
            if (scenario.printLF){lf_file.close();}
            if (scenario.printURF){urf_file.close();}
            if (scenario.printBTC){btc_file.close();}
            if (scenario.printWellBTC){well_btc_file.close();}
        }
        return true;
    }

}


#endif //MANTISSERVER_REGION_H
