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

namespace po = boost::program_options;

namespace mantisServer{
    class Region{
    public:
        Region(){}
        bool readRegionData(std::string inputfile);

    private:
        std::string path;
        BackroundRaster raster;
        BMapCollection Bmaps;
        LinearData Unsat;
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
            double Xorig = vm_ro["Raster.Xorig"].as<int>();
            double Yorig = vm_ro["Raster.Yorig"].as<int>();
            double cs = vm_ro["Raster.CellSize"].as<int>();
            std::string rasterfile =  path + vm_ro["Raster.File"].as<std::string>();
            std::cout << "      Reading base Raster ..." << std::endl;
            raster.setGridLocation(Xorig, Yorig, cs);
            bool tf = raster.readData(rasterfile, Nr, Nc, Ncells);
            if (!tf)
                return false;
        }

        {// Read Background Maps
            std::cout << "      Reading Background map data  ..." << std::endl;
            std::string bmapsfile =  path + vm_ro["Data.BMAPS"].as<std::string>();
            bool tf = Bmaps.readData(bmapsfile);
            if (!tf)
                return false;
        }

        {//Read Unsaturated data
            std::cout << "      Reading Unsaturated zone data..." << std::endl;
            std::string unsatfile =  path + vm_ro["Data.UNSAT"].as<std::string>();
            Unsat.setNoDataValue(0.0);
            bool tf = Unsat.readData(unsatfile,raster.Ncell());
            if (!tf)
                return false;
        }

        {//Read Recharge
            std::cout << "      Reading data for Recharge ..." << std::endl;
            std::string rchfile =  vm_ro["Data.RCH"].as<std::string>();
            bool tf = Rch.readData(path,rchfile, raster.Ncell());
            if (!tf)
                return false;
        }


        {// Read the wells
            std::cout << "      Reading well data ..." << std::endl;
            std::string wellfile =  vm_ro["Data.WELLS"].as<std::string>();
            bool tf = FWC.readMainfile(path, wellfile, true, Bmaps);
            if (!tf)
                return false;
        }

        {// Read the Urfs
            std::cout << "      Reading URF data ..." << std::endl;
            std::string urfsfile =  vm_ro["Data.URFS"].as<std::string>();
            bool tf = FWC.readMainfile(path, urfsfile, false, Bmaps);
            if (!tf)
                return false;
            //FWC.FlowScenarios["Modesto"].Wells[1].streamlines[0].print();
        }

        {// Run postprocess routines for the wells
            FWC.calcWellWeights();
            FWC.calcWellSourceArea(raster, Rch);
        }


        {// Read the loading maps
            std::string nLoadfile =  vm_ro["Data.NO3"].as<std::string>();
            bool tf = NLL.readData(path, nLoadfile);
            if (!tf)
                return false;
        }

        return true;
    }

}


#endif //MANTISSERVER_REGION_H
