#include <iostream>
#include <boost/mpi.hpp>

#include <chrono>

bool PrintMatrices = false;

#include "MS_input.h"
#include "SWAT_data.h"
#include "MS_raster.h"
#include "MS_mpi_utils.h"
#include "NPSAT_data.h"
#include "MS_unsat.h"
#include "MS_HRU_raster.h"
#include "MS_loading.h"
#include "MSdebug.h"



int main(int argc, char* argv[]) {
    boost::mpi::environment env( argc, argv );
    boost::mpi::communicator world;
    //{
    //    bool tf = false;
    //    std::cout << static_cast<int>(tf) << std::endl;
    //    tf = true;
    //    std::cout << static_cast<int>(tf) << std::endl;
    //    MS::testConvolution();
    //    MS::testBroadcast(world);
    //    MS::testMatrixBroadcast<double>(world);
    //    MS::testScalarBroadcast(world);
    //    return 0;
    //}



    MS::UserInput UI(world);

    bool tf = UI.read(argc,argv);
    if (!tf){
        return 0;
    }
    std::map<int, int> swatYRIDmap;
    MS::yearMap(swatYRIDmap, UI.simOptions.StartYear, UI.simOptions.StartYear + UI.simOptions.Nyears, UI.swatOptions.StartYear, UI.swatOptions.Nyears);

    auto startTotalSimulation = std::chrono::high_resolution_clock::now();
    if (world.rank() == 0){
        std::time_t result = std::time(nullptr);
        std::cout << "MantisSA started at " << std::asctime(std::localtime(&result)) << std::endl;
        std::cout << "MantisSA will run using " << world.size() << " processors" << std::endl;
        std::cout << "MantisSA version: " + UI.Version << std::endl;
    }
    world.barrier();


    std::ofstream dbg_file;
    if (UI.doDebug){
        dbg_file.open(UI.dbg_file.c_str());
        dbg_file << "Time, Eid, Sid, hruidx, gw_ratio, m_gw, v_gw, m_npsat, gw_conc, c_swat, UNshift" << std::endl;
    }


    //std::vector<int> VIids, VDids;
    std::map<int,MS::SelectedWellsGroup> SWGmap;
    if (UI.outputOptions.printSelectedWells){
        tf = MS::readSelectedWells(UI.outputOptions.SelectedWells, UI.outputOptions.SelectedWellGroups, SWGmap, world);
    }


    MS::BackgroundRaster backRaster;
    tf = backRaster.readData(UI.rasterOptions.File, UI.rasterOptions.Nrows, UI.rasterOptions.Ncols, UI.rasterOptions.Ncells, world);
    if (!tf){
        return 0;
    }

    MS::WELLS VI;
    MS::WELLS VD;
    {
        tf = MS::readNPSATdata(UI.npsatOptions.VIdataFile, VI, UI.simOptions.Porosity, UI.simOptions.Nyears,
                               backRaster, UI.riverOptions, world);
        if (!tf){return 0;}
        tf = MS::readNPSATdata(UI.npsatOptions.VDdataFile, VD, UI.simOptions.Porosity, UI.simOptions.Nyears,
                               backRaster, UI.riverOptions, world);
        if (!tf){return 0;}
    }

    MS::UNSAT UZ;
    tf = UZ.readdata(UI.unsatOptions.Depth_file, UI.unsatOptions.Depth_name, UI.unsatOptions.Rch_file,
                     UI.unsatOptions.wc, UI.unsatOptions.minDepth,UI.unsatOptions.minRch, world);
    if (!tf){ return 0;}

    MS::HRU_Raster hru_raster;
    tf = hru_raster.read(UI.swatOptions.HRU_raster_file, world);
    if (!tf){
        return 0;
    }
    MS::SWAT_data swat;
    tf = swat.read_HRU_idx_Map(UI.swatOptions.HRU_index_file, world);
    if (!tf){
        return 0;
    }
    tf = swat.read(UI.swatOptions.Data_file, UI.swatOptions.Nyears, UI.swatOptions.version, world);
    if (!tf){
        return 0;
    }

    MS::HistoricLoading HISTLOAD;
    if (!UI.historicOptions.filename.empty()){
        tf = HISTLOAD.Setup(UI.historicOptions.filename,
                            UI.historicOptions.ext,
                            world);
        if (!tf){
            return 0;
        }
    }


    //This is a map between Eid and cell that receiving water from this well
    // Only root processor knows this
    MS::WELL_CELLS well_cells;
    if (world.rank() == 0){
        bool tf = MS::readDistribPumping(UI.npsatOptions.DistribuPumpFile, well_cells);
        if (!tf){
            return 0;
        }
    }
    world.barrier();



    //std::vector<std::vector<int>> Eid_proc;
    //sendWellEids(VI, Eid_proc, world);

    if (UI.npsatOptions.bUseInitConcVI) {
        tf = MS::readInitSaltConc(UI.npsatOptions.InitSaltVIFile, VI, world.rank());
        if (!tf){return 0;}
    }
    if (UI.npsatOptions.bUseInitConcVD) {
        tf = MS::readInitSaltConc(UI.npsatOptions.InitSaltVDFile, VD, world.rank());
        if (!tf){return 0;}
    }


    std::cout << "Proc: " << world.rank() << " has " << VI.wells.size() << " VI wells and " << VD.wells.size() << " VD wells" << std::endl;

    //std::cout << "Here1" << std::endl;

    std::vector<double> ConcFromPump(backRaster.Ncell(), 0);
    std::vector<int> WellEidFromPump(backRaster.Ncell(), 0);
    // Initialize Concentration from pumping with the initial concentration
    world.barrier();

    std::map<int, MS::WELL>::iterator itw;
    if (UI.simOptions.EnableFeedback){
        if (world.rank() == 0){
            std::cout << "Initialize Pumping distribution matrix..." << std::endl;
        }
        std::vector<double> wellInitConc;
        for (itw = VI.wells.begin(); itw != VI.wells.end(); ++itw){
            wellInitConc.push_back(static_cast<double>(itw->first));
            wellInitConc.push_back(itw->second.initConc);
        }
        std::vector<std::vector<double>> AllwellsInitConc(world.size());
        MS::sendVec2Root<double>(wellInitConc, AllwellsInitConc, world);
        // Put the well concentrations in the right cells
        if (world.rank() == 0){
            MS::WELL_CELLS::iterator it; // well_cells
            for (int i = 0; i < world.size(); ++i){
                for (int j = 0; j < static_cast<int>(AllwellsInitConc[i].size()); j = j + 2){
                    int eid = static_cast<int>(AllwellsInitConc[i][j]);
                    it = well_cells.find(eid);
                    if (it != well_cells.end()){
                        for (int k = 0; k < it->second.size(); ++k){
                            ConcFromPump[it->second[k]] = AllwellsInitConc[i][j+1];
                            WellEidFromPump[it->second[k]] = eid;
                        }
                    }
                }
            }
        }
        world.barrier();
        MS::sendVectorFromRoot2AllProc(ConcFromPump, world);
        world.barrier();
        MS::sendVectorFromRoot2AllProc(WellEidFromPump, world);
        for (itw = VI.wells.begin(); itw != VI.wells.end(); ++itw){
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                int IJ = backRaster.IJ(itw->second.strml[i].urfI, itw->second.strml[i].urfJ);
                if (IJ >= 0){
                    itw->second.strml[i].wellSourceId = WellEidFromPump[IJ];
                }
            }
        }
        for (itw = VD.wells.begin(); itw != VD.wells.end(); ++itw){
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                int IJ = backRaster.IJ(itw->second.strml[i].urfI, itw->second.strml[i].urfJ);
                if (IJ >= 0){
                    itw->second.strml[i].wellSourceId = WellEidFromPump[IJ];
                }
            }
        }
    }
    world.barrier();


    // Main simulation loop

    //double min_neg_m = 0;

    //int hruidx;
    //std::vector<double> totMfeed(UI.NsimYears, 0);
    auto startTotal = std::chrono::high_resolution_clock::now();
    int YYYY = UI.simOptions.StartYear;
    int iswat = swatYRIDmap.find(YYYY)->second;
    for (int iyr = 0; iyr < UI.simOptions.Nyears; ++iyr){
        //std::cout << "Here1" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> wellConc;
        world.barrier();
        bool printThis = false;
        for (itw = VI.wells.begin(); itw != VI.wells.end(); ++itw){
            if (UI.doDebug){
                if (UI.dbg_id == itw->first){
                    printThis = true;
                }
            }
            //std::cout << itw->first << std::endl;
            std::vector<double> wellbtc(iyr + 1,0.0);
            double sumW = 0;
            int cntS = 0;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){


                //std::cout << i << " " << std::flush;
                //Calculate the C_SWAT
                double c_tot = 0.0;
                double c_gw = 0.0;


                MS::CreateLoadingForStreamline(c_tot, c_gw, iyr, iswat, YYYY, VI.version,
                                               itw->second.strml[i], UI, itw->second.initConc,
                                               HISTLOAD,backRaster, UZ, hru_raster, swat,
                                               ConcFromPump);


                itw->second.strml[i].gw_conc.push_back(c_gw);
                itw->second.strml[i].lf_conc.push_back(c_tot);

                std::vector<double> btc(itw->second.strml[i].lf_conc.size(), 0);
                std::vector<double> prebtc(itw->second.strml[i].lf_conc.size(), 0);
                MS::convolute(itw->second.strml[i].urf, itw->second.strml[i].lf_conc, btc, prebtc, itw->second.initConc);


                for (unsigned int j = 0; j < btc.size(); ++j){
                    wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * (btc[j] + prebtc[j]);
                    if (iyr == UI.simOptions.Nyears-1){
                        itw->second.strml[i].btc.push_back(btc[j] + prebtc[j]);
                    }
                }

                sumW = sumW + itw->second.strml[i].W;
                cntS = cntS + 1;

                //if (printThis){
                //    dbg_file << iyr << ", " << itw->first << ", " << itw->second.strml[i].Sid << ", " << hruidx << ", "
                //            //<< gw_ratio << ", " << m_gw << ", " << v_gw << ", " << m_npsat << ", " << gw_conc << ", "
                //            << c_swat << ", " << shift << std::endl;
                //}
            }// Loop streamlines

            if (cntS == 0){
                if (iyr == UI.simOptions.Nyears - 1 && UI.simOptions.OutofAreaConc < 0){
                    // if we ignore out of area streamlines we have to ignore the wells
                    // with streamlines that originate from outside. To avoid breaking the printing procedures
                    // we will print -9
                    for (unsigned int i = 0; i < wellbtc.size(); ++i){
                        wellbtc[i] = -9.0;
                    }
                }
                sumW = 1;
            }

            //std::cout << std::endl;
            if (iyr == UI.simOptions.Nyears - 1){
                for (unsigned int i = 0; i < wellbtc.size(); ++i){
                    itw->second.wellBtc.push_back(wellbtc[i] / sumW);
                }
            }
            else{
                wellConc.push_back(static_cast<double>(itw->first));
                wellConc.push_back(wellbtc.back()/sumW);
            }
            if (printThis){printThis = false;}
        }// Loop wells

        // Go through the domestic wells but built only this year concentration
        for (itw = VD.wells.begin(); itw != VD.wells.end(); ++itw){
            //std::cout << itw->first << std::endl;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){

                double c_tot = 0.0;
                double c_gw = 0.0;


                MS::CreateLoadingForStreamline(c_tot, c_gw, iyr, iswat, YYYY, VD.version,
                                               itw->second.strml[i], UI, itw->second.initConc,
                                               HISTLOAD,backRaster, UZ, hru_raster, swat,
                                               ConcFromPump);

                itw->second.strml[i].gw_conc.push_back(c_gw);
                itw->second.strml[i].lf_conc.push_back(c_tot);
            }
        }

        // Send this year pumped concentration to root processor
        if (iyr < UI.simOptions.Nyears-1){
            if (UI.simOptions.EnableFeedback) {
                std::vector<std::vector<double>> AllwellsConc(world.size());
                MS::sendVec2Root<double>(wellConc, AllwellsConc, world);
                // Put the well concentrations in the right cells
                if (world.rank() == 0){
                    MS::WELL_CELLS::iterator it; // well_cells
                    for (int i = 0; i < world.size(); ++i){
                        for (int j = 0; j < static_cast<int>(AllwellsConc[i].size()); j = j + 2){
                            int eid = static_cast<int>(AllwellsConc[i][j]);
                            it = well_cells.find(eid);
                            if (it != well_cells.end()){
                                for (int k = 0; k < it->second.size(); ++k){
                                    ConcFromPump[it->second[k]] = AllwellsConc[i][j+1];
                                }
                            }
                        }
                    }
                }
                world.barrier();
                MS::sendVectorFromRoot2AllProc(ConcFromPump, world);
                world.barrier();
            }
            world.barrier();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            if (world.rank() == 0){
                std::cout << "Year: " << iyr << " in " << elapsed.count() << std::endl;
            }
        }
        else{ // if this is the last year of the simulation calculate the domestic wells
            world.barrier();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            if (world.rank() == 0){
                std::cout << "Year: " << iyr << " in " << elapsed.count() << std::endl;
            }

            start = std::chrono::high_resolution_clock::now();
            if (world.rank() == 0){
                std::cout << "Simulating VD BTCs ..." << std::endl;
            }
            world.barrier();
            for (itw = VD.wells.begin(); itw != VD.wells.end(); ++itw){
                std::vector<double> wellbtc(UI.simOptions.Nyears, 0.0);
                double sumW = 0;
                int cntS = 0;
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    if (itw->second.strml[i].lf_conc.size() == 0){
                        continue;
                    }
                    std::vector<double> btc(itw->second.strml[i].lf_conc.size(), 0);
                    std::vector<double> prebtc(itw->second.strml[i].lf_conc.size(), 0);
                    MS::convolute(itw->second.strml[i].urf, itw->second.strml[i].lf_conc, btc, prebtc, itw->second.initConc);
                    //int shift = UZ.traveltime(itw->second.strml[i].IJ);

                    for (unsigned int j = 0; j < btc.size(); ++j){
                        wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * (btc[j] + prebtc[j]);
                        itw->second.strml[i].btc.push_back(btc[j] + prebtc[j]);
                    }
                    sumW = sumW + itw->second.strml[i].W;
                    cntS = cntS + 1;
                }
                if (cntS == 0){
                    if (iyr == UI.simOptions.Nyears-1 && UI.simOptions.OutofAreaConc < 0){
                        // if we ignore out of area streamlines we have to ignore the wells
                        // with streamlines that originate from outside. To avoid breaking the printing procedures
                        // we will print -9
                        for (unsigned int i = 0; i < wellbtc.size(); ++i){
                            wellbtc[i] = -9.0;
                        }
                    }
                    sumW = 1;
                }
                for (unsigned int i = 0; i < wellbtc.size(); ++i){
                    itw->second.wellBtc.push_back(wellbtc[i] / sumW);
                }
            }

            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish - start;
            if (world.rank() == 0){
                std::cout << "VD simulation time : " << elapsed.count() << std::endl;
            }

        }
        iswat = iswat + 1;
        if (iswat >= UI.swatOptions.Nyears){
            iswat = 0;
        }
        YYYY = YYYY + 1;
    }



    auto finishTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTotal = finishTotal - startTotal;
    if (world.rank() == 0){
        std::cout << "Total simulation for " << UI.simOptions.Nyears << " : " << elapsedTotal.count() / 60.0 << " min" << std::endl;
    }

    //{
    //    if (world.rank() == 0){
    //        std::string tmpname = "ConcPump1.dat";
    //        MS::printConcFromPump(tmpname, ConcFromPump);
    //    }
    //}

    { // Print the VI
        world.barrier();
        std::vector<double> thisProcBTCVI;
        MS::linearizeWellBTCs(VI, thisProcBTCVI);
        std::vector<std::vector<double>> AllProcBTCVI;
        MS::sendVec2Root<double>(thisProcBTCVI, AllProcBTCVI, world);
        if (world.rank() == 0){
            std::cout << "Printing VI BTCs ..." << std::endl;
            MS::printWELLSfromAllProc(AllProcBTCVI,UI.outputOptions.OutFile + "_VI.dat", UI.simOptions.Nyears);
        }
    }

    { // Print the VD
        world.barrier();
        std::vector<double> thisProcBTCVD;
        MS::linearizeWellBTCs(VD, thisProcBTCVD);
        std::vector<std::vector<double>> AllProcBTCVD;
        MS::sendVec2Root<double>(thisProcBTCVD, AllProcBTCVD, world);
        if (world.rank() == 0){
            std::cout << "Printing VD BTCs ..." << std::endl;
            MS::printWELLSfromAllProc(AllProcBTCVD,UI.outputOptions.OutFile + "_VD.dat",UI.simOptions.Nyears);
        }
    }

    //if (!VIids.empty()){ // Printing output for selected VI
        if (UI.outputOptions.printLoad){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                std::vector<double> thisProcMfeed;
                MS::BundleDetailData(VI, thisProcDATA, thisProcMfeed, it->second.idVI ,UZ, hru_raster, swat, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                std::vector<std::vector<double>> AllProcMfeed;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                MS::sendVec2Root<double>(thisProcMfeed, AllProcMfeed, world);
                if (world.rank() == 0){
                    std::cout << "Printing VI Salt load data for " << it->second.groupName << "..." << std::endl;
                    std::string lf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_lf.dat";
                    std::string mf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_mf.dat";
                    MS::printDetailOutputFromAllProc(AllProcDATA, lf_file_name, UI.simOptions.Nyears);
                    MS::printMfeedFromAllProc(AllProcMfeed, mf_file_name, UI.simOptions.Nyears);
                }
            }
        }

        if (UI.outputOptions.printURFs){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                MS::linearizeURFS(VI, thisProcDATA, it->second.idVI, swat, hru_raster, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                if (world.rank() == 0){
                    std::cout << "Printing VI URFs " << it->second.groupName << "..." << std::endl;
                    std::string urf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_urf.dat";
                    MS::printURFsFromAllProc(AllProcDATA, urf_file_name, UI.simOptions.Nyears);
                }
            }
        }

        if (UI.outputOptions.printBTCs){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                MS::linearizeBTCs(VI, thisProcDATA, it->second.idVI, swat, hru_raster, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                if (world.rank() == 0){
                    std::cout << "Printing VI BTCs " << it->second.groupName << "..." << std::endl;
                    std::string btc_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_btc.dat";
                    MS::printBTCssFromAllProc(AllProcDATA, btc_file_name, UI.simOptions.Nyears);
                }
            }
        }
    //}

    //if (!VDids.empty()){ // Printing detailed output for VD
        if (UI.outputOptions.printLoad){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                std::vector<double> thisProcMfeed;
                MS::BundleDetailData(VD, thisProcDATA, thisProcMfeed, it->second.idVD, UZ, hru_raster, swat, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                std::vector<std::vector<double>> AllProcMfeed;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                MS::sendVec2Root<double>(thisProcMfeed, AllProcMfeed, world);
                if (world.rank() == 0){
                    std::cout << "Printing VD Salt load data for " << it->second.groupName << "..." << std::endl;
                    std::string lf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_lf.dat";
                    std::string mf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_mf.dat";
                    MS::printDetailOutputFromAllProc(AllProcDATA,lf_file_name, UI.simOptions.Nyears);
                    MS::printMfeedFromAllProc(AllProcMfeed, mf_file_name, UI.simOptions.Nyears);
                }
            }
        }

        if (UI.outputOptions.printURFs){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                MS::linearizeURFS(VD, thisProcDATA, it->second.idVD, swat, hru_raster, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                if (world.rank() == 0){
                    std::cout << "Printing VD URFs " << it->second.groupName << "..." << std::endl;
                    std::string urf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_urf.dat";
                    MS::printURFsFromAllProc(AllProcDATA, urf_file_name, UI.simOptions.Nyears);
                }
            }
        }

        if (UI.outputOptions.printBTCs){
            world.barrier();
            std::map<int,MS::SelectedWellsGroup>::iterator it;
            for (it = SWGmap.begin(); it != SWGmap.end(); ++it){
                std::vector<double> thisProcDATA;
                MS::linearizeBTCs(VD, thisProcDATA, it->second.idVD, swat, hru_raster, UI.simOptions.Nyears);
                std::vector<std::vector<double>> AllProcDATA;
                MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
                if (world.rank() == 0){
                    std::cout << "Printing VD BTCs " << it->second.groupName << "..." << std::endl;
                    std::string btc_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_btc.dat";
                    MS::printBTCssFromAllProc(AllProcDATA, btc_file_name, UI.simOptions.Nyears);
                }

            }
        }
    //}

    if (UI.doDebug){
        dbg_file.close();
    }

    world.barrier();
    auto finishTotalSimulation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTotalSimulation = finishTotalSimulation - startTotalSimulation;

    world.barrier();
    if (world.rank() == 0){
        std::time_t result = std::time(nullptr);
        std::cout << "MantisSA Finished at " << std::asctime(std::localtime(&result)) << std::endl;
        std::cout << "Total simulation for : " << elapsedTotalSimulation.count() / 60.0 << " min" << std::endl;
    }

    return 0;
}
