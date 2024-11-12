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


    MS::HRU_Raster hru_raster;
    tf = hru_raster.read(UI.HRU_raster_file, world);
    if (!tf){
        return 0;
    }


    std::ofstream dbg_file;
    if (UI.doDebug){
        dbg_file.open(UI.dbg_file.c_str());
        dbg_file << "Time, Eid, Sid, hruidx, gw_ratio, m_gw, v_gw, m_npsat, Mfeed, c_swat, UNshift" << std::endl;
    }


    std::vector<int> VIids, VDids;
    if (!UI.SelectedWells_file.empty()){
        MS::readSelectedWells(UI.SelectedWells_file, VIids, VDids, world);
    }


    MS::BackgroundRaster backRaster;
    tf = backRaster.readData(UI.rasteroptions.File,UI.rasteroptions.Nrows,UI.rasteroptions.Ncols,UI.rasteroptions.Ncells, world);
    if (!tf){
        return 0;
    }

    MS::UNSAT UZ;
    tf = UZ.readdata(UI.depth_input_file,UI.depth_name,UI.rch_input_file,
                     UI.wc, UI.minDepth,UI.minRch, world);
    if (!tf){ return 0;}

    MS::WELLS VI;
    MS::WELLS VD;
    {
        tf = MS::readNPSATdata(UI.npsat_VI_file, VI, UI.porosity, UI.NsimYears, backRaster, world);
        if (!tf){return 0;}
        tf = MS::readNPSATdata(UI.npsat_VD_file, VD, UI.porosity, UI.NsimYears, backRaster, world);
        if (!tf){return 0;}
    }

    MS::SWAT_data swat;
    tf = swat.read_HRU_idx_Map(UI.hru_idx_file, world);
    tf = swat.read(UI.swat_input_file, UI.NswatYears, world);
    if (!tf){
        return 0;
    }



    //This is a map between Eid and cell that receiving water from this well
    // Only root processor knows this
    MS::WELL_CELLS well_cells;
    if (world.rank() == 0){
        bool tf = MS::readDistribPumping(UI.cell_well_file, well_cells);
        if (!tf){
            return 0;
        }
    }
    world.barrier();



    std::vector<std::vector<int>> Eid_proc;
    sendWellEids(VI, Eid_proc, world);


    tf = MS::readInitSaltConc(UI.init_salt_VI_file, VI);
    if (!tf){return 0;}
    tf = MS::readInitSaltConc(UI.init_salt_VD_file, VD);
    if (!tf){return 0;}

    std::cout << "Proc: " << world.rank() << " has " << VI.size() << " VI wells and " << VD.size() << " VD wells" << std::endl;

    //std::cout << "Here1" << std::endl;

    std::vector<double> ConcFromPump(backRaster.Ncell(), 0);

    // Main simulation loop
    MS::WELLS ::iterator itw;
    //int hruidx;
    //std::vector<double> totMfeed(UI.NsimYears, 0);
    auto startTotal = std::chrono::high_resolution_clock::now();
    int iswat = 0;
    for (int iyr = 0; iyr < UI.NsimYears; ++iyr){
        //std::cout << "Here1" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> wellConc;
        world.barrier();
        bool printThis = false;
        for (itw = VI.begin(); itw != VI.end(); ++itw){
            if (UI.doDebug){
                if (UI.dbg_id == itw->first){
                    printThis = true;
                }
            }
            //std::cout << itw->first << std::endl;
            std::vector<double> wellbtc(iyr + 1,0.0);
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                //std::cout << i << " " << std::flush;
                //hruidx = itw->second.strml[i].hru_idx - 1;
                int hru = hru_raster.getHRU(itw->second.strml[i].IJ);
                int hruidx = swat.hru_index(hru);
                if (hruidx < 0){
                    continue;
                }

                // Calculate the ratio of groundwater to the total
                double gw_ratio = 0.0;
                if (std::abs(swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]) > 0.00000001){
                    gw_ratio = swat.irrGW_mm[iswat][hruidx] / (swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]);
                }

                // irrigated Salt Mass. The Salts that come from the pumped water
                double m_gw = swat.irrsaltGW_Kgha[iswat][hruidx];

                //Calculate the C_SWAT
                double c_swat = 0.0;

                if (swat.perc_mm[iswat][hruidx] < 0.01){
                    itw->second.strml[i].Mfeed.push_back(0.0);
                    c_swat = swat.Salt_perc_ppm[iswat][hruidx];
                }
                else{
                    //Calculate the percolated groundwater volume
                    double v_gw = swat.perc_mm[iswat][hruidx] * gw_ratio;

                    // Calculate concentration from NPSAT spreading (Feedback)
                    double m_npsat = 0.0;
                    if (itw->second.strml[i].IJ > 0){
                        m_npsat = ConcFromPump[itw->second.strml[i].IJ] * v_gw / 100.0;
                    }

                    double Mfeed = m_npsat - m_gw;
                    if (Mfeed < 0) {
                        Mfeed = 0.0;
                    }
                    itw->second.strml[i].Mfeed.push_back(Mfeed);

                    c_swat = (swat.totpercsalt_kgha[iswat][hruidx] + Mfeed) * 100 / swat.perc_mm[iswat][hruidx];
                }

                if (c_swat > UI.maxConc){
                    c_swat = UI.maxConc;
                }

                if (UI.SurfConcValue > 0){
                    double surf_perc = UZ.getSurfPerc(itw->second.strml[i].IJ);
                    if (surf_perc > 0){
                        c_swat = c_swat*(1-surf_perc) + surf_perc*UI.SurfConcValue;
                    }
                }

                itw->second.strml[i].lf.push_back(c_swat);

                std::vector<double> btc(itw->second.strml[i].lf.size(), 0);
                std::vector<double> prebtc(itw->second.strml[i].lf.size(), 0);
                MS::convolute(itw->second.strml[i].urf, itw->second.strml[i].lf, btc, prebtc, itw->second.initConc);
                int shift = UZ.traveltime(itw->second.strml[i].IJ);

                for (unsigned int j = 0; j < btc.size(); ++j){
                    if (j <= shift){
                        wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * itw->second.initConc;
                    }
                    else{
                        wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * (btc[j-shift] + prebtc[j-shift]);
                    }
                }

                if (printThis){
                    dbg_file << iyr << ", " << itw->first << ", " << itw->second.strml[i].Sid << ", " << hruidx << ", "
                             //<< gw_ratio << ", " << m_gw << ", " << v_gw << ", " << m_npsat << ", " << Mfeed << ", "
                             << c_swat << ", " << shift << std::endl;
                }



            }// Loop streamlines
            //std::cout << std::endl;
            if (iyr == UI.NsimYears-1){
                for (unsigned int i = 0; i < wellbtc.size(); ++i){
                    itw->second.btc.push_back(wellbtc[i]);
                }
            }
            else{
                wellConc.push_back(wellbtc.back());
            }
            if (printThis){printThis = false;}


        }// Loop wells

        // Go through the domestic wells but built only this year concentration
        for (itw = VD.begin(); itw != VD.end(); ++itw){
            //std::cout << itw->first << std::endl;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                //std::cout << i << " " << std::flush;
                //hruidx = itw->second.strml[i].hru_idx - 1;
                int hru = hru_raster.getHRU(itw->second.strml[i].IJ);
                int hruidx = swat.hru_index(hru);
                if (hruidx < 0){
                    continue;
                }

                // Calculate the groundwater concentration
                double gw_ratio = 0.0;
                if (std::abs(swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]) > 0.00000001){
                    gw_ratio = swat.irrGW_mm[iswat][hruidx] / (swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]);
                }

                // irrigated Salt Mass
                double m_gw = swat.irrsaltGW_Kgha[iswat][hruidx];

                //Calculate the C_SWAT
                double c_swat = 0.0;

                if (swat.perc_mm[iswat][hruidx] < 0.01){
                    itw->second.strml[i].Mfeed.push_back(0.0);
                    c_swat = swat.Salt_perc_ppm[iswat][hruidx];
                }
                else{
                    //Calculate the percolated groundwater volume
                    double v_gw = swat.perc_mm[iswat][hruidx] * gw_ratio;

                    // Calculate concentration from NPSAT spreading (Feedback)
                    double m_npsat = 0.0;
                    if (itw->second.strml[i].IJ > 0){
                        m_npsat = ConcFromPump[itw->second.strml[i].IJ] * v_gw / 100.0;
                    }

                    double Mfeed = m_npsat - m_gw;
                    if (Mfeed < 0) {
                        Mfeed = 0.0;
                    }
                    itw->second.strml[i].Mfeed.push_back(Mfeed);

                    c_swat = (swat.totpercsalt_kgha[iswat][hruidx] + Mfeed) * 100 / swat.perc_mm[iswat][hruidx];
                }

                if (c_swat > UI.maxConc){
                    c_swat = UI.maxConc;
                }

                if (UI.SurfConcValue > 0){
                    double surf_perc = UZ.getSurfPerc(itw->second.strml[i].IJ);
                    if (surf_perc > 0){
                        c_swat = c_swat*(1-surf_perc) + surf_perc*UI.SurfConcValue;
                    }
                }

                itw->second.strml[i].lf.push_back(c_swat);
            }
        }

        // Send this year pumped concentration to root processor
        if (iyr < UI.NsimYears-1){
            std::vector<std::vector<double>> AllwellsConc(world.size());
            MS::sendVec2Root<double>(wellConc, AllwellsConc, world);
            // Put the well concentrations in the right cells
            if (world.rank() == 0){
                MS::WELL_CELLS::iterator it; // well_cells
                for (int i = 0; i < world.size(); ++i){
                    for (int j = 0; j < static_cast<int>(AllwellsConc[i].size()); ++j){
                        it = well_cells.find(Eid_proc[i][j]);
                        if (it != well_cells.end()){
                            for (int k = 0; k < it->second.size(); ++k){
                                ConcFromPump[it->second[k]] = AllwellsConc[i][j];
                            }
                        }
                    }
                }
            }
            world.barrier();
            MS::sendVectorFromRoot2AllProc(ConcFromPump, world);
            world.barrier();
        }
        else{ // if this is the last year of the simulation calculate the domestic wells
            if (world.rank() == 0){
                std::cout << "Simulating VD BTCs ..." << std::endl;
            }
            world.barrier();
            for (itw = VD.begin(); itw != VD.end(); ++itw){
                std::vector<double> wellbtc(UI.NsimYears, 0.0);
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    std::vector<double> btc(itw->second.strml[i].lf.size(), 0);
                    std::vector<double> prebtc(itw->second.strml[i].lf.size(), 0);
                    MS::convolute(itw->second.strml[i].urf, itw->second.strml[i].lf, btc, prebtc, itw->second.initConc);
                    int shift = UZ.traveltime(itw->second.strml[i].IJ);

                    for (unsigned int j = 0; j < btc.size(); ++j){
                        if (j <= shift){
                            wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * itw->second.initConc;
                        }
                        else{
                            wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * (btc[j-shift] + prebtc[j-shift]);
                        }
                    }
                }
                for (unsigned int i = 0; i < wellbtc.size(); ++i){
                    itw->second.btc.push_back(wellbtc[i]);
                }
            }
        }

        world.barrier();
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        if (world.rank() == 0){
            std::cout << "Year: " << iyr << " in " << elapsed.count() << std::endl;
        }
        iswat = iswat + 1;
        if (iswat >= UI.NswatYears){
            iswat = 0;
        }
    }



    auto finishTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTotal = finishTotal - startTotal;
    if (world.rank() == 0){
        std::cout << "Total simulation for " << UI.NsimYears << " : " << elapsedTotal.count() / 60.0 << " min" << std::endl;
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
        MS::linearizeBTC(VI, thisProcBTCVI);
        std::vector<std::vector<double>> AllProcBTCVI;
        MS::sendVec2Root<double>(thisProcBTCVI, AllProcBTCVI, world);
        if (world.rank() == 0){
            std::cout << "Printing VI BTCs ..." << std::endl;
            MS::printWELLSfromAllProc(AllProcBTCVI,UI.outfileVI,UI.NsimYears);
        }
    }

    { // Print the VD
        world.barrier();
        std::vector<double> thisProcBTCVD;
        MS::linearizeBTC(VD, thisProcBTCVD);
        std::vector<std::vector<double>> AllProcBTCVD;
        MS::sendVec2Root<double>(thisProcBTCVD, AllProcBTCVD, world);
        if (world.rank() == 0){
            std::cout << "Printing VD BTCs ..." << std::endl;
            MS::printWELLSfromAllProc(AllProcBTCVD,UI.outfileVD,UI.NsimYears);
        }
    }

    if (!VIids.empty()){ // Printing detailed output for VI
        world.barrier();
        std::vector<double> thisProcDATA;
        std::vector<double> thisProcMfeed;
        MS::BundleDetailData(VI, thisProcDATA, thisProcMfeed, VIids,UZ, hru_raster, swat, UI.NsimYears);
        std::vector<std::vector<double>> AllProcDATA;
        std::vector<std::vector<double>> AllProcMfeed;
        MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
        MS::sendVec2Root<double>(thisProcMfeed, AllProcMfeed, world);
        if (world.rank() == 0){
            std::cout << "Printing VI Detailed data ..." << std::endl;
            MS::printDetailOutputFromAllProc(AllProcDATA,UI.outfileVIdetail,UI.NsimYears);
            MS::printMfeedFromAllProc(AllProcMfeed,UI.outfileVImfeed,UI.NsimYears);
        }
    }

    if (!VDids.empty()){ // Printing detailed output for VD
        world.barrier();
        std::vector<double> thisProcDATA;
        std::vector<double> thisProcMfeed;
        MS::BundleDetailData(VD, thisProcDATA, thisProcMfeed, VDids,UZ, hru_raster, swat, UI.NsimYears);
        std::vector<std::vector<double>> AllProcDATA;
        std::vector<std::vector<double>> AllProcMfeed;
        MS::sendVec2Root<double>(thisProcDATA, AllProcDATA, world);
        MS::sendVec2Root<double>(thisProcMfeed, AllProcMfeed, world);
        if (world.rank() == 0){
            std::cout << "Printing VD Detailed data ..." << std::endl;
            MS::printDetailOutputFromAllProc(AllProcDATA,UI.outfileVDdetail,UI.NsimYears);
            MS::printMfeedFromAllProc(AllProcMfeed,UI.outfileVDmfeed,UI.NsimYears);
        }
    }

    if (UI.doDebug){
        dbg_file.close();
    }


    return 0;
}
