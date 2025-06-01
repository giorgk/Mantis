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

    MS::UNSAT UZ;
    tf = UZ.readdata(UI.unsatOptions.Depth_file, UI.unsatOptions.Depth_name, UI.unsatOptions.Rch_file,
                     UI.unsatOptions.wc, UI.unsatOptions.minDepth,UI.unsatOptions.minRch, world);
    if (!tf){ return 0;}

    MS::WELLS VI;
    MS::WELLS VD;
    {
        tf = MS::readNPSATdata(UI.npsatOptions.VIdataFile, VI, UI.simOptions.Porosity, UI.simOptions.Nyears,
                               UI.npsatOptions.version, backRaster, world);
        if (!tf){return 0;}
        tf = MS::readNPSATdata(UI.npsatOptions.VDdataFile, VD, UI.simOptions.Porosity, UI.simOptions.Nyears,
                               UI.npsatOptions.version, backRaster, world);
        if (!tf){return 0;}
    }

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

    MS::HistoricLoading HIST;
    if (!UI.historicOptions.filename.empty()){
        tf = HIST.Setup(UI.historicOptions.StartYear,
                       UI.historicOptions.EndYear,
                       UI.historicOptions.Interval,
                       UI.historicOptions.filename,
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


    std::cout << "Proc: " << world.rank() << " has " << VI.size() << " VI wells and " << VD.size() << " VD wells" << std::endl;

    //std::cout << "Here1" << std::endl;

    std::vector<double> ConcFromPump(backRaster.Ncell(), 0);
    std::vector<int> WellEidFromPump(backRaster.Ncell(), 0);
    // Initialize Concentration from pumping with the initial concentration
    world.barrier();

    MS::WELLS ::iterator itw;
    if (UI.simOptions.EnableFeedback){
        if (world.rank() == 0){
            std::cout << "Initialize Pumping distribution matrix..." << std::endl;
        }
        std::vector<double> wellInitConc;
        for (itw = VI.begin(); itw != VI.end(); ++itw){
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
        for (itw = VI.begin(); itw != VI.end(); ++itw){
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                int IJ = backRaster.IJ(itw->second.strml[i].urfI, itw->second.strml[i].urfJ);
                if (IJ >= 0){
                    itw->second.strml[i].wellSourceId = WellEidFromPump[IJ];
                }
            }
        }
        for (itw = VD.begin(); itw != VD.end(); ++itw){
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
    int iswat = 0;
    int YYYY = UI.simOptions.StartYear;
    for (int iyr = 0; iyr < UI.simOptions.Nyears; ++iyr){
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
            double sumW = 0;
            int cntS = 0;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){


                //std::cout << i << " " << std::flush;
                //Calculate the C_SWAT
                double c_tot = 0.0;
                double c_gw = 0.0;
                double c_fut = 0.0;

                MS::CreateLoadingForStreamline(c_fut, c_gw, iyr, iswat,
                                               itw->second.strml[i], UI, itw->second.initConc,
                                               backRaster, UZ, hru_raster, swat,
                                               ConcFromPump);

                if (YYYY < UI.historicOptions.BlendEnd){
                    int IJ = backRaster.IJ(itw->second.strml[i].urfI, itw->second.strml[i].urfJ);
                    if (IJ >= 0){
                        double c_hist = HIST.calculateConc(IJ, YYYY, world);
                        double u = 0.0;
                        MS::BlendLoading(c_tot, u, UI.historicOptions, c_hist, c_fut, YYYY);
                        c_gw = c_gw * u;
                    }
                }
                else{
                    c_tot = c_fut;
                }

                /*
                int hruidx = -9;
                double Mfeed = 0.0;
                bool addThis = false;

                if (itw->second.strml[i].inRiv){
                    c_swat = UI.riverOptions.ConcValue;
                    addThis = true;
                }
                else if (itw->second.strml[i].a > UI.simOptions.MaxAge){
                    c_swat = itw->second.initConc;
                    addThis = true;
                }
                else{
                    double n_cswat = 0.0;
                    for (int ii = itw->second.strml[i].urfI - UI.simOptions.nBuffer; ii <= itw->second.strml[i].urfI + UI.simOptions.nBuffer; ++ii){
                        for(int jj = itw->second.strml[i].urfJ - UI.simOptions.nBuffer; jj <= itw->second.strml[i].urfJ + UI.simOptions.nBuffer; ++jj){
                            double cell_cswat = 0.0;
                            double cell_mfeed = 0.0;
                            double cell_cgw = 0.0;
                            int IJ = backRaster.IJ(ii, jj);
                            if (IJ < 0){
                                if (UI.simOptions.OutofAreaUseInitConc){
                                    cell_cswat = itw->second.initConc;
                                }
                                else if (UI.simOptions.OutofAreaConc < 0){
                                    continue;
                                }
                                else{
                                    cell_cswat = UI.simOptions.OutofAreaConc;
                                }
                            }
                            else{
                                int shift = UZ.traveltime(IJ);
                                if (shift > iyr){
                                    cell_cswat = itw->second.initConc;
                                }
                                else{
                                    int hru = hru_raster.getHRU(IJ);
                                    hruidx = swat.hru_index(hru);
                                    if (hruidx < 0){
                                        if (UI.simOptions.OutofAreaUseInitConc){
                                            cell_cswat = itw->second.initConc;
                                        }
                                        else if (UI.simOptions.OutofAreaConc < 0){
                                            continue;
                                        }
                                        else{
                                            cell_cswat = UI.simOptions.OutofAreaConc;
                                        }
                                    }
                                    else{
                                        if (swat.perc_mm[iswat][hruidx] < 0.01){
                                            // How to modify concentration if perc is zero?
                                            itw->second.strml[i].gw_conc.push_back(0.0);
                                            cell_cswat = swat.Salt_perc_ppm[iswat][hruidx];
                                        }
                                        else{
                                            if (UI.swatOptions.version == 1){// Version 1 of Feedback loop
                                                if (UI.simOptions.EnableFeedback == false) {
                                                    cell_cgw = 0.0;
                                                    cell_cswat = swat.Salt_perc_ppm[iswat][hruidx];
                                                }
                                                else{
                                                    double m_total = 0;
                                                    // What happens if irrGW_mm is 0. the m_npsat becomes zero.
                                                    cell_cgw = ConcFromPump[IJ];
                                                    double m_npsat = ConcFromPump[IJ] * swat.irrGW_mm[iswat][hruidx] / 100.0;
                                                    if (m_npsat < 0.001){
                                                        cell_cgw = 0.0;
                                                    }

                                                    m_total = swat.irrsaltSW_kgha[iswat][hruidx] +
                                                              m_npsat +
                                                              swat.fertsalt_kgha[iswat][hruidx] +
                                                              swat.dssl_kgha[iswat][hruidx] -
                                                              swat.Qsalt_kgha[iswat][hruidx] -
                                                              swat.uptk_kgha[iswat][hruidx] -
                                                              swat.dSoilSalt_kgha[iswat][hruidx];
                                                    m_total = m_total * swat.pGW[iswat][hruidx];
                                                    if (m_total < 0){
                                                        m_total = 0;
                                                        //bool check_this = true;
                                                        //if (m_total < min_neg_m){
                                                        //    min_neg_m = m_total;
                                                        //}
                                                    }


                                                    //cell_mfeed = m_npsat - swat.irrsaltGW_Kgha[iswat][hruidx];
                                                    //if (cell_mfeed > 0){
                                                    //    cell_mfeed = cell_mfeed * swat.pGW[iswat][hruidx];
                                                    //}

                                                    //if (UI.EnableFeedback) {
                                                    //    cell_mfeed = 0.0;
                                                    //}

                                                    //m_total = swat.totpercsalt_kgha[iswat][hruidx] + cell_mfeed;
                                                    //if (m_total < 0){
                                                    // This can happend when totpercsalt_kgha < irrsaltGW_Kgha
                                                    //    m_total = 0;
                                                    //}

                                                    cell_cswat = m_total * 100 / swat.perc_mm[iswat][hruidx];
                                                }


                                            }
                                            else if (UI.swatOptions.version == 0){
                                                // Calculate the ratio of groundwater to the total
                                                double gw_ratio = 0.0;
                                                if (std::abs(swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]) > 0.00000001){
                                                    gw_ratio = swat.irrGW_mm[iswat][hruidx] / (swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]);
                                                }

                                                // irrigated Salt Mass. The Salts that come from the pumped water
                                                double m_gw = swat.irrsaltGW_Kgha[iswat][hruidx];
                                                //Calculate the percolated groundwater volume
                                                double v_gw = swat.perc_mm[iswat][hruidx] * gw_ratio;

                                                // Calculate concentration from NPSAT spreading (Feedback)
                                                double m_npsat = ConcFromPump[IJ] * v_gw / 100.0;

                                                cell_mfeed = m_npsat - m_gw;
                                                if (cell_mfeed < 0) {
                                                    cell_mfeed = 0.0;
                                                }
                                                if (UI.simOptions.EnableFeedback == false) {
                                                    cell_mfeed = 0.0;
                                                }

                                                cell_cswat = (swat.totpercsalt_kgha[iswat][hruidx] + cell_mfeed) * 100 / swat.perc_mm[iswat][hruidx];
                                            }
                                        }
                                    }
                                    if (UI.simOptions.SurfConcValue > 0){
                                        double surf_perc = UZ.getSurfPerc(IJ);
                                        if (surf_perc > 0){
                                            cell_cswat = cell_cswat*(1-surf_perc) + surf_perc*UI.simOptions.SurfConcValue;
                                        }
                                    }
                                }
                            }

                            c_swat = c_swat + cell_cswat;
                            Mfeed = Mfeed + cell_cgw;
                            n_cswat = n_cswat + 1.0;
                            addThis = true;
                        }// loop jj
                    }// loop ii
                    c_swat = c_swat / n_cswat;
                    Mfeed = Mfeed / n_cswat;

                    if (c_swat > UI.simOptions.MaxConc){
                        c_swat = UI.simOptions.MaxConc;
                    }
                }
                if (!addThis){
                    continue;
                }

                */

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
        for (itw = VD.begin(); itw != VD.end(); ++itw){
            //std::cout << itw->first << std::endl;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                double c_swat = 0.0;
                double Mfeed = 0.0;
                int hruidx = -9;
                bool addThis = false;

                if (itw->second.strml[i].inRiv){
                    c_swat = UI.riverOptions.ConcValue;
                    addThis = true;
                }
                else if(itw->second.strml[i].a > UI.simOptions.MaxAge){
                    c_swat = itw->second.initConc;
                    addThis = true;
                }
                else {
                    double n_cswat = 0.0;
                    for (int ii = itw->second.strml[i].urfI - UI.simOptions.nBuffer; ii <= itw->second.strml[i].urfI + UI.simOptions.nBuffer; ++ii) {
                        for (int jj = itw->second.strml[i].urfJ - UI.simOptions.nBuffer; jj <= itw->second.strml[i].urfJ + UI.simOptions.nBuffer; ++jj) {
                            double cell_cswat = 0.0;
                            double cell_mfeed = 0.0;
                            double cell_cgw = 0.0;
                            int IJ = backRaster.IJ(ii, jj);
                            if (IJ < 0) {
                                if (UI.simOptions.OutofAreaUseInitConc){
                                    cell_cswat = itw->second.initConc;
                                }
                                else if (UI.simOptions.OutofAreaConc < 0){
                                    continue;
                                }
                                else{
                                    cell_cswat = UI.simOptions.OutofAreaConc;
                                }
                            }
                            else{
                                int shift = UZ.traveltime(IJ);
                                if (shift > iyr) {
                                    cell_cswat = itw->second.initConc;
                                }
                                else {
                                    int hru = hru_raster.getHRU(IJ);
                                    hruidx = swat.hru_index(hru);
                                    if (hruidx < 0) {
                                        if (UI.simOptions.OutofAreaUseInitConc){
                                            cell_cswat = itw->second.initConc;
                                        }
                                        else if (UI.simOptions.OutofAreaConc < 0) {
                                            continue;
                                        }
                                        else {
                                            cell_cswat = UI.simOptions.OutofAreaConc;
                                        }
                                    } else {
                                        if (swat.perc_mm[iswat][hruidx] < 0.01) {
                                            itw->second.strml[i].gw_conc.push_back(0.0);
                                            cell_cswat = swat.Salt_perc_ppm[iswat][hruidx];
                                        }
                                        else {
                                            if (UI.swatOptions.version == 1) {// Version 1 of Feedback loop
                                                if (UI.simOptions.EnableFeedback == false) {
                                                    cell_cgw = 0.0;
                                                    cell_cswat = swat.Salt_perc_ppm[iswat][hruidx];
                                                }
                                                else{
                                                    double m_total = 0;
                                                    // What happens if irrGW_mm is 0. the m_npsat becomes zero.
                                                    cell_cgw = ConcFromPump[IJ];
                                                    double m_npsat = ConcFromPump[IJ] * swat.irrGW_mm[iswat][hruidx] / 100.0;
                                                    if (m_npsat < 0.001){
                                                        cell_cgw = 0.0;
                                                    }

                                                    m_total = swat.irrsaltSW_kgha[iswat][hruidx] +
                                                              m_npsat +
                                                              swat.fertsalt_kgha[iswat][hruidx] +
                                                              swat.dssl_kgha[iswat][hruidx] -
                                                              swat.Qsalt_kgha[iswat][hruidx] -
                                                              swat.uptk_kgha[iswat][hruidx] -
                                                              swat.dSoilSalt_kgha[iswat][hruidx];
                                                    m_total = m_total * swat.pGW[iswat][hruidx];
                                                    if (m_total < 0){
                                                        m_total = 0;
                                                        //bool check_this = true;
                                                        //if (m_total < min_neg_m){
                                                        //    min_neg_m = m_total;
                                                        //}
                                                    }

                                                    //cell_mfeed = m_npsat - swat.irrsaltGW_Kgha[iswat][hruidx];
                                                    //if (cell_mfeed > 0){
                                                    //    cell_mfeed = cell_mfeed * swat.pGW[iswat][hruidx];
                                                    //}

                                                    //if (UI.EnableFeedback) {
                                                    //    cell_mfeed = 0.0;
                                                    //}

                                                    //m_total = swat.totpercsalt_kgha[iswat][hruidx] + cell_mfeed;
                                                    //if (m_total < 0){
                                                    // This can happend when totpercsalt_kgha < irrsaltGW_Kgha
                                                    //    m_total = 0;
                                                    //}

                                                    cell_cswat = m_total * 100 / swat.perc_mm[iswat][hruidx];
                                                }
                                            }
                                            else if (UI.swatOptions.version == 0){
                                                // Calculate the ratio of groundwater to the total
                                                double gw_ratio = 0.0;
                                                if (std::abs(swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]) >
                                                    0.00000001) {
                                                    gw_ratio = swat.irrGW_mm[iswat][hruidx] /
                                                               (swat.irrGW_mm[iswat][hruidx] + swat.irrSW_mm[iswat][hruidx]);
                                                }

                                                // irrigated Salt Mass. The Salts that come from the pumped water
                                                double m_gw = swat.irrsaltGW_Kgha[iswat][hruidx];

                                                //Calculate the percolated groundwater volume
                                                double v_gw = swat.perc_mm[iswat][hruidx] * gw_ratio;

                                                // Calculate concentration from NPSAT spreading (Feedback)
                                                double m_npsat = ConcFromPump[IJ] * v_gw / 100.0;
                                                cell_mfeed = m_npsat - m_gw;
                                                if (cell_mfeed < 0) {
                                                    cell_mfeed = 0.0;
                                                }
                                                if (UI.simOptions.EnableFeedback == false) {
                                                    cell_mfeed = 0.0;
                                                }
                                                cell_cswat = (swat.totpercsalt_kgha[iswat][hruidx] + cell_mfeed) * 100 /
                                                             swat.perc_mm[iswat][hruidx];
                                            }

                                        }
                                    }

                                    if (UI.npsatOptions.version == 1){
                                        double u = MS::calcRiverInfluence(itw->second.strml[i].rivRist, UI.riverOptions);
                                        cell_cswat = u * cell_cswat + (1-u) * UI.riverOptions.ConcValue;
                                    }

                                    if (UI.simOptions.SurfConcValue > 0) {
                                        double surf_perc = UZ.getSurfPerc(IJ);
                                        if (surf_perc > 0) {
                                            cell_cswat = cell_cswat * (1 - surf_perc) + surf_perc * UI.simOptions.SurfConcValue;
                                        }
                                    }
                                }
                            }

                            c_swat = c_swat + cell_cswat;
                            Mfeed = Mfeed + cell_cgw;
                            n_cswat = n_cswat + 1.0;
                            addThis = true;
                        }// loop jj
                    }// loop ii

                    c_swat = c_swat / n_cswat;
                    Mfeed = Mfeed / n_cswat;

                    if (c_swat > UI.simOptions.MaxConc) {
                        c_swat = UI.simOptions.MaxConc;
                    }
                }
                if (addThis){
                    itw->second.strml[i].gw_conc.push_back(Mfeed);
                    itw->second.strml[i].lf_conc.push_back(c_swat);
                }
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
            for (itw = VD.begin(); itw != VD.end(); ++itw){
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
        if (iswat >= UI.simOptions.Nyears){
            iswat = 0;
        }
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
