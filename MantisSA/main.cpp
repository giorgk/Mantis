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
#include "MS_well_cell.h"
#include "MSdebug.h"
#include "MS_tests.h"



int main(int argc, char* argv[]) {
    boost::mpi::environment env( argc, argv );
    boost::mpi::communicator world;


    //MS::TESTS::sendVec2Root(world);
    //MS::TESTS::sendScalarFromRoot2AllProc(world);
    //MS::TESTS::sendVectorFromRoot2AllProc(world);
    //MS::TESTS::sentMatrixFromRoot2AllProc(world, false);
    //MS::TESTS::sentMatrixFromRoot2AllProc(world, true);


    //{
    //    bool tf = false;
    //    std::cout << static_cast<int>(tf) << std::endl;
    //    tf = true;
    //    std::cout << static_cast<int>(tf) << std::endl;
//        MS::testConvolution();
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
    MS::yearMap(swatYRIDmap, UI.simOptions.StartYear, UI.simOptions.StartYear + UI.simOptions.Nyears, UI.swatOptions.StartYear, UI.swatOptions.Nyears);09

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
        tf = MS::readSelectedWells(UI.outputOptions.SelectedWells, UI.outputOptions.SelectedWellsGroups, SWGmap, world);
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

    world.barrier();
    if (world.rank() == 0){
        std::cout << "The unsaturated data are using the old matlab format and should be updated if the depth data are updated to the new format" << std::endl;
    }
    world.barrier();
    MS::UNSAT UZ;
    tf = UZ.readdata(UI.unsatOptions.Depth_file, UI.unsatOptions.Depth_name, UI.unsatOptions.Rch_file,
                     UI.unsatOptions.wc, UI.unsatOptions.minDepth,UI.unsatOptions.minRch, world);
    if (!tf){ return 0;}

    MS::HRU_Raster hru_raster;
    tf = hru_raster.read(UI.swatOptions.HRU_raster_file, world, backRaster.Ncell());
    if (!tf){
        return 0;
    }

    MS::SWAT_data swat(UI.swatOptions.nhrus, UI.swatOptions.Nyears);
    tf = swat.read_HRU_idx_Map(UI.swatOptions.HRU_index_file, world);
    if (!tf){
        return 0;
    }
    tf = swat.read(UI.swatOptions.Data_file, UI.swatOptions.version, world);
    if (!tf){
        return 0;
    }

    MS::HistoricLoading HISTLOAD;
    if (!UI.historicOptions.filename.empty()){
        HISTLOAD.n_reserve = backRaster.Ncell();
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
    if (!well_cells.read_build(UI.npsatOptions.DistribuPumpFile, world)) {
        return 0;
    }
    world.barrier();

    //std::vector<std::vector<int>> Eid_proc;
    //sendWellEids(VI, Eid_proc, world);

    if (UI.npsatOptions.bUseInitConcVI) {
        tf = MS::readInitSaltConc(UI.npsatOptions.InitSaltVIFile, VI, world);
        if (!tf){return 0;}
    }
    if (UI.npsatOptions.bUseInitConcVD) {
        tf = MS::readInitSaltConc(UI.npsatOptions.InitSaltVDFile, VD, world);
        if (!tf){return 0;}
    }

    world.barrier();
    std::cout << "Proc: " << world.rank() << " has " << VI.wells.size() << " VI wells and " << VD.wells.size() << " VD wells" << std::endl;


    std::vector<double> ConcFromPump(backRaster.Ncell(), 0);
    std::vector<int> WellEidFromPump(backRaster.Ncell(), 0);
    MS::Matrix<double> mass_removed;
    MS::Matrix<double> volume_removed;
    if (world.rank() == 0) {
        mass_removed.allocate(UI.swatOptions.nhrus, UI.simOptions.Nyears);
        volume_removed.allocate(UI.swatOptions.nhrus, UI.simOptions.Nyears);
    }
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
        std::vector<std::vector<double>> AllwellsInitConc;
        MS::sendVec2Root<double>(wellInitConc, AllwellsInitConc, world);
        // Put the well concentrations in the right cells
        if (world.rank() == 0){
            //MS::WELL_CELLS::iterator it; // well_cells
            for (int i = 0; i < world.size(); ++i){
                for (int j = 0; j < static_cast<int>(AllwellsInitConc[i].size()); j = j + 2){
                    const int eid = static_cast<int>(AllwellsInitConc[i][j]);
                    const double conc = AllwellsInitConc[i][j + 1];

                    const int* begin = nullptr;
                    const int* end   = nullptr;

                    if (well_cells.get_cells_ptr(eid, begin, end)) {
                        for (const int* p = begin; p != end; ++p) {
                            ConcFromPump[*p] = conc;
                            WellEidFromPump[*p] = eid;
                        }
                    }

                    // it = well_cells.find(eid);
                    // if (it != well_cells.end()){
                    //     for (int k = 0; k < it->second.size(); ++k){
                    //         ConcFromPump[it->second[k]] = AllwellsInitConc[i][j+1];
                    //         WellEidFromPump[it->second[k]] = eid;
                    //     }
                    // }
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

    //return 0;


    // Main simulation loop

    //double min_neg_m = 0;

    //int hruidx;
    //std::vector<double> totMfeed(UI.NsimYears, 0);
    if (world.rank() == 0) {
        std::cout << std::endl;
        std::cout << "==============================================================" << std::endl;
        std::cout << "Start simulation..." << std::endl;
        std::cout << "==============================================================" << std::endl;
        std::cout << std::endl;
    }

    auto startTotal = std::chrono::high_resolution_clock::now();

    int YYYY = UI.simOptions.StartYear;
    int iswat = swatYRIDmap.find(YYYY)->second;

    for (int iyr = 0; iyr < UI.simOptions.Nyears; ++iyr){
        if (world.rank() == 0) {
            std::cout << "Year " << iyr << " ..."  << std::flush;
        }
        //std::cout << "Here1" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        const bool is_last_year = (iyr == UI.simOptions.Nyears - 1);

        // Packed as [well_id, concentration, well_id, concentration, ...]
        // Used only for feedback to update ConcFromPump on root.
        std::vector<double> wellConc;

        world.barrier();

        bool printThis = false;

        // ============================================================
        // 1. Simulate loading / BTC for VI wells for this year
        // ============================================================
        for (itw = VI.wells.begin(); itw != VI.wells.end(); ++itw){
            if (UI.doDebug){
                if (UI.dbg_id == itw->first){
                    printThis = true;
                }
            }

            auto &well = itw->second;

            // For the current year, wellbtc has size (iyr + 1), because the
            // convolution is based on lf_conc accumulated up to this year.
            std::vector<double> wellbtc(iyr + 1, 0.0);

            double sumW = 0;
            int cntS = 0;

            for (unsigned int i = 0; i < well.strml.size(); ++i) {
                auto &strm = well.strml[i];

                // Compute loading concentration for this streamline at the current year.
                double c_tot = 0.0;
                double c_gw = 0.0;

                MS::CreateLoadingForStreamline(c_tot, c_gw, iyr, iswat, YYYY, VI.version,
                                               strm, UI, well.initConc,
                                               HISTLOAD, backRaster, UZ, hru_raster, swat,
                                               ConcFromPump);

                // Store the yearly loading history for this streamline.
                strm.gw_conc.push_back(c_gw);
                strm.lf_conc.push_back(c_tot);

                // Convolve the streamline loading history with the streamline URF.
                std::vector<double> btc(strm.lf_conc.size(), 0);
                std::vector<double> prebtc(strm.lf_conc.size(), 0);

                MS::convolute(strm.urf, strm.lf_conc, btc, prebtc, well.initConc);

                // Add weighted streamline contribution to the well BTC.
                for (unsigned int j = 0; j < btc.size(); ++j){
                    const double c_here = btc[j] + prebtc[j];
                    wellbtc[j] += strm.W * c_here;

                    // On the last simulation year, store the full BTC per streamline.
                    if (is_last_year){
                        strm.btc.push_back(c_here);
                    }
                }

                sumW += strm.W;
                cntS += 1;
                //well.m_rmv[i] += m_rmv;

                //if (printThis){
                //    dbg_file << iyr << ", " << itw->first << ", " << itw->second.strml[i].Sid << ", " << hruidx << ", "
                //            //<< gw_ratio << ", " << m_gw << ", " << v_gw << ", " << m_npsat << ", " << gw_conc << ", "
                //            << c_swat << ", " << shift << std::endl;
                //}
            }// Loop streamlines

            // No active streamlines for this well
            if (cntS == 0){
                if (is_last_year && UI.simOptions.OutofAreaConc < 0) {
                    // if we ignore out of area streamlines we have to ignore the wells
                    // with streamlines that originate from outside. To avoid breaking the printing procedures
                    // we will print -9
                    for (unsigned int i = 0; i < wellbtc.size(); ++i) {
                        wellbtc[i] = -9.0;
                    }
                }
                sumW = 1.0;
            }

            // On the last year, store the full well BTC.
            // Otherwise store only the current-year value for feedback.
            if (is_last_year) {
                for (unsigned int i = 0; i < wellbtc.size(); ++i) {
                    well.wellBtc.push_back(wellbtc[i] / sumW);
                }
            }
            else{
                wellConc.push_back(static_cast<double>(itw->first));
                wellConc.push_back(wellbtc.back() / sumW);
            }
            if (printThis){printThis = false;}
        }// end loop VI wells

        // ============================================================
        // 2. For VD wells, only build yearly loading histories here.
        //    Full BTC is computed only on the last year.
        // ============================================================
        for (itw = VD.wells.begin(); itw != VD.wells.end(); ++itw) {
            auto &well = itw->second;

            for (unsigned int i = 0; i < well.strml.size(); ++i) {
                auto &strm = well.strml[i];

                double c_tot = 0.0;
                double c_gw = 0.0;
                //double m_rmv = 0.0;

                MS::CreateLoadingForStreamline(c_tot, c_gw, iyr, iswat, YYYY, VD.version,
                                               strm, UI, well.initConc,
                                               HISTLOAD, backRaster, UZ, hru_raster, swat,
                                               ConcFromPump);

                strm.gw_conc.push_back(c_gw);
                strm.lf_conc.push_back(c_tot);
            }
        }

        // ============================================================
        // 3. Feedback step for VI wells (not needed on last year)
        // ============================================================
        if (!is_last_year) {
            if (UI.simOptions.EnableFeedback) {
                // Gather current-year well concentrations to root.
                std::vector<std::vector<double>> AllwellsConc(world.size());
                MS::sendVec2Root<double>(wellConc, AllwellsConc, world);

                // Root updates ConcFromPump in all cells associated with each well.
                if (world.rank() == 0) {
                    for (int i = 0; i < world.size(); ++i){
                        for (int j = 0; j < static_cast<int>(AllwellsConc[i].size()); j = j + 2){
                            const int eid = static_cast<int>(AllwellsConc[i][j]);
                            const double conc = AllwellsConc[i][j + 1];

                            const int* begin = nullptr;
                            const int* end   = nullptr;

                            if (well_cells.get_cells_ptr(eid, begin, end)) {
                                for (const int* p = begin; p != end; ++p) {
                                    ConcFromPump[*p] = conc;
                                }
                            }
                        }
                    }

                    // Salt removal
                    if (UI.saltRemoveOptions.enable) {
                        for (int i = 0; i < ConcFromPump.size(); ++i) {
                            constexpr double cell_area = 2500.0;
                            double mass_remove_cell = 0.0;
                            const int hru = hru_raster.getHRU(i);// This is the hru id for the cell i
                            if (hru < 0) {
                                continue;
                            }
                            const int hru_idx = swat.hru_index(hru); // This is the row index of the hru in the swat data
                            if (hru_idx < 0) {
                                continue;
                            }

                            if (hru_idx >= mass_removed.num_rows() || iyr >= mass_removed.num_cols()) {
                                std::cout << "mass_removed index error: hru_idx=" << hru_idx
                                          << " iyr=" << iyr
                                          << " rows=" << mass_removed.num_rows()
                                          << " cols=" << mass_removed.num_cols() << std::endl;
                                __debugbreak();
                            }

                            if (hru_idx >= swat.irrSW_mm.num_rows() || iswat >= swat.irrSW_mm.num_cols()) {
                                std::cout << "irrSW_mm index error: hru_idx=" << hru_idx
                                          << " iswat=" << iswat << std::endl;
                                __debugbreak();
                            }

                            double volume_SW_cell = swat.irrSW_mm(hru_idx,iswat) * cell_area / 1000.0;
                            if (volume_SW_cell > 0.01) {
                                bool tf1 = true;
                            }
                            double volume_GW_cell = swat.irrGW_mm(hru_idx,iswat) * cell_area / 1000.0;
                            if (volume_GW_cell > 0.01) {
                                bool tf2 = true;
                            }
                            double mass_SW_cell = swat.irrsaltSW_kgha(hru_idx,iswat) * cell_area/10000.0;
                            double mass_GW_cell = 0.001 * ConcFromPump[i] * volume_GW_cell;
                            double mass_cell = mass_GW_cell + mass_SW_cell;
                            double volume_cell = volume_GW_cell + volume_SW_cell;
                            double conc_trgt = UI.saltRemoveOptions.Trgt_AW_ppm;
                            if (UI.saltRemoveOptions.usefield) {
                                conc_trgt = swat.Trgt_AW_ppm(hru_idx,iswat);
                            }
                            double mass_trgt_cell = 0.001 * conc_trgt * volume_cell;
                            if (mass_trgt_cell > mass_cell) {
                                bool tf3 = true;
                            }
                            mass_remove_cell = std::max(0.0, mass_cell - mass_trgt_cell);
                            if (mass_remove_cell > 0.01) {
                                bool tf4 = true;
                            }
                            if (volume_GW_cell < 0.01) {
                                ConcFromPump[i] = 0.0;
                            }
                            else {
                                const double mass_GW_new = std::max(0.0, mass_cell - mass_remove_cell - mass_SW_cell);
                                ConcFromPump[i] = 1000.0 * mass_GW_new / volume_GW_cell;
                            }
                            mass_removed(hru_idx,iyr) += mass_remove_cell;
                            volume_removed(hru_idx,iyr) += volume_cell;
                        }
                        //std::cout << "Mass removed: " << mass_removed.sum_column(iyr)/1000000.0 << " kt"  << std::endl;

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
                std::cout << " in " << elapsed.count() << " sec" << "(Mass removed: " << mass_removed.sum_column(iyr)/1000000.0 << " kt)"  << std::endl;
            }
        }
        else{
            // ========================================================
            // 4. Last year: compute full BTCs for VD wells
            // ========================================================
            world.barrier();

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;

            if (world.rank() == 0){
                std::cout << "Year: " << iyr << " in " << elapsed.count() << std::endl;
                std::cout << "Simulating VD BTCs ..." << std::endl;
            }

            start = std::chrono::high_resolution_clock::now();
            world.barrier();

            for (itw = VD.wells.begin(); itw != VD.wells.end(); ++itw) {
                auto &well = itw->second;

                std::vector<double> wellbtc(UI.simOptions.Nyears, 0.0);
                double sumW = 0;
                int cntS = 0;

                for (unsigned int i = 0; i < well.strml.size(); ++i) {
                    auto &strm = well.strml[i];

                    if (strm.lf_conc.empty()) {
                        continue;
                    }

                    std::vector<double> btc(strm.lf_conc.size(), 0.0);
                    std::vector<double> prebtc(strm.lf_conc.size(), 0.0);

                    MS::convolute(strm.urf, strm.lf_conc, btc, prebtc, well.initConc);

                    for (unsigned int j = 0; j < btc.size(); ++j){
                        const double c_here = btc[j] + prebtc[j];
                        wellbtc[j] += strm.W * c_here;
                        strm.btc.push_back(c_here);
                    }

                    sumW += strm.W;
                    cntS += 1;
                }

                if (cntS == 0){
                    if (UI.simOptions.OutofAreaConc < 0) {
                        // Mark invalid / ignored wells with -9 so downstream printing logic works.
                        for (unsigned int i = 0; i < wellbtc.size(); ++i){
                            wellbtc[i] = -9.0;
                        }
                    }
                    sumW = 1.0;
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

        ++iswat;
        if (iswat >= UI.swatOptions.Nyears){
            iswat = 0;
        }

        // Advance calendar year
        ++YYYY;
    }

    auto finishTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTotal = finishTotal - startTotal;
    if (world.rank() == 0){
        std::cout << std::endl;
        std::cout << "==============================================================" << std::endl;
        std::cout << "Total simulation for " << UI.simOptions.Nyears << " : " << elapsedTotal.count() / 60.0 << " min" << std::endl;
        std::cout << "==============================================================" << std::endl;
        std::cout << std::endl;
        std::cout << "Printing simulation results..." << std::endl;
    }

    {// Print mass removal
        if (UI.saltRemoveOptions.enable) {
            if (world.rank() == 0) {
                {
                    std::cout << "Printing mass removed ..." << std::endl;
                    std::string fn;
                    if (UI.outputOptions.compress){
                        fn = UI.outputOptions.OutFile + "_Mrmv.dat.gz";
                    }
                    else{
                        fn = UI.outputOptions.OutFile + "_Mrmv.dat";
                    }
                    MS::printHRUMassVolRemoved(swat.hru_idx_map,mass_removed,fn, "mrmv", UI.outputOptions.compress);
                }
                {
                    std::cout << "Printing Volume removed ..." << std::endl;
                    std::string fn;
                    if (UI.outputOptions.compress){
                        fn = UI.outputOptions.OutFile + "_Vrmv.dat.gz";
                    }
                    else{
                        fn = UI.outputOptions.OutFile + "_Vrmv.dat";
                    }
                    MS::printHRUMassVolRemoved(swat.hru_idx_map, volume_removed,fn, "vrmv", UI.outputOptions.compress);
                }
            }
        }
    }

    { // Print the VI
        world.barrier();
        // Flatten this rank's well BTC data into a single vector:
        // [Nwells, eid1, btc1_1, ..., btc1_Nyears, eid2, ...]
        std::vector<double> thisProcBTCVI;
        MS::linearizeWellBTCs(VI, thisProcBTCVI);

        // Gather all ranks' flattened BTC vectors to root
        std::vector<std::vector<double>> AllProcBTCVI;
        MS::sendVec2Root<double>(thisProcBTCVI, AllProcBTCVI, world);

        if (world.rank() == 0){
            std::cout << "Printing VI BTCs ..." << std::endl;
            std::string fn;
            if (UI.outputOptions.compress){
                fn = UI.outputOptions.OutFile + "_VI.dat.gz";
            }
            else{
                fn = UI.outputOptions.OutFile + "_VI.dat";
            }
            MS::printWELLSfromAllProc(AllProcBTCVI, fn, UI.simOptions.Nyears, UI.outputOptions.compress);
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
            std::string fn = UI.outputOptions.OutFile + "_VD.dat";
            if (UI.outputOptions.compress){
                fn = fn + ".gz";
            }
            MS::printWELLSfromAllProc(AllProcBTCVD, fn, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string lf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_RSL.dat";
                    std::string mf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_RSF.dat";
                    if (UI.outputOptions.compress){
                        lf_file_name += ".gz";
                        mf_file_name += ".gz";
                    }
                    MS::printDetailOutputFromAllProc(AllProcDATA, lf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
                    MS::printMfeedFromAllProc(AllProcMfeed, mf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string urf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_RSU.dat";
                    if (UI.outputOptions.compress){
                        urf_file_name += ".gz";
                    }
                    MS::printURFsFromAllProc(AllProcDATA, urf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string btc_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VI_RSC.dat";
                    if (UI.outputOptions.compress){
                        btc_file_name += ".gz";
                    }
                    MS::printBTCssFromAllProc(AllProcDATA, btc_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string lf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_RSL.dat";
                    std::string mf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_RSF.dat";
                    if (UI.outputOptions.compress){
                        lf_file_name += ".gz";
                        mf_file_name += ".gz";
                    }
                    MS::printDetailOutputFromAllProc(AllProcDATA,lf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
                    MS::printMfeedFromAllProc(AllProcMfeed, mf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string urf_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_RSU.dat";
                    if (UI.outputOptions.compress){
                        urf_file_name += ".gz";
                    }
                    MS::printURFsFromAllProc(AllProcDATA, urf_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
                    std::string btc_file_name = UI.outputOptions.OutFile + "_" + it->second.groupName + "_VD_RSC.dat";
                    if (UI.outputOptions.compress){
                        btc_file_name += ".gz";
                    }
                    MS::printBTCssFromAllProc(AllProcDATA, btc_file_name, UI.simOptions.Nyears, UI.outputOptions.compress);
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
        std::cout << "==============================================================" << std::endl;
        std::cout << "MantisSA Finished at " << std::asctime(std::localtime(&result)) << std::endl;
        std::cout << "Total simulation for : " << elapsedTotalSimulation.count() / 60.0 << " min" << std::endl;
        std::cout << "==============================================================" << std::endl;
    }

    return 0;
}
