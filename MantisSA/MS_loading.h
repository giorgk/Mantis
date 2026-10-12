//
// Created by giorg on 5/30/2025.
//

#ifndef MANTISSA_MS_LOADING_H
#define MANTISSA_MS_LOADING_H

#include <random>
#include "MSdebug.h"

namespace MS{
    class HistoricLoading{
    public:
        HistoricLoading(){}
        //int YearA;
        //int YearB;
        //int idx;
        int n_reserve = 0;
        std::string Prefix;
        std::string Ext;
        std::vector<std::vector<double>> HistLoad;
        std::vector<int> Years;
        bool Setup(std::string pref, std::string ext, boost::mpi::communicator &world);
        bool ReadLoading(boost::mpi::communicator &world);
        double calculateConc(int IJ, int year) const;

        bool readInfo(std::string filename, boost::mpi::communicator &world);
        bool readRasterLoad(const std::string &filename, std::vector<double> &V,boost::mpi::communicator &world);
        bool RootReadsHDF_HistFileDistrib(const std::string &filename,
                                         std::vector<double> &V,
                                         boost::mpi::communicator &world);
    };

    bool HistoricLoading::ReadLoading(boost::mpi::communicator &world) {
        for (unsigned int i = 0; i < Years.size(); ++i ){
            std::string fn = Prefix + std::to_string(Years[i]) + "." + Ext;
            std::vector<double> tmp;
            bool tf = readRasterLoad(fn, tmp, world);
            if (!tf){
                return tf;
            }
            HistLoad.push_back(tmp);
        }
        return true;
    }

    double HistoricLoading::calculateConc(int IJ, int year) const {

        if (year <= Years[0]) {
            return HistLoad[0][IJ];
        }

        const std::size_t n = Years.size();
        if (year >= Years[n - 1]) {
            return HistLoad[n - 1][IJ];
        }



        for (std::size_t i = 0; i < n - 1; ++i){
            if (year >= Years[i] && year <= Years[i+1]){
                double u = (static_cast<double>(year) - static_cast<double>(Years[i]))/
                           (static_cast<double>(Years[i+1]) - static_cast<double>(Years[i]));
                 return (1-u) * HistLoad[i][IJ] + u * HistLoad[i+1][IJ];
            }
        }

        // Should never happen if Years is sorted and covers the cases above.
        return HistLoad[n - 1][IJ];
    }

    bool HistoricLoading::readRasterLoad(const std::string &filename, std::vector<double> &V, boost::mpi::communicator &world) {
        std::string ext = getExtension(filename);
        std::vector<double> tmp;
        bool tf = false;
        if (world.rank() == 0) {
            std::cout << "Reading " << filename << std::endl;
        }
         if (ext == "h5" || ext == "H5"){
             tf = RootReadsHDF_HistFileDistrib(filename, tmp, world);
        }
        else{
            tf = RootReadsVectorFileDistrib<double>(filename, tmp,  world, 5000000, n_reserve);
        }

        if (!tf) {
            V.clear();
            return false;
        }

        if (PrintMatrices) {
            printVectorForAllProc(tmp, world, 345, 353);
        }
        if (tmp.empty()) {
            std::cout << "Error: no data read from " << filename << std::endl;
            V.clear();
            return false;
        }
        V = tmp;
        return true;

    }

    bool HistoricLoading::readInfo(std::string filename, boost::mpi::communicator &world) {
        bool tf = RootReadsVectorFileDistrib<int>(filename, Years,  world);
        if (!tf){ return false;}
        return true;
    }

    bool HistoricLoading::Setup(std::string pref, std::string ext,
                                boost::mpi::communicator &world) {
        Prefix = pref;
        Ext = ext;
        std::string finfo = Prefix + "info.dat";
        bool tf = readInfo(finfo,world);
        if (!tf){ return false;}
        tf = ReadLoading(world);
        return tf;
    }

    bool HistoricLoading::RootReadsHDF_HistFileDistrib(const std::string &filename,
                                                       std::vector<double> &V,
                                                       boost::mpi::communicator &world) {
#if _USEHF > 0
        int readSuccess = 1;
        if (world.rank() == 0) {
            try {
                const std::string NameSet("HIST");
                HighFive::File hdfFile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet dataset = hdfFile.getDataSet(NameSet);

                dataset.read(V);

                if (V.empty()) {
                    std::cout << "Error: empty dataset in " << filename << std::endl;
                    readSuccess = 0;
                }
            }
            catch (const std::exception &e) {
                std::cout << "Error reading HDF5 historic load file "
                          << filename << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
            catch (...) {
                std::cout << "Unknown error reading HDF5 historic load file "
                          << filename << std::endl;
                readSuccess = 0;
            }
        }

        sendScalarFromRoot2AllProc(readSuccess, world);

        if (readSuccess == 0) {
            V.clear();
            return false;
        }

        sendVectorFromRoot2AllProc(V, world);
        return true;

#else
        (void)filename;
        (void)V;
        (void)world;
        std::cout << "Error: HDF5 support is not enabled." << std::endl;
        return false;
#endif
    }


    // Computes the future-loading concentration for one raster cell.
    //
    // Outputs:
    //   c_cell     = final concentration assigned to this cell for the future-loading regime
    //   gw_cell    = groundwater-irrigation concentration contribution diagnostic
    //   water_cell = effective percolating-water depth used to weight this cell in the
    //                streamline concentration calculation
    //   return     = true when water_cell comes from a valid SWAT HRU; false when the
    //                caller must use a fallback water-volume estimate
    //
    // Inputs:
    //   IJ            = active-cell index in the background raster
    //   iswat         = SWAT time index
    //   ver           = option flag; if ver == 1, blend with river concentration using RivInfl
    //   initConc      = fallback concentration if HRU information is unavailable
    //   concfromPump  = concentration associated with pumped groundwater used for irrigation
    //   RivInfl       = fraction of river influence [0,1]
    //   rivConcVal    = prescribed river concentration
    //   SurfConcValue = optional surface-source concentration mixed by surf_perc
    //   surf_perc     = fraction of concentration attributed to surface source [0,1]
    //   hru_raster    = maps raster cell IJ to SWAT HRU
    //   swat          = SWAT outputs for percolation, salts, uptake, etc.
    //   dp_mult       = dilution / deep percolation multiplier used in final concentration conversion
    //
    // Logic:
    //   1. Find the SWAT HRU for this cell.
    //   2. If no HRU exists, use initConc.
    //   3. If percolation is nearly zero, use the SWAT percolate salt concentration directly.
    //   4. Otherwise, form a simple salt mass balance and convert to concentration.
    //   5. Optionally blend with river concentration.
    //   6. Optionally blend with a prescribed surface concentration.
    bool calculateFutureLoading(double &c_cell, double &gw_cell, double &water_cell,
                                const int &IJ, const int &iswat,
                                const int &ver, const double &initConc, const double &concfromPump, const double &RivInfl,
                                const double &rivConcVal, const double &SurfConcValue, const double &surf_perc,
                                const MS::HRU_Raster &hru_raster,
                                const MS::SWAT_data &swat,
                                const double &dp_mult,
                                CellDebugData *debug_data
                                //const MS::MiscOptions &miscopt
                                ){
        // Default outputs
        c_cell = 0.0;
        gw_cell = 0.0;
        water_cell = 0.0;

        // Map model cell to SWAT HRU
        const int hru = hru_raster.getHRU(IJ);
        const int hruidx = swat.hru_index(hru);

        if (debug_data != nullptr) {
            debug_data->hru = hru;
            debug_data->hruidx = hruidx;
            debug_data->treated_conc_from_pump = concfromPump;
            debug_data->surf_perc = surf_perc;
        }

        // If no valid HRU exists, fall back to the initial concentration
        if (hruidx < 0){
            c_cell = initConc;
            if (debug_data != nullptr) {
                debug_data->c_before_river = c_cell;
                debug_data->c_after_river = c_cell;
                debug_data->c_after_surface = c_cell;
                debug_data->c_future = c_cell;
            }
            return false;
        }

        // CHANGE (v0.1.16): SWAT percolation is the only weight used for buffer cells.
        // All raster cells have the same area, so percolation depth is proportional to
        // water volume and the common cell area cancels from the weighted average.
        const double perc_cell = swat.perc_mm(hruidx, iswat);
        // OLD (pre-v0.1.16): water_cell = perc_cell * dp_mult;
        water_cell = perc_cell;
        if (water_cell < 0.0){
            water_cell = 0.0;
        }

        if (debug_data != nullptr) {
            debug_data->irrGW_mm = swat.irrGW_mm(hruidx, iswat);
            debug_data->irrSW_mm = swat.irrSW_mm(hruidx, iswat);
            debug_data->irrsaltSW_kgha = swat.irrsaltSW_kgha(hruidx, iswat);
            debug_data->fertsalt_kgha = swat.fertsalt_kgha(hruidx, iswat);
            debug_data->dssl_kgha = swat.dssl_kgha(hruidx, iswat);
            debug_data->Qsalt_kgha = swat.Qsalt_kgha(hruidx, iswat);
            debug_data->uptk_kgha = swat.uptk_kgha(hruidx, iswat);
            debug_data->pGW = swat.pGW(hruidx, iswat);
            debug_data->perc_mm = perc_cell;
            debug_data->Salt_perc_ppm = swat.Salt_perc_ppm(hruidx, iswat);
            debug_data->water_weight = water_cell;
        }

        // Very small percolation:
        // avoid unstable division by a tiny number and use SWAT's percolate concentration directly.
        if (perc_cell < 0.01){
            c_cell = swat.Salt_perc_ppm(hruidx, iswat);
            gw_cell = 0.0;
        }
        else{
            // Groundwater-irrigation concentration diagnostic.
            // This tracks the pumped-groundwater contribution used by the future-loading calculation.
            gw_cell = concfromPump;

            // Mass contribution from NPSAT groundwater irrigation term.
            // Unit consistency depends on the SWAT/NPSAT convention already used in the rest of the model.
            double m_npsat = concfromPump * swat.irrGW_mm(hruidx, iswat) / 100.0;

            // Small or negative groundwater contribution:
            // set the diagnostic to zero, and do not allow negative mass.
            if (m_npsat < 0.01){
                gw_cell = 0.0;
                if (m_npsat < 0.0){
                    m_npsat = 0;
                }
            }

            if (debug_data != nullptr) {
                debug_data->m_npsat = m_npsat;
            }

            // Total salt mass available to deep percolation.
            // Positive sources:
            //   - irrigation surface-water salt
            //   - groundwater-irrigation salt (m_npsat)
            //   - fertilizer salt
            //   - dissolution / other source term
            // Negative sinks:
            //   - runoff/export
            //   - plant uptake
            const double m_total_raw = swat.irrsaltSW_kgha(hruidx, iswat) +
                      m_npsat +
                      swat.fertsalt_kgha(hruidx, iswat) +
                      swat.dssl_kgha(hruidx, iswat) -
                      swat.Qsalt_kgha(hruidx, iswat) -
                      swat.uptk_kgha(hruidx, iswat);// -
                    //swat.dSoilSalt_kgha[hruidx][iswat];

            double m_total = m_total_raw;

            // Apply the groundwater partition / weighting term used by the model.
            m_total = m_total * swat.pGW(hruidx, iswat);

            // Do not allow negative total mass
            if (m_total < 0){
                m_total = 0.0;
            }

            if (debug_data != nullptr) {
                debug_data->m_total_raw = m_total_raw;
                debug_data->m_total = m_total;
            }

            // Convert mass loading to concentration.
            // perc_mm is the deep percolation denominator,
            // dp_mult is an additional scaling / dilution factor.
            //double c_cell_tmp = m_total * 100 / (swat.perc_mm[hruidx][iswat]);
            c_cell = m_total * 100 / (perc_cell * dp_mult);
        }

        if (debug_data != nullptr) {
            debug_data->c_before_river = c_cell;
        }

        // Optional river-influence blending.
        // ver == 1 activates mixing between computed cell concentration and river concentration.
        if (ver == 1){
            c_cell = (1 - RivInfl) * c_cell + RivInfl * rivConcVal;
        }

        if (debug_data != nullptr) {
            debug_data->c_after_river = c_cell;
        }

        // Optional surface-source blending.
        if (SurfConcValue > 0 && surf_perc > 0) {
            c_cell = c_cell * (1 - surf_perc) + surf_perc * SurfConcValue;
        }

        if (debug_data != nullptr) {
            debug_data->c_after_surface = c_cell;
            debug_data->c_future = c_cell;
        }

        return true;
    }

    void CreateLoadingForStreamline(double &c_strml, double &gw_strml,
                                    const int &iyr, // simulation year index (relative index, not calendar year)
                                    const int &iswat, // SWAT time index for future loading
                                    const int &yyyy, // calendar year corresponding to iyr
                                    const int &ver, // option flag controlling river influence blending
                                    const STRML &strm, // streamline data (location, age, river influence, etc.)
                                    const UserInput &UI,
                                    const double &initConc, // fallback / initial concentration
                                    const HistoricLoading &histload,
                                    const MS::BackgroundRaster &backRaster,
                                    const MS::UNSAT &UZ, // unsaturated travel time / surface fraction info
                                    const MS::HRU_Raster &hru_raster,
                                    const MS::SWAT_data &swat,
                                    const std::vector<double> &ConcFromPump,
                                    const std::vector<int> &WellEidFromPump,
                                    std::ofstream *debug_file = nullptr,
                                    int debug_eid = -1){

        // Output defaults
        c_strml = 0.0;
        gw_strml = 0.0;

        // Case 1:
        // If the streamline is in a river, use the prescribed river concentration directly.
        // No neighborhood averaging is needed.
        if (strm.inRiv){
            c_strml = UI.riverOptions.ConcValue;
            return;
        }

        // Case 2:
        // If the streamline age exceeds the model's maximum age,
        // assign the initial concentration and stop.
        if (strm.a > UI.simOptions.MaxAge){
            c_strml = initConc;
            return;
        }

        // CHANGE (v0.1.12): concentrations are now weighted by the water volume represented
        // by each cell. This is equivalent to summing mass and volume, while also preserving
        // the concentration-level historic, river, surface-water, and transition blends.
        double weighted_c_sum = 0.0;
        double weighted_gw_sum = 0.0;
        double total_water = 0.0;

        // OLD (pre-v0.1.12): arithmetic averaging gave dry and wet cells equal influence.
        // double n_cells = 0.0;

        // CHANGE (v0.1.17): zero-based counter identifying each buffer position
        // for this streamline. For nBuffer == 1 the counter runs from 0 to 8.
        int buffer_cell = 0;

        for (int ii = strm.urfI - UI.simOptions.nBuffer; ii <= strm.urfI + UI.simOptions.nBuffer; ++ii){
            for (int jj = strm.urfJ - UI.simOptions.nBuffer; jj <= strm.urfJ + UI.simOptions.nBuffer; ++jj){
                double c_cell = 0.0; // concentration assigned to this neighborhood cell
                double gw_cell = 0.0; // groundwater contribution diagnostic for this cell
                double water_cell = 0.0; // effective water depth/volume used as this cell's weight
                CellDebugData cell_debug;
                CellDebugData *debug_data = nullptr;
                if (debug_file != nullptr) {
                    debug_data = &cell_debug;
                }

                // Convert raster coordinates to active-cell index.
                int IJ = backRaster.IJ(ii, jj);

                // If the raster cell is inactive / outside domain, use initial concentration.
                if (IJ < 0){
                    c_cell = initConc;
                    cell_debug.regime = "inactive";
                    // CHANGE (v0.1.12): an inactive cell has no modeled recharge volume, so it
                    // does not bias a streamline near the edge of the active raster.
                    water_cell = 0.0;
                }
                else{
                    cell_debug.well_source_id = WellEidFromPump[IJ];
                    cell_debug.hru = hru_raster.getHRU(IJ);
                    cell_debug.hruidx = swat.hru_index(cell_debug.hru);
                    cell_debug.treated_conc_from_pump = ConcFromPump[IJ];
                    cell_debug.surf_perc = UZ.getSurfPerc(IJ);

                    // CHANGE (v0.1.16): a cell without a valid SWAT HRU is excluded
                    // from the average. Do not substitute UZ recharge or another weight.
                    if (cell_debug.hruidx < 0) {
                        c_cell = initConc;
                        water_cell = 0.0;
                        cell_debug.regime = "missing_hru";
                    }
                    else {
                        // SWAT percolation is always the cell weight, including initial,
                        // historic, transition, and future loading regimes.
                        water_cell = swat.perc_mm(cell_debug.hruidx, iswat);
                        if (water_cell < 0.0) {
                            water_cell = 0.0;
                        }
                        cell_debug.perc_mm = swat.perc_mm(cell_debug.hruidx, iswat);

                        // OLD (pre-v0.1.16): historic and initial cells used UZ recharge.
                        // water_cell = UZ.getRch(IJ);

                    // Unsaturated-zone travel time shift (in years/index units).
                    // The loading that arrives now at the water table comes from an earlier time.
                    int shift = UZ.traveltime(IJ);
                    cell_debug.shift = shift;

                    // If travel time is longer than elapsed simulation time,
                    // then historic/future loading has not yet arrived; use initial concentration.
                    if (shift > iyr){
                        c_cell = initConc;
                        cell_debug.regime = "initial";
                    }
                    else{
                        // After unsaturated travel time has passed, we use historic loading from year zero
                        // Time index after removing unsaturated travel time.
                        int shifted_iyr = iyr - shift;
                        cell_debug.shifted_iyr = shifted_iyr;

                        // Period 1: historic loading only
                        if (shifted_iyr < UI.simOptions.nYears_historic){
                            cell_debug.regime = "historic";
                            // We are in the historic loading period
                            c_cell = histload.calculateConc(IJ, yyyy - shift);
                            cell_debug.c_hist = c_cell;

                            // Optional river-influence blending
                            if (ver == 1){
                                c_cell = (1 - strm.rivInfl) * c_cell + strm.rivInfl * UI.riverOptions.ConcValue;
                            }
                        }

                        // Period 2: transition/blending from historic to future loading
                        else if (shifted_iyr < UI.simOptions.nYears_blendEnd){
                            cell_debug.regime = "blend";
                            // We are in the blending period
                            // OLD (pre-v0.1.16): transition weighting blended UZ recharge
                            // and future percolation. The weight now remains SWAT perc_mm.
                            // const double historic_water_cell = water_cell;

                            // Historic value at shifted time
                            double c_cell_hist = histload.calculateConc(IJ, yyyy - shift);

                            if (ver == 1){
                                c_cell_hist = (1 - strm.rivInfl) * c_cell_hist + strm.rivInfl * UI.riverOptions.ConcValue;
                            }

                            // Future value based on SWAT + pumping-derived groundwater term
                            double c_cell_fut = 0.0;
                            double gw_cell_fut = 0.0;
                            double future_water_cell = 0.0;
                            double cfP = ConcFromPump[IJ];
                            double surf_perc = UZ.getSurfPerc(IJ);

                            const bool has_swat_water = calculateFutureLoading(
                                                   c_cell_fut, gw_cell_fut, future_water_cell,
                                                   IJ, iswat, ver, initConc, cfP, strm.rivInfl,
                                                   UI.riverOptions.ConcValue, UI.simOptions.SurfConcValue,
                                                   surf_perc, hru_raster, swat, UI.miscOptions.dp_mult,
                                                   debug_data);

                            // OLD draft of v0.1.12: switching immediately to the future water
                            // volume made the weighting discontinuous at the start of blending.
                            // if (has_swat_water){
                            //     water_cell = future_water_cell;
                            // }

                            // Linear blend factor from 0 (start of blend) to 1 (end of blend)
                            double u = (static_cast<double>(shifted_iyr - UI.simOptions.nYears_historic)) /
                                    (static_cast<double>(UI.simOptions.nYears_blendEnd - UI.simOptions.nYears_historic));

                            cell_debug.blend_u = u;
                            cell_debug.c_hist = c_cell_hist;
                            cell_debug.c_future = c_cell_fut;

                            c_cell = (1-u) * c_cell_hist + u * c_cell_fut;
                            gw_cell = gw_cell_fut * u;

                            // CHANGE (v0.1.16): never blend or replace the SWAT perc_mm
                            // weight. A failed HRU lookup excludes the cell.
                            if (!has_swat_water){
                                water_cell = 0.0;
                            }

                            // OLD (pre-v0.1.16):
                            // water_cell = (1-u) * historic_water_cell + u * future_water_cell;
                        }

                        // Period 3: future loading only
                        else{
                            cell_debug.regime = "future";
                            // We are in the future loading period
                            double cfP = ConcFromPump[IJ];
                            double surf_perc = UZ.getSurfPerc(IJ);

                            double future_water_cell = 0.0;
                            const bool has_swat_water = calculateFutureLoading(
                                                   c_cell, gw_cell, future_water_cell,
                                                   IJ, iswat, ver, initConc, cfP, strm.rivInfl,
                                                   UI.riverOptions.ConcValue, UI.simOptions.SurfConcValue,
                                                   surf_perc, hru_raster, swat, UI.miscOptions.dp_mult,
                                                   debug_data);

                            // CHANGE (v0.1.16): keep the SWAT perc_mm weight initialized
                            // above. If the HRU lookup fails, exclude the cell.
                            if (!has_swat_water){
                                water_cell = 0.0;
                            }

                            // OLD (pre-v0.1.16): water_cell = future_water_cell;
                        }
                    }
                    } // valid SWAT HRU
                }

                // CHANGE (v0.1.12): accumulate concentration times water volume. For future
                // loading without concentration-level mixing this reduces exactly to summing
                // m_total over cells and dividing by their total percolating-water volume.
                if (water_cell > 0.0){
                    weighted_c_sum = weighted_c_sum + c_cell * water_cell;
                    weighted_gw_sum = weighted_gw_sum + gw_cell * water_cell;
                    total_water = total_water + water_cell;
                }

                if (debug_file != nullptr) {
                    cell_debug.water_weight = water_cell;
                    cell_debug.c_cell = c_cell;
                    cell_debug.c_gw = gw_cell;
                    writeDebugCellRow(*debug_file, iyr, yyyy, iswat, debug_eid,
                                      strm.Sid, buffer_cell, ii, jj, IJ,
                                      strm.W, initConc,
                                      strm.rivInfl, cell_debug);
                }

                buffer_cell = buffer_cell + 1;

                // OLD (pre-v0.1.12): unweighted arithmetic accumulation.
                // c_strml = c_strml + c_cell;
                // gw_strml = gw_strml + gw_cell;
                // n_cells = n_cells + 1.0;
            }// loop jj
        }// loop ii

        if (total_water > 0.0){
            c_strml = weighted_c_sum / total_water;
            gw_strml = weighted_gw_sum / total_water;
        }
        else{
            // No cell in the neighborhood represents a positive water volume.
            c_strml = initConc;
            gw_strml = 0.0;
        }

        // OLD (pre-v0.1.12): unweighted arithmetic average.
        // c_strml = c_strml / n_cells;
        // gw_strml = gw_strml / n_cells;

        // CHANGE (v0.1.15): a negative MaxConc disables concentration capping.
        // OLD (pre-v0.1.15): if (c_strml > UI.simOptions.MaxConc){
        if (UI.simOptions.MaxConc >= 0.0 && c_strml > UI.simOptions.MaxConc){
            c_strml = UI.simOptions.MaxConc;
        }
    }

    void BlendLoading(double &c_final, double &u, HistoricOptions &h, double &c_hist, double &c_fut, int YYYY){
        if (YYYY <= h.BlendStart){
            c_final = c_hist;
            u = 0.0;
        }
        else if (YYYY >= h.BlendEnd){
            c_final = c_fut;
            u = 1.0;
        }
        else{
            u = (YYYY - h.BlendStart) / (h.BlendEnd - h.BlendStart);
            c_final = (1-u) * c_hist + u * c_fut;
        }
    }






}

#endif //MANTISSA_MS_LOADING_H
