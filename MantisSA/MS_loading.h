//
// Created by giorg on 5/30/2025.
//

#ifndef MANTISSA_MS_LOADING_H
#define MANTISSA_MS_LOADING_H

#include <random>

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
    //   c_cell  = final concentration assigned to this cell for the future-loading regime
    //   gw_cell = groundwater-irrigation concentration contribution diagnostic
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
    void calculateFutureLoading(double &c_cell, double &gw_cell, const int &IJ, const int &iswat,
                                const int &ver, const double &initConc, const double &concfromPump, const double &RivInfl,
                                const double &rivConcVal, const double &SurfConcValue, const double &surf_perc,
                                const MS::HRU_Raster &hru_raster,
                                const MS::SWAT_data &swat,
                                const double &dp_mult
                                //const MS::MiscOptions &miscopt
                                ){
        // Default outputs
        c_cell = 0.0;
        gw_cell = 0.0;

        // Map model cell to SWAT HRU
        const int hru = hru_raster.getHRU(IJ);
        const int hruidx = swat.hru_index(hru);

        // If no valid HRU exists, fall back to the initial concentration
        if (hruidx < 0){
            c_cell = initConc;
            return;
        }

        // Very small percolation:
        // avoid unstable division by a tiny number and use SWAT's percolate concentration directly.
        if (swat.perc_mm(hruidx, iswat) < 0.01){
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

            // Total salt mass available to deep percolation.
            // Positive sources:
            //   - irrigation surface-water salt
            //   - groundwater-irrigation salt (m_npsat)
            //   - fertilizer salt
            //   - dissolution / other source term
            // Negative sinks:
            //   - runoff/export
            //   - plant uptake
            double m_total = swat.irrsaltSW_kgha(hruidx, iswat) +
                      m_npsat +
                      swat.fertsalt_kgha(hruidx, iswat) +
                      swat.dssl_kgha(hruidx, iswat) -
                      swat.Qsalt_kgha(hruidx, iswat) -
                      swat.uptk_kgha(hruidx, iswat);// -
                    //swat.dSoilSalt_kgha[hruidx][iswat];

            // Apply the groundwater partition / weighting term used by the model.
            m_total = m_total * swat.pGW(hruidx, iswat);

            // Do not allow negative total mass
            if (m_total < 0){
                m_total = 0.0;
            }

            // Convert mass loading to concentration.
            // perc_mm is the deep percolation denominator,
            // dp_mult is an additional scaling / dilution factor.
            //double c_cell_tmp = m_total * 100 / (swat.perc_mm[hruidx][iswat]);
            c_cell = m_total * 100 / (swat.perc_mm(hruidx, iswat) * dp_mult);
        }

        // Optional river-influence blending.
        // ver == 1 activates mixing between computed cell concentration and river concentration.
        if (ver == 1){
            c_cell = (1 - RivInfl) * c_cell + RivInfl * rivConcVal;
        }

        // Optional surface-source blending.
        if (SurfConcValue > 0 && surf_perc > 0) {
            c_cell = c_cell * (1 - surf_perc) + surf_perc * SurfConcValue;
        }
    }

    void CreateLoadingForStreamline(double &c_strml, double &gw_strml, double &m_rmv_strml,
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
                                    const std::vector<double> &ConcFromPump){

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

        // We will average concentration over a square neighborhood centered at (urfI, urfJ).
        double n_cells = 0.0;

        for (int ii = strm.urfI - UI.simOptions.nBuffer; ii <= strm.urfI + UI.simOptions.nBuffer; ++ii){
            for (int jj = strm.urfJ - UI.simOptions.nBuffer; jj <= strm.urfJ + UI.simOptions.nBuffer; ++jj){
                double c_cell = 0.0; // concentration assigned to this neighborhood cell
                double gw_cell = 0.0; // groundwater contribution diagnostic for this cell
                double m_remove_cell = 0.0;

                // Convert raster coordinates to active-cell index.
                int IJ = backRaster.IJ(ii, jj);

                // If the raster cell is inactive / outside domain, use initial concentration.
                if (IJ < 0){
                    c_cell = initConc;
                }
                else{
                    // Unsaturated-zone travel time shift (in years/index units).
                    // The loading that arrives now at the water table comes from an earlier time.
                    int shift = UZ.traveltime(IJ);

                    // If travel time is longer than elapsed simulation time,
                    // then historic/future loading has not yet arrived; use initial concentration.
                    if (shift > iyr){
                        c_cell = initConc;
                    }
                    else{
                        // After unsaturated travel time has passed, we use historic loading from year zero
                        // Time index after removing unsaturated travel time.
                        int shifted_iyr = iyr - shift;

                        // Period 1: historic loading only
                        if (shifted_iyr < UI.simOptions.nYears_historic){
                            // We are in the historic loading period
                            c_cell = histload.calculateConc(IJ, yyyy - shift);

                            // Optional river-influence blending
                            if (ver == 1){
                                c_cell = (1 - strm.rivInfl) * c_cell + strm.rivInfl * UI.riverOptions.ConcValue;
                            }
                        }

                        // Period 2: transition/blending from historic to future loading
                        else if (shifted_iyr < UI.simOptions.nYears_blendEnd){
                            // We are in the blending period
                            // Historic value at shifted time
                            double c_cell_hist = histload.calculateConc(IJ, yyyy - shift);

                            if (ver == 1){
                                c_cell_hist = (1 - strm.rivInfl) * c_cell_hist + strm.rivInfl * UI.riverOptions.ConcValue;
                            }

                            // Future value based on SWAT + pumping-derived groundwater term
                            double c_cell_fut = 0.0;
                            double gw_cell_fut = 0.0;
                            double m_remove = 0.0;
                            double cfP = ConcFromPump[IJ];
                            double surf_perc = UZ.getSurfPerc(IJ);

                            calculateFutureLoading(c_cell_fut, gw_cell_fut, IJ, iswat, ver, initConc, cfP,
                                                   strm.rivInfl, UI.riverOptions.ConcValue, UI.simOptions.SurfConcValue,
                                                   surf_perc, hru_raster, swat, UI.miscOptions.dp_mult);

                            // Linear blend factor from 0 (start of blend) to 1 (end of blend)
                            double u = (static_cast<double>(shifted_iyr - UI.simOptions.nYears_historic)) /
                                    (static_cast<double>(UI.simOptions.nYears_blendEnd - UI.simOptions.nYears_historic));

                            c_cell = (1-u) * c_cell_hist + u * c_cell_fut;
                            gw_cell = gw_cell_fut * u;
                        }

                        // Period 3: future loading only
                        else{
                            // We are in the future loading period
                            double cfP = ConcFromPump[IJ];
                            double surf_perc = UZ.getSurfPerc(IJ);

                            calculateFutureLoading(c_cell, gw_cell, IJ, iswat, ver, initConc, cfP,
                                                   strm.rivInfl, UI.riverOptions.ConcValue, UI.simOptions.SurfConcValue,
                                                   surf_perc, hru_raster, swat, UI.miscOptions.dp_mult);
                        }
                    }
                }

                // Accumulate for neighborhood averaging
                c_strml = c_strml + c_cell;
                gw_strml = gw_strml + gw_cell;
                m_rmv_strml = m_rmv_strml + m_remove_cell;
                n_cells = n_cells + 1.0;
            }// loop jj
        }// loop ii
        c_strml = c_strml / n_cells;
        gw_strml = gw_strml / n_cells;
        m_rmv_strml = m_rmv_strml / n_cells;

        if (c_strml > UI.simOptions.MaxConc){
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
