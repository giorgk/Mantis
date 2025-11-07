//
// Created by giorg on 5/30/2025.
//

#ifndef MANTISSA_MS_LOADING_H
#define MANTISSA_MS_LOADING_H

namespace MS{
    class HistoricLoading{
    public:
        HistoricLoading(){}
        //int YearA;
        //int YearB;
        //int idx;
        std::string Prefix;
        std::string Ext;
        std::vector<std::vector<double>> HistLoad;
        std::vector<int> Years;
        bool Setup(std::string pref, std::string ext, boost::mpi::communicator &world);
        bool ReadLoading(boost::mpi::communicator &world);
        double calculateConc(int IJ, int year) const;

        bool readInfo(std::string filename, boost::mpi::communicator &world);
        bool readRasterLoad(std::string filename, std::vector<double> &V,boost::mpi::communicator &world);
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
        double conc = 0;
        if (year <= Years[0]){
            conc = HistLoad[0][IJ];
        }
        else if(year >= Years[Years.size()-1]){
            conc = HistLoad[Years.size()-1][IJ];
        }
        else{
            for (unsigned int i = 0; i < Years.size()-1; ++i ){
                if (year >= Years[i] && year <= Years[i+1]){
                    double u = (static_cast<double>(year) - static_cast<double>(Years[i]))/
                               (static_cast<double>(Years[i+1]) - static_cast<double>(Years[i]));
                    conc = (1-u) * HistLoad[i][IJ] + u * HistLoad[i+1][IJ];
                    break;
                }
            }
        }
        return conc;
    }

    bool HistoricLoading::readRasterLoad(std::string filename, std::vector<double> &V, boost::mpi::communicator &world) {
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
#if _USEHF>0
            const std::string NameSet("HIST");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
            dataset.read(V);
            return true;
#endif
        }
        else{
            std::vector<std::vector<double>> tmp;
            bool tf = RootReadsMatrixFileDistrib<double>(filename, tmp, 1, true, world, 5000000);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(tmp,world,0,1,345,353);}
            V = tmp[0];
            return true;
        }
    }

    bool HistoricLoading::readInfo(std::string filename, boost::mpi::communicator &world) {
        std::vector<std::vector<int>> tmp;
        bool tf = RootReadsMatrixFileDistrib<int>(filename, tmp, 1, true, world, 5000000);
        if (!tf){ return false;}
        Years = tmp[0];
        return true;
    }

    bool HistoricLoading::Setup(std::string pref, std::string ext,
                                boost::mpi::communicator &world) {
        Prefix = pref;
        Ext = ext;
        std::string finfo = Prefix + "info.dat";
        bool tf = readInfo(finfo,world);
        tf = ReadLoading(world);
        return tf;
    }

    void calculateFutureLoading(double &c_cell, double &gw_cell,
                                const int &IJ, const int &iswat, const int &ver,
                                const double &initConc, const double &concfromPump, const double &RivInfl,
                                const double &rivConcVal, const double &SurfConcValue, const double &surf_perc,
                                const MS::HRU_Raster &hru_raster,
                                const MS::SWAT_data &swat,
                                const double &dp_mult
                                //const MS::MiscOptions &miscopt
                                ){
        c_cell = 0.0;
        gw_cell = 0.0;
        int hru = hru_raster.getHRU(IJ);
        int hruidx = swat.hru_index(hru);
        if (hruidx < 0){
            c_cell = initConc;
        }
        else{
            if (swat.perc_mm[hruidx][iswat] < 0.01){
                c_cell = swat.Salt_perc_ppm[hruidx][iswat];
                gw_cell = 0.0;
                if (ver == 1){
                    c_cell = (1 - RivInfl) * c_cell + RivInfl * rivConcVal;
                }
            }
            else{
                double m_total = 0;
                gw_cell = concfromPump;
                double m_npsat = concfromPump * swat.irrGW_mm[hruidx][iswat] / 100.0;
                if (m_npsat < 0.01){
                    gw_cell = 0.0;
                    if (m_npsat < 0.0){
                        m_npsat = 0;
                    }
                }
                m_total = swat.irrsaltSW_kgha[hruidx][iswat] +
                          m_npsat +
                          swat.fertsalt_kgha[hruidx][iswat] +
                          swat.dssl_kgha[hruidx][iswat] -
                          swat.Qsalt_kgha[hruidx][iswat] -
                          swat.uptk_kgha[hruidx][iswat];// -
                //swat.dSoilSalt_kgha[hruidx][iswat];
                m_total = m_total * swat.pGW[hruidx][iswat];
                if (m_total < 0){
                    m_total = 0.0;
                }
                c_cell = m_total * 100 / (swat.perc_mm[hruidx][iswat] * dp_mult);
            }
            if (ver == 1){
                c_cell = (1 - RivInfl) * c_cell + RivInfl * rivConcVal;
            }

            if (SurfConcValue > 0 & surf_perc > 0) {
                c_cell = c_cell * (1 - surf_perc) + surf_perc * SurfConcValue;
            }
        }
    }

    void CreateLoadingForStreamline(double &c_strml, double &gw_strml, const int &iyr, const int &iswat, const int &yyyy, const int &ver,
                                      const STRML &strm, const UserInput &UI, const double &initConc,
                                      const HistoricLoading &histload,
                                      const MS::BackgroundRaster &backRaster,
                                      const MS::UNSAT &UZ, const MS::HRU_Raster &hru_raster,
                                      const MS::SWAT_data &swat,
                                      const std::vector<double> &ConcFromPump){
        c_strml = 0.0;
        gw_strml = 0.0;

        if (strm.inRiv){
            c_strml = UI.riverOptions.ConcValue;
            return;
        }
        if (strm.a > UI.simOptions.MaxAge){
            c_strml = initConc;
            return;
        }
        double n_cells = 0.0;
        for (int ii = strm.urfI - UI.simOptions.nBuffer; ii <= strm.urfI + UI.simOptions.nBuffer; ++ii){
            for (int jj = strm.urfJ - UI.simOptions.nBuffer; jj <= strm.urfJ + UI.simOptions.nBuffer; ++jj){
                double c_cell = 0.0;
                double gw_cell = 0.0;
                int IJ = backRaster.IJ(ii, jj);
                if (IJ < 0){
                    c_cell = initConc;
                }
                else{
                    int shift = UZ.traveltime(IJ);
                    if (shift > iyr){
                        c_cell = initConc;
                    }
                    else{
                        // After unsaturated travel time has passed, we use historic loading from year zero
                        int shifted_iyr = iyr - shift;
                        if (shifted_iyr < UI.simOptions.nYears_historic){
                            // We are in the historic loading period
                            c_cell = histload.calculateConc(IJ, yyyy - shift);
                            if (ver == 1){
                                c_cell = (1 - strm.rivInfl) * c_cell + strm.rivInfl * UI.riverOptions.ConcValue;
                            }
                        }
                        else if (shifted_iyr < UI.simOptions.nYears_blendEnd){
                            // We are in the blending period
                            double c_cell_hist = histload.calculateConc(IJ, yyyy - shift);
                            if (ver == 1){
                                c_cell_hist = (1 - strm.rivInfl) * c_cell_hist + strm.rivInfl * UI.riverOptions.ConcValue;
                            }

                            double c_cell_fut, gw_cell_fut;
                            double cfP = ConcFromPump[IJ];
                            double surf_perc = UZ.getSurfPerc(IJ);
                            calculateFutureLoading(c_cell_fut, gw_cell_fut, IJ, iswat, ver, initConc, cfP,
                                                   strm.rivInfl, UI.riverOptions.ConcValue, UI.simOptions.SurfConcValue,
                                                   surf_perc, hru_raster, swat, UI.miscOptions.dp_mult);
                            double u = (static_cast<double>(shifted_iyr - UI.simOptions.nYears_historic)) /
                                    (static_cast<double>(UI.simOptions.nYears_blendEnd - UI.simOptions.nYears_historic));
                            c_cell = (1-u) * c_cell_hist + u * c_cell_fut;
                            gw_cell = gw_cell_fut * u;
                        }
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
                c_strml = c_strml + c_cell;
                gw_strml = gw_strml + gw_cell;
                n_cells = n_cells + 1.0;
            }// loop jj
        }// loop ii
        c_strml = c_strml / n_cells;
        gw_strml = gw_strml / n_cells;

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
