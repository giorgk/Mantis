//
// Created by giorg on 5/30/2025.
//

#ifndef MANTISSA_MS_LOADING_H
#define MANTISSA_MS_LOADING_H

namespace MS{
    class HistoricLoading{
    public:
        HistoricLoading(){}
        int StartYear;
        int EndYear;
        int Interval;
        int YearA;
        int YearB;
        std::string Prefix;
        std::string Ext;
        std::vector<double> A;
        std::vector<double> B;
        bool Setup(int sy, int ey, int dy, std::string pref, std::string ext, boost::mpi::communicator &world);
        bool InitialLoading(boost::mpi::communicator &world);
        bool LoadNext(boost::mpi::communicator &world);
        double calculateConc(int IJ, int year,boost::mpi::communicator &world);

        bool readRasterLoad(std::string filename, std::vector<double> &V,boost::mpi::communicator &world);
    };

    bool HistoricLoading::InitialLoading(boost::mpi::communicator &world) {
        std::string fileA = Prefix + std::to_string(StartYear) + "." + Ext;
        std::string fileB = Prefix + std::to_string(StartYear) + "." + Ext;

        bool tf = readRasterLoad(fileA, A, world);
        if (!tf){
            return tf;
        }
        tf = readRasterLoad(fileB, B, world);
        return tf;
    }

    double HistoricLoading::calculateConc(int IJ, int year, boost::mpi::communicator &world) {
        double conc = 0;
        if (year >= YearA && year <= YearB){
            double u = (static_cast<double>(year) - static_cast<double>(YearA))/
                       (static_cast<double>(YearB) - static_cast<double>(YearA));
            conc = (1-u) * A[IJ] + u * B[IJ];
        }
        else if (year < YearA && YearA == StartYear){
            conc = A[IJ];
        }
        else if (year > YearB && YearB == EndYear){
            conc = B[IJ];
        }
        else{
            LoadNext(world);
            conc = calculateConc(IJ, year, world);
        }
        return conc;
    }

    bool HistoricLoading::readRasterLoad(std::string filename, std::vector<double> &V, boost::mpi::communicator &world) {
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
#if _USEHF>0
            const std::string NameSet("HRUs");
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

    bool HistoricLoading::LoadNext(boost::mpi::communicator &world) {
        A.clear();
        A = std::move(B);
        YearA = YearB;
        YearB = YearB + Interval;
        std::string fileB = Prefix + std::to_string(StartYear) + "." + Ext;
        bool tf = readRasterLoad(fileB, B, world);

        return tf;

    }

    bool HistoricLoading::Setup(int sy, int ey, int dy,
                                std::string pref, std::string ext,
                                boost::mpi::communicator &world) {
        StartYear = sy;
        EndYear = ey;
        Interval = dy;
        Prefix = pref;
        Ext = ext;
        YearA = StartYear;
        YearB = YearA + Interval;
        bool tf = InitialLoading(world);
        return tf;
    }

    void CreateLoadingForStreamline(double &c_strml, double &gw_strml, const int &iyr, const int &iswat,
                                      const STRML &strm, const UserInput &UI, const double &initConc,
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
            for (int jj = strm.urfI - UI.simOptions.nBuffer; jj <= strm.urfI + UI.simOptions.nBuffer; ++jj){
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
                        int hru = hru_raster.getHRU(IJ);
                        int hruidx = swat.hru_index(hru);
                        if (hruidx < 0){
                            c_cell = initConc;
                        }
                        else{
                            if (swat.perc_mm[iswat][hruidx] < 0.01){
                                c_cell = swat.Salt_perc_ppm[iswat][hruidx];
                                gw_cell = 0.0;
                            }
                            else{
                                double m_total = 0;
                                gw_cell = ConcFromPump[IJ];
                                double m_npsat = ConcFromPump[IJ] * swat.irrGW_mm[iswat][hruidx] / 100.0;
                                if (m_npsat < 0.01){
                                    gw_cell = 0.0;
                                    if (m_npsat < 0.0){
                                        m_npsat = 0;
                                    }
                                }
                                m_total = swat.irrsaltSW_kgha[iswat][hruidx] +
                                          m_npsat +
                                          swat.fertsalt_kgha[iswat][hruidx] +
                                          swat.dssl_kgha[iswat][hruidx] -
                                          swat.Qsalt_kgha[iswat][hruidx] -
                                          swat.uptk_kgha[iswat][hruidx];// -
                                          //swat.dSoilSalt_kgha[iswat][hruidx];
                                m_total = m_total * swat.pGW[iswat][hruidx];
                                if (m_total < 0){
                                    m_total = 0.0;
                                }
                                c_cell = m_total * 100 / swat.perc_mm[iswat][hruidx];
                            }

                            if (UI.npsatOptions.version == 1){
                                double u = calcRiverInfluence(strm.rivRist, UI.riverOptions);
                                c_cell = u * c_cell + (1-u) * UI.riverOptions.ConcValue;
                            }

                            if (UI.simOptions.SurfConcValue > 0) {
                                double surf_perc = UZ.getSurfPerc(IJ);
                                if (surf_perc > 0) {
                                    c_cell = c_cell * (1 - surf_perc) + surf_perc * UI.simOptions.SurfConcValue;
                                }
                            }
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
