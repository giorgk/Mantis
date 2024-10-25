//
// Created by giorg on 2/25/2023.
//

#ifndef MANTISSERVER_WELLS_H
#define MANTISSERVER_WELLS_H

#include <algorithm>

#include "MShelper.h"
#include "BMaps.h"
#include "BRaster.h"
#include "Rch.h"

namespace mantisServer{

    // These are tmp for deletion
    int n_noSource = 0;
    int n_total = 0;

    /**
	 * @brief Stores data for each streamline.
	 *
	 * Although this is a class, it is used more like struct container.
	 *
	 */
    class Streamline{
    public:
        Streamline(){};
        /*!
         *
         * @param row_ind URF ending point row index
         * @param col_ind URF ending point column index
         * @param riv_in True if source comes from river
         * @param npxl_in Number of source area pixels
         * @param w_in weight
         * @param len_in lenth of streamline
         * @param mean_in a list of fitted means
         * @param std_in a list of fitted stds
         * @param age_in a list of travel times
         */
        void setParameters(int row_ind, int col_ind, int riv_in, int npxl_in, double w_in,  double len_in,
                           std::vector<double> &mean_in, std::vector<double> &std_in, std::vector<double> &age_in);

        //void addSourceAreaCell(int lin_ind, int row, int col);
        void addSourceAreaCell(cell c);
        void clearSourceArea();
        void log();

        void print(int i);

        //! the row number of the pixel where this streamline starts from near the land surface.
        int row;

        //! the column number of the pixel where this streamline starts from near the land surface.
        int col;

        //! the mean value of the fitted unit response function.
        std::vector<double> mu;

        //! the standard deviation of the fitted unit response function.
        std::vector<double> std;

        std::vector<double> age;

        //!the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
        double w;

        //! The streamline length
        double len;

        //! the mean velocity along the streamline
        //std::vector<double> vel;

        //! This is true for the streamlines that start from rivers
        bool inRiver;

        int Npxl;


        bool bisZeroLoading = false;

        //std::map<int,cell> SourceArea;
        std::vector<cell> SourceArea;
    };

    void Streamline::clearSourceArea() {
        SourceArea.clear();
    }

    //void Streamline::addSourceAreaCell(int lin_ind, int row, int col) {
    //    addSourceAreaCell(lin_ind, cell(row,col));
    //}

    void Streamline::addSourceAreaCell(cell c) {
        for (unsigned int i = 0; i < SourceArea.size(); ++i){
            if (SourceArea[i].lin_ind == c.lin_ind)
                return;
        }
        SourceArea.push_back(c);
        //std::map<int,cell>::iterator it;
        //it = SourceArea.find(lin_ind);
        //if (it == SourceArea.end()){
        //    SourceArea.insert(std::pair<int, cell>(lin_ind, c));
        //}
    }

    void Streamline::print(int i) {
        std::cout << "r: " << row << ", c: " << col << ", m: " << mu[i] << ", s: " << std[i] << ", w: " << w << ", riv: " << inRiver << ", N: " << Npxl << std::endl;
    }

    void Streamline::setParameters(int row_ind, int col_ind, int riv_in, int npxl_in, double w_in,  double len_in,
                                   std::vector<double> &mean_in, std::vector<double> &std_in, std::vector<double> &age_in) {
        row = row_ind;
        col = col_ind;
        inRiver = riv_in == 1;
        Npxl = npxl_in;
        w = w_in;
        len = len_in;
        mu = mean_in;
        std = std_in;
        age = age_in;
    }

    void Streamline::log() {
        std::cout << row << ", " << col << ", "
        << inRiver << ", " << Npxl << ", " << w  << ", "
        << len << std::endl;
    }

    class Well{
    public:
        //Sid, row_ind, col_ind, riv, npxl, w, len,mean, std, age
        void addStreamline(int Sid, int row_ind, int col_ind, int riv, int npxl, double w, double len,
                           std::vector<double> &mean, std::vector<double> &std, std::vector<double> &age);
        void setAdditionalData(double x, double y, double u, double w, double s, double q, double r, double a);
        void calculateSourceArea(BackgroundRaster &braster, RechargeScenario &rch, int doCalc, std::ofstream &WSAstrm, bool debug = false, bool printWSA = false);
        void calculateSourceAreaType0(BackgroundRaster &braster, bool debug = false);
        void calculateSourceAreaType1(BackgroundRaster &braster, RechargeScenario &rch, bool debug = false);
        void calculateSourceAreaType2(BackgroundRaster &braster, RechargeScenario &rch, std::ofstream &WSAstrm, bool debug = false, bool printWSA = false);
        void calculateSourceAreaType3(BackgroundRaster &braster, RechargeScenario &rch, std::ofstream &WSAstrm, bool debug = false, bool printWSA = false);
        void calculateWeights();
        bool bsimulateThis(Scenario &scenario);
        void log();

        std::map<int, Streamline> streamlines;
        double xcoord;
        double ycoord;
        double unsat_depth;
        double wt2t;
        double screenLength;
        double pumpingRate;
        double ratio;
        double angle;
    };

    void Well::calculateWeights() {
        std::map<int, Streamline>::iterator it;
        double sumW = 0.0;
        for (it = streamlines.begin(); it != streamlines.end(); ++it){
            sumW = it->second.w + sumW;
        }
        for (it = streamlines.begin(); it != streamlines.end(); ++it){
            it->second.w = it->second.w/sumW;
        }
    }

    void Well::addStreamline(int Sid, int row_ind, int col_ind, int riv, int npxl, double w, double len,
                             std::vector<double> &mean, std::vector<double> &std, std::vector<double> &age) {
        Streamline s;
        s.setParameters(row_ind, col_ind, riv,npxl, w, len, mean, std, age);
        streamlines.insert(std::pair<int, Streamline>(Sid, s));
    }

    void Well::setAdditionalData(double x, double y, double u, double w, double s, double q, double r, double a){
        xcoord = x;
        ycoord = y;
        unsat_depth = u;
        wt2t = w;
        screenLength = s;
        pumpingRate = q;
        ratio = r;
        angle = a;
    }

    void Well::log() {
        std::cout << xcoord << ", " << ycoord << ", " << wt2t << ", "
        << screenLength << ", " << pumpingRate << ", " << ratio << ", " << angle <<std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::map<int, Streamline>::iterator it;
        for (it = streamlines.begin(); it != streamlines.end();++it){
            it->second.log();
        }
    }

    void Well::calculateSourceArea(BackgroundRaster &braster, RechargeScenario &rch, int doCalc, std::ofstream &WSAstrm, bool debug, bool printWSA){
        if (doCalc == 0){
            calculateSourceAreaType0(braster, debug);
        }
        else if (doCalc == 1){
            calculateSourceAreaType2(braster, rch, WSAstrm, debug, printWSA);
        }
        //else if (doCalc == 2){
        //    calculateSourceAreaType2(braster, rch, WSAstrm, debug, printWSA);
        //}
        //else if (doCalc == 3){
        //    calculateSourceAreaType3(braster, rch, WSAstrm, debug, printWSA);
        //}
    }

    void Well::calculateSourceAreaType0(BackgroundRaster &braster, bool debug){

        std::map<int, Streamline>::iterator itstrml;
        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            int rasterValue = braster.IJ(itstrml->second.row, itstrml->second.col);
            if (rasterValue == -1){
                // The source area is outside the active area therefore the loading will be zero
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                continue;
            }
            if (itstrml->second.inRiver){
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                continue;
            }
            cell c;
            c.row = itstrml->second.row;
            c.col = itstrml->second.col;
            c.lin_ind = rasterValue;
            c.dist = 0.0;
            itstrml->second.addSourceAreaCell(c);
            if (debug){
                std::cout << itstrml->first << std::endl;
            }

            // This is used only for debugging
            for (int i = -1; i < 2; ++i){
                for (int j = -1; j < 2; ++j){
                    cell tmp;
                    tmp.row = itstrml->second.row + i;
                    tmp.col = itstrml->second.col + j;
                    int rval = braster.IJ(tmp.row, tmp.col);
                    if (rval == -1){
                        continue;
                    }
                    tmp.lin_ind = rval;
                    tmp.dist = 0.0;
                    itstrml->second.addSourceAreaCell(tmp);
                }
            }
        }

    }

    /*void Well::calculateSourceAreaType1(BackgroundRaster &braster, RechargeScenario &rch, bool debug) {
        bool tf;
        int rs, cs, rn, cn, rasterValue;
        double Xorig, Yorig, cellSize, xs, ys, dNr, side1, side2;
        //double SAxmin, SAxmax, SAymin, SAymax;
        double x1, x2, x3, x4, y1, y2, y3, y4;
        double Qtmp, Qtarget, xtmp, ytmp, rch_v, dummy;
        braster.getGridLocation(Xorig, Yorig, cellSize);
        double cellArea = cellSize*cellSize;
        int Nrow = braster.Nr();
        dNr = static_cast<double>(Nrow);
        double radAngle = angle*pi/180.0;
        double cosd = std::cos(radAngle);
        double sind = std::sin(radAngle);
        std::map<int, Streamline>::iterator itstrml;
        std::map< int, cell>::iterator itcell1, itcell2;
        std::vector<cell> sp = SearchPattern();

        std::vector<boost_point> srcPnts;
        boost_poly srcPoly;


        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            if (debug){
                std::cout << itstrml->first << std::endl;
            }
            // The value of the raster points to the right row on the data that have linear form
            rasterValue = braster.IJ(itstrml->second.row, itstrml->second.col);
            if (rasterValue == -1){
                // The source area is outside the active area therefore the loading will be zero
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                continue;
            }

            if (itstrml->second.inRiver){
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                continue;
            }



            int niter = 0;
            bool sourceFound = false;
            while (niter < 20){
                rs = itstrml->second.row;
                cs = itstrml->second.col;
                double drs = static_cast<double>(rs);
                double dcs = static_cast<double>(cs);

                braster.cellCoords(rs, cs, xs, ys);
                //xs = Xorig + cellSize/2 + cellSize*(cs);
                // For the Y the row numbers start from the top
                //ys =  Yorig + cellSize*dNr - cellSize/2 - cellSize*(rs);
                side1 = cellSize + cellSize * (itstrml->second.Npxl+niter);
                side2 = std::max(2.0*cellSize, side1/ratio/2.0);
                double SAxmin = -side2 - 25;
                double SAxmax =  side2 + 25;
                double SAymin = -side1/2 - 25;
                double SAymax =  side1/2 + 25;

                x1 = (cosd * SAxmin + sind * SAymin) + xs;
                y1 = (-sind * SAxmin + cosd * SAymin) + ys;
                srcPnts.push_back(boost_point(x1, y1));

                x2 = (cosd * SAxmin + sind * SAymax) + xs;
                y2 = (-sind * SAxmin + cosd * SAymax) + ys;
                srcPnts.push_back(boost_point(x2, y2));

                x3 = (cosd * SAxmax + sind * SAymax) + xs;
                y3 = (-sind * SAxmax + cosd * SAymax) + ys;
                srcPnts.push_back(boost_point(x3, y3));

                x4 = (cosd * SAxmax + sind * SAymin) + xs;
                y4 = (-sind * SAxmax + cosd * SAymin) + ys;
                srcPnts.push_back(boost_point(x4, y4));
                boost::geometry::assign_points(srcPoly, srcPnts);
                boost::geometry::correct(srcPoly);

                if (debug){
                    std::cout << "pp=[" << x1 << " " << y1 << "; " << x2 << " " << y2 << "; "
                              << x3 << " " << y3 << "; " << x4 << " " << y4 << "]; " << std::endl;
                }

                std::map< int, cell> for_test;
                std::map< int, cell> tested;
                std::map< int, cell> next_round;

                //lin_ind = braster.linear_index(itstrml->second.row, itstrml->second.col);
                for_test.insert(std::pair<int, cell>(rasterValue, cell(rs,cs)));

                Qtmp = 0.0;
                Qtarget = std::abs(pumpingRate * itstrml->second.w);

                while (true){
                    for (itcell1 = for_test.begin(); itcell1 != for_test.end(); ++itcell1){
                        tested.insert(std::pair<int,cell>(itcell1->first, itcell1->second));
                        braster.cellCoords(itcell1->second.row, itcell1->second.col, xtmp, ytmp);
                        if (debug){
                            std::cout << "plot(" << xtmp << "," << ytmp << ",'.k');" << std::endl;
                        }

                        // Check if the point is in the source area by testing the barycentric coordinates
                        tf = boost::geometry::within(boost_point(xtmp, ytmp), srcPoly);
                        if (!tf){
                            continue;
                        }
                        //tf = isInTriangle(x1, y1, x2, y2, x3, y3, xtmp, ytmp);
                        //if (!tf){
                        //    tf = isInTriangle(x1, y1, x3, y3, x4, y4, xtmp, ytmp);
                        //    if (!tf){
                        //        continue;
                        //    }
                        //}
                        rasterValue = braster.IJ(itcell1->second.row, itcell1->second.col);
                        if (rasterValue != -1){
                            tf = rch.getValues(rasterValue, rch_v, dummy);
                            if (tf){
                                if (rch_v > 10.0){
                                    itcell1->second.lin_ind = rasterValue;
                                    itcell1->second.dist = (drs - static_cast<double>(itcell1->second.row))*(drs - static_cast<double>(itcell1->second.row)) +
                                                           (dcs - static_cast<double>(itcell1->second.col))*(dcs - static_cast<double>(itcell1->second.col));
                                    itstrml->second.addSourceAreaCell(itcell1->second);
                                    Qtmp += (rch_v/365/1000)*cellArea;
                                }
                            }
                            if (Qtmp >= Qtarget){
                                sourceFound = true;
                                break;
                            }
                            next_round.insert(std::pair<int,cell>(itcell1->first, itcell1->second));
                        }
                    }
                    if (sourceFound){
                        break;
                    }
                    else{
                        if (next_round.empty()){
                            break;
                        }
                        else{
                            for_test.clear();
                            for (itcell1 = next_round.begin(); itcell1 != next_round.end(); ++itcell1){
                                for (unsigned int i = 0; i < sp.size(); ++i){
                                    rn = itcell1->second.row + sp[i].row;
                                    cn = itcell1->second.col + sp[i].col;
                                    rasterValue = braster.IJ(rn,cn);
                                    if (rasterValue == -1)
                                        continue;
                                    itcell2 = tested.find(rasterValue);
                                    if (itcell2 == tested.end()){
                                        for_test.insert(std::pair<int,cell>(rasterValue, cell(rn,cn)));
                                    }
                                }
                            }
                            next_round.clear();
                        }
                    }
                }
                if (sourceFound){
                    break;
                }
                if (Qtmp < 0.0001){
                    // If the source area of this iteration comes from a
                    // zero recharge area we assume the water is clean
                    // if the source area is sufficiently large
                    if (itstrml->second.Npxl + niter > 21){
                        break;
                    }
                }
                niter = niter + 1;
                itstrml->second.clearSourceArea();
            }

        }
    }*/

    void Well::calculateSourceAreaType2(BackgroundRaster &braster, RechargeScenario &rch, std::ofstream &WSAstrm, bool debug, bool printWSA){
        std::map<int, Streamline>::iterator itstrml;
        double Xcntr, Ycntr, xcell, ycell, rch_v, dummy;
        //double SAxmin, SAxmax, SAymin, SAymax, minSrcXrt, minSrcYrt;
        double Xorig, Yorig, cellSize;
        braster.getGridLocation(Xorig, Yorig, cellSize);
        double cellArea = cellSize*cellSize;
        int Ncol = braster.Nc();
        int Nrow = braster.Nr();

        double M11, M12, M21, M22;//, Minv11, Minv12, Minv21, Minv22, Dinv;
        M11 = cosd(angle); M12 = sind(angle);
        M21 = -M12; M22 = M11;
        //Dinv = 1/(M11*M22 - M12*M21);
        //Minv11 = Dinv*M22; Minv12 = -Dinv*M12;
        //Minv21 = -Dinv*M21; Minv22 = Dinv*M11;
        std::vector<double> srcX, srcY, srcXrt, srcYrt;
        //std::vector<boost_point> srcPnts;
        //boost_poly srcPoly;

        //int maxcells = 0;

        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            //if (itstrml->second.inRiver || std::abs(itstrml->second.mu) < 0.0001){
            if (itstrml->second.inRiver){
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                //int rasterValue = braster.IJ(itstrml->second.row, itstrml->second.col);
                //if (rasterValue >= 0){
                //    cell c;
                //    c.row = itstrml->second.row;
                //    c.col = itstrml->second.col;
                //    c.lin_ind = rasterValue;
                //    c.dist = 0.0;
                //    itstrml->second.addSourceAreaCell(c);
                //    if (itstrml->second.SourceArea.size() > maxcells){
                //        maxcells = itstrml->second.SourceArea.size();
                //    }
                //}
                continue;
            }

            double Qtarget = std::abs(pumpingRate * itstrml->second.w);
            double Qtmp = 0.0;

            double minSrcXrt = 999999999999.9; double minSrcYrt = 999999999999.9;
            int row = itstrml->second.row;
            int col = itstrml->second.col;
            int count_iter = 0;
            braster.cellCoords(row, col, Xcntr, Ycntr);
            double Npxl = static_cast<double>(itstrml->second.Npxl);
            if (Npxl == 0){
                Npxl = 1;
            }
            while (count_iter < 1){
                itstrml->second.SourceArea.clear();
                Qtmp = 0.0;
                double lngth = cellSize + cellSize * Npxl;
                double wdth = std::max(lngth/ratio/1.95, cellSize/1.95);
                double SAxmin = -wdth;
                double SAxmax =  wdth;
                double SAymin = -lngth/1.95 - cellSize/1.95;
                double SAymax =  lngth/1.95 + cellSize/1.95;
                srcX.clear(); srcY.clear(); srcXrt.clear(), srcYrt.clear(); //srcPnts.clear();
                srcX.push_back(SAxmin); srcY.push_back(SAymin);
                srcX.push_back(SAxmin); srcY.push_back(SAymax);
                srcX.push_back(SAxmax); srcY.push_back(SAymax);
                srcX.push_back(SAxmax); srcY.push_back(SAymin);

                for (unsigned int i = 0; i < srcX.size(); ++i){
                    srcXrt.push_back(M11 * srcX[i] + M12 * srcY[i] + Xcntr);
                    srcYrt.push_back(M21 * srcX[i] + M22 * srcY[i] + Ycntr);
                    if (srcXrt[i] < minSrcXrt){
                        minSrcXrt = srcXrt[i];
                    }
                    if (srcYrt[i] < minSrcYrt){
                        minSrcYrt = srcYrt[i];
                    }
                    //srcPnts.push_back(boost_point(srcXrt[i], srcYrt[i]));

                }
                //boost::geometry::assign_points(srcPoly, srcPnts);
                //boost::geometry::correct(srcPoly);

                //std::cout << "srcX = [" << srcXrt[0] << "," << srcXrt[1] << "," << srcXrt[2] << "," << srcXrt[3] << "]" << std::endl;
                //std::cout << "scrY = [" << srcYrt[0] << "," << srcYrt[1] << "," << srcYrt[2] << "," << srcYrt[3] << "]" << std::endl;

                int dJ = static_cast<int>(std::ceil((Xcntr - minSrcXrt)/cellSize));
                int dI = static_cast<int>(std::ceil((Ycntr - minSrcYrt)/cellSize));
                int JgridStart = std::max(1, col - dJ);
                int JgridEnd = std::min(Ncol, col + dJ);
                int IgridStart = std::max(1, row - dI);
                int IgridEnd = std::min(Nrow, row + dI);

                //Debug only
                int cnt_cells = 0;
                bool done_that = false;

                for (int i = IgridStart; i <= IgridEnd; ++i){
                    for (int j = JgridStart; j <= JgridEnd; ++j){
                        braster.cellCoords(i, j, xcell, ycell);
                        //std::cout << xcell << "," << ycell << std::endl;
                        //bool tf = boost::geometry::within(boost_point(xcell, ycell), srcPoly);
                        bool tf = isInTriangle(srcXrt[0], srcYrt[0], srcXrt[1], srcYrt[1], srcXrt[2], srcYrt[2], xcell, ycell);
                        if (!tf){
                            tf = isInTriangle(srcXrt[0], srcYrt[0], srcXrt[2], srcYrt[2], srcXrt[3], srcYrt[3], xcell, ycell);
                            if (!tf){
                                continue;
                            }
                        }

                        if (tf){
                            int rasterValue = braster.IJ(i, j);
                            if (rasterValue != -1){
                                tf = rch.getValues(rasterValue, rch_v, dummy);
                                if (rch_v > 10.0){
                                    cell c;
                                    c.row = i;
                                    c.col = j;
                                    c.lin_ind = rasterValue;
                                    c.dist = static_cast<double>((i-row)*(i-row) + (j-col)*(j-col));
                                    itstrml->second.addSourceAreaCell(c);
                                    Qtmp += (rch_v/365/1000)*cellArea;

                                }
                            }
                        }
                    }
                }
                if (Qtmp >= Qtarget){
                    break;
                }
                // Double the source area and try again
                Npxl = Npxl * 2;
                count_iter = count_iter + 1;
            }
            if (itstrml->second.SourceArea.size() == 0){
                n_noSource++;
            }
            n_total++;
            std::sort(itstrml->second.SourceArea.begin(), itstrml->second.SourceArea.end(), compareCellByDistance);
            if (printWSA){
                WSAstrm << itstrml->first << " " << std::setprecision(2) << std::fixed
                        << srcXrt[0] << " " << srcYrt[0] << " "
                        << srcXrt[1] << " " << srcYrt[1] << " "
                        << srcXrt[2] << " " << srcYrt[2] << " "
                        << srcXrt[3] << " " << srcYrt[3] << std::endl;
            }
        }
    }

    void Well::calculateSourceAreaType3(BackgroundRaster &braster, RechargeScenario &rch, std::ofstream &WSAstrm, bool debug, bool printWSA){
        std::map<int, Streamline>::iterator itstrml;
        double Xcntr, Ycntr, xcell, ycell, rch_v, dummy;
        double Xorig, Yorig, cellSize;
        braster.getGridLocation(Xorig, Yorig, cellSize);
        int Ncol = braster.Nc();
        int Nrow = braster.Nr();

        double M11, M12, M21, M22;
        M11 = cosd(angle); M12 = sind(angle);
        M21 = -M12; M22 = M11;

        std::vector<double> srcX, srcY, srcXrt, srcYrt;

        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            //std::cout << itstrml->first << std::endl;
            //if (itstrml->second.inRiver || std::abs(itstrml->second.mu) < 0.0001){
            if (itstrml->second.inRiver){
                itstrml->second.bisZeroLoading = true;
                //itstrml->second.mu = 0.0;
                //itstrml->second.std = 0.0;
                continue;
            }
            double minSrcXrt = 999999999999.9; double minSrcYrt = 999999999999.9;

            int row = itstrml->second.row;
            int col = itstrml->second.col;
            braster.cellCoords(row, col, Xcntr, Ycntr);
            double Npxl = static_cast<double>(itstrml->second.Npxl);
            itstrml->second.SourceArea.clear();
            double lngth = cellSize + cellSize * Npxl;
            double wdth = std::max(lngth/ratio, cellSize);
            double SAxmin = -wdth;
            double SAxmax =  wdth;
            double SAymin = -lngth - cellSize;
            double SAymax =  lngth + cellSize;
            srcX.clear(); srcY.clear(); srcXrt.clear(), srcYrt.clear(); //srcPnts.clear();
            srcX.push_back(SAxmin); srcY.push_back(SAymin);
            srcX.push_back(SAxmin); srcY.push_back(SAymax);
            srcX.push_back(SAxmax); srcY.push_back(SAymax);
            srcX.push_back(SAxmax); srcY.push_back(SAymin);

            for (unsigned int i = 0; i < srcX.size(); ++i){
                srcXrt.push_back(M11 * srcX[i] + M12 * srcY[i] + Xcntr);
                srcYrt.push_back(M21 * srcX[i] + M22 * srcY[i] + Ycntr);
                if (srcXrt[i] < minSrcXrt){
                    minSrcXrt = srcXrt[i];
                }
                if (srcYrt[i] < minSrcYrt){
                    minSrcYrt = srcYrt[i];
                }
            }
            int dJ = static_cast<int>(std::ceil((Xcntr - minSrcXrt)/cellSize));
            int dI = static_cast<int>(std::ceil((Ycntr - minSrcYrt)/cellSize));
            int JgridStart = std::max(1, col - dJ);
            int JgridEnd = std::min(Ncol, col + dJ);
            int IgridStart = std::max(1, row - dI);
            int IgridEnd = std::min(Nrow, row + dI);
            for (int i = IgridStart; i <= IgridEnd; ++i){
                for (int j = JgridStart; j <= JgridEnd; ++j){
                    braster.cellCoords(i, j, xcell, ycell);
                    bool tf = isInTriangle(srcXrt[0], srcYrt[0], srcXrt[1], srcYrt[1], srcXrt[2], srcYrt[2], xcell, ycell);
                    if (!tf){
                        tf = isInTriangle(srcXrt[0], srcYrt[0], srcXrt[2], srcYrt[2], srcXrt[3], srcYrt[3], xcell, ycell);
                        if (!tf){
                            continue;
                        }
                    }

                    if (tf){
                        int rasterValue = braster.IJ(i, j);
                        if (rasterValue != -1){
                            tf = rch.getValues(rasterValue, rch_v, dummy);
                            if (rch_v > 0.0){
                                cell c;
                                c.row = i;
                                c.col = j;
                                c.lin_ind = rasterValue;
                                c.dist = static_cast<double>((i-row)*(i-row) + (j-col)*(j-col));
                                itstrml->second.addSourceAreaCell(c);
                            }
                        }
                    }
                }
            }
            //if (itstrml->second.SourceArea.size() == 0){
            //    std::cout << "Streamline Sid: " << itstrml->first << "has no source area" << std::endl;
            //}
        }
    }

    bool Well::bsimulateThis(Scenario &scenario) {
        if (scenario.bNarrowSelection){
            if (scenario.useRadSelect) {
                if (!scenario.RadSelect.isPointIn(xcoord, ycoord)) {
                    return false;
                }
            }
            if (scenario.useRectSelect){
                if (!scenario.RectSelect.isPointIn(xcoord, ycoord)){
                    return false;
                }
            }
            if (scenario.DepthRange.canUse()){
                if (!scenario.DepthRange.isInRange(unsat_depth + wt2t + screenLength)){
                    return false;
                }
            }
            if (scenario.unsatRange.canUse()){
                if (!scenario.unsatRange.isInRange(unsat_depth)){
                    return false;
                }
            }
            if (scenario.wt2tRange.canUse()){
                if (!scenario.wt2tRange.isInRange(wt2t)){
                    return false;
                }
            }
            if (scenario.ScreenLengthRange.canUse()){
                if (!scenario.ScreenLengthRange.isInRange(screenLength)){
                    return false;
                }
            }
        }
        return true;
    }

    class WellList{
    public:
        WellList(){}
        void addWell(int Eid, Well w);
        bool addstreamline(int Eid,int Sid, int row_ind, int col_ind, int riv, int npxl, double w, double len,
                           std::vector<double> mean, std::vector<double> std, std::vector<double> age);
        void setPorosityIndices(std::vector<int> &por);
        int calcSourceArea = false;
        int printSourceArea = false;
        std::string rch_map;
        std::map<int, Well> Wells;
        std::map<int,int> PorosityIndex;
    };
    //wellid, Sid ,row_ind, col_ind, riv, npxl, w, len, mean, std, age
    bool WellList::addstreamline(int Eid,int Sid, int row_ind, int col_ind, int riv, int npxl, double w, double len,
                                 std::vector<double> mean, std::vector<double> std, std::vector<double> age){

        std::map<int, Well>::iterator eidit;
        eidit = Wells.find(Eid);
        if (eidit == Wells.end()){
            std::cout << "I can't find a well with id [ " << Eid << " ]" << "in the list" << std::endl;
            return false;
        }
        eidit->second.addStreamline(Sid, row_ind, col_ind, riv, npxl, w, len,mean, std, age);
        return true;
    }
    void WellList::setPorosityIndices(std::vector<int> &por) {
        for (unsigned int i = 0; i < por.size(); ++i){
            PorosityIndex.insert(std::pair<int, int>(por[i],i));
        }
    }

    void WellList::addWell(int Eid, Well w) {
        Wells.insert(std::pair<int, Well>(Eid, w));
    }

    class FlowWellCollection{
    public:
        FlowWellCollection(){}

        bool readMainfile(std::string path, std::string filename,
                          bool bReadWells, BMapCollection &Bmaps);
        void calcWellWeights();
        void calcWellSourceArea(BackgroundRaster &braster, RechargeScenarioList &rchList);
        std::map<std::string ,WellList> FlowScenarios;
        bool hasFlowScenario(std::string &flowScen, int por, int &porind);
    private:
        bool readWells(std::string filename, BMapCollection &Bmaps);
        bool readURFs(std::string filename);
        //addStreamline(name, eid, sid, r, c, riv, npxl, w, len, M, S, A)
        bool addStreamline(std::string flowScenName, int wellid,
                           int Sid, int row_ind, int col_ind, int riv, int npxl,double w, double len,
                            std::vector<double> &mean, std::vector<double> &std, std::vector<double> &age);
        bool addPorosityScenarios(std::string flowScenName, std::vector<int> &por);
    };

    bool FlowWellCollection::hasFlowScenario(std::string &flowScen, int por, int &porind) {
        std::map<std::string ,WellList>::iterator flowit;
        flowit = FlowScenarios.find(flowScen);
        if (flowit != FlowScenarios.end()){
            std::map<int,int>::iterator porit;
            porit = flowit->second.PorosityIndex.find(por);
            if (porit != flowit->second.PorosityIndex.end()){
                porind = porit->second;
                return true;
            }
            else{
                return false;
            }
        }
        else{
            return false;
        }
    }

    void FlowWellCollection::calcWellWeights() {
        std::map<std::string ,WellList>::iterator flowit;
        std::map<int, Well>::iterator wellit;
        // First Calculate the streamline weights to sum to 1
        for (flowit = FlowScenarios.begin(); flowit != FlowScenarios.end(); ++flowit){
            for (wellit = flowit->second.Wells.begin(); wellit != flowit->second.Wells.end(); ++ wellit){
                wellit->second.calculateWeights();
            }
        }
    }

    void FlowWellCollection::calcWellSourceArea(BackgroundRaster &braster, RechargeScenarioList &rchList) {



        std::map<std::string ,WellList>::iterator flowit;
        std::map<int, Well>::iterator wellit;
        std::map<std::string, RechargeScenario>::iterator rchit;
        for (flowit = FlowScenarios.begin(); flowit != FlowScenarios.end(); ++flowit){
            rchit = rchList.RechargeList.find(flowit->second.rch_map);
            if (rchit == rchList.RechargeList.end()){
                std::cout << "I can't find a recharge scenario with name " << flowit->second.rch_map << std::endl;
                continue;
            }
            int count = 0;
            bool dbg = false;
            std::string WSAfn = flowit->first + "WSA.dat";
            std::ofstream WSAstrm;
            if (flowit->second.calcSourceArea == 2 && flowit->second.printSourceArea == 1){
                WSAstrm.open(WSAfn.c_str());
            }
            n_noSource = 0;
            n_total = 0;

            for (wellit = flowit->second.Wells.begin(); wellit != flowit->second.Wells.end(); ++ wellit){
                //std::cout << wellit->first << std::endl;
                //if (wellit->first == 9341){
                //    dbg = true;
                //    std::cout << "Stop here" << std::endl;
                //}
                //if (count == 5851){
                //    std::cout << "Stop here" << std::endl;
                //}
                if (flowit->second.calcSourceArea == 2 && flowit->second.printSourceArea == 1){
                    WSAstrm << -9 << " " << wellit->first <<  " 0 0 0 0 0 0 0" << std::endl;
                }
                wellit->second.calculateSourceArea(braster, rchit->second, flowit->second.calcSourceArea, WSAstrm, dbg, flowit->second.printSourceArea);
                //dbg = false;
                count = count + 1;
                if (count % 1000 == 0){
                    std::cout << "----" << count << "----" << std::endl;
                }
            }
            if (flowit->second.calcSourceArea == 2 && flowit->second.printSourceArea == 1){
                WSAstrm.close();
            }
            std::cout << n_noSource << " | " << n_total << std::endl;
        }


        //flowit->second.calcSourceArea
    }

    bool FlowWellCollection::addStreamline(std::string flowScenName, int wellid,
                                           int Sid, int row_ind, int col_ind, int riv, int npxl,double w, double len,
                                           std::vector<double> &mean, std::vector<double> &std, std::vector<double> &age){
        std::map<std::string ,WellList>::iterator flowit;
        flowit = FlowScenarios.find(flowScenName);
        if (flowit == FlowScenarios.end()){
            std::cout << "I can't find wells under the [ " << flowScenName << " ] flow scenario" << std::endl;
            return false;
        }
        //addStreamline(eid, sid, r, c, riv, npxl, w, len, M, S, A)
        flowit->second.addstreamline(wellid, Sid ,row_ind, col_ind, riv, npxl, w, len, mean, std, age);
        return true;
    }

    bool FlowWellCollection::addPorosityScenarios(std::string flowScenName, std::vector<int> &por) {
        std::map<std::string ,WellList>::iterator flowit;
        flowit = FlowScenarios.find(flowScenName);
        if (flowit == FlowScenarios.end()){
            std::cout << "I can't find wells under the [ " << flowScenName << " ] flow scenario" << std::endl;
            return false;
        }
        flowit->second.setPorosityIndices(por);
        return true;
    }

    bool FlowWellCollection::readMainfile(std::string path, std::string filename,
                                          bool bReadWells, BMapCollection &Bmaps) {
        std::ifstream mainfile;
        std::string fullfilename = path + filename;
        mainfile.open(fullfilename);
        if (!mainfile.is_open()) {
            std::cout << "Cant open file: " << fullfilename << std::endl;
            return false;
        }
        else{
            std::string line, filename1;
            while (getline(mainfile, line)){
                filename1.clear();
                if (line.empty())
                    continue;
                if (line.front() == '#')
                    continue;

                std::istringstream inp(line.c_str());
                inp >> filename1;

                if (filename1.empty())
                    continue;
                if (filename1.front() == '#')
                    continue;
                filename1 = path + filename1;
                bool tf;
                if (bReadWells){
                    tf = readWells(filename1,Bmaps);
                }
                else{
                    tf = readURFs(filename1);
                }

                if (!tf) {
                    std::cout << "An error occurred while reading " << filename1 << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

    bool FlowWellCollection::readWells(std::string filename, BMapCollection &Bmaps) {
        auto start = std::chrono::high_resolution_clock::now();
        std::ifstream Welldatafile;
        Welldatafile.open(filename);
        if (!Welldatafile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            int Nwells, Eid, calcSource, printSource;
            std::string setName, line, rch_scen;
            WellList wList;
            {
                getline(Welldatafile, line);
                std::istringstream inp(line.c_str());
                inp >> Nwells;
                inp >> setName;
                inp >> rch_scen;
                inp >> calcSource;
                inp >> printSource;
                wList.calcSourceArea = calcSource;
                wList.printSourceArea = printSource;
                wList.rch_map = rch_scen;
            }

            double xw, yw, unsat, wt2t, SL, Q, ratio, angle;
            std::string regionCode;
            bool tf;
            for (int i = 0; i < Nwells; ++i){
                getline(Welldatafile, line);
                std::istringstream inp(line.c_str());
                inp >> Eid;
                //std::cout << Eid << std::endl;
                inp >> xw;
                inp >> yw;
                inp >> unsat;
                inp >> wt2t;
                inp >> SL;
                inp >> Q;
                inp >> ratio;
                inp >> angle;
                Well w;
                w.setAdditionalData(xw,yw,unsat, wt2t,SL,Q,ratio,angle);
                wList.addWell(Eid, w);
                for (int j = 0; j < Bmaps.Nbmaps(); ++j){
                    inp >> regionCode;
                    tf = Bmaps.addwell(regionCode,setName,Eid);
                    if (!tf){
                        std::cout << "I cannot find a modelArea with name [ " << regionCode << " ]" << std::endl;
                        return false;
                    }
                }
            }
            FlowScenarios.insert(std::pair<std::string, WellList>(setName, wList));
        }


        Welldatafile.close();
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Wells in " << elapsed.count() << std::endl;
        return true;
    }

    bool FlowWellCollection::readURFs(std::string filename) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Reading " << filename << std::endl;
#if _USEHF>0
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
            const std::string NamesNameSet("Names");
            const std::string IntsNameSet("ESIJRN");
            const std::string FloatNameSet("WLMSA");
            const std::string PORNameSet("POR");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
            HighFive::DataSet datasetInts = HDFfile.getDataSet(IntsNameSet);
            HighFive::DataSet datasetFloat = HDFfile.getDataSet(FloatNameSet);
            HighFive::DataSet datasetPor = HDFfile.getDataSet(PORNameSet);
            std::vector<std::string> names;
            std::vector<std::vector<int>> IDS;
            std::vector<std::vector<double>> DATA;
            std::vector<int> POR;
            datasetNames.read(names);
            datasetInts.read(IDS);
            datasetFloat.read(DATA);
            datasetPor.read(POR);

            if (IDS[0].size() != DATA[0].size()){
                std::cout << "The rows of integer and float data do not match" << std::endl;
                return false;
            }

            bool tf1 = addPorosityScenarios(names[0], POR);
            if (!tf1){
                std::cout << "Error while inserting the porosity indices for flow scenario [ "
                          << names[0] << " ]" << std::endl;
                return false;
            }

            int eid, sid, r, c, riv, npxl;
            double w, len;
            bool tf;
            for (unsigned int i = 0; i < IDS[0].size(); ++i){
                eid = IDS[0][i];
                sid = IDS[1][i];
                r = IDS[2][i];
                c = IDS[3][i];
                riv = IDS[4][i];
                npxl = IDS[5][i];
                w = DATA[0][i];
                len = DATA[1][i];
                int idx = 2;
                std::vector<double> M,S,A;
                for (unsigned int j = 0; j < POR.size(); ++j){
                    M.push_back(DATA[idx][i]);
                    S.push_back(DATA[idx+1][i]);
                    A.push_back(DATA[idx+2][i]);
                    idx = idx + 3;
                }
                tf = addStreamline(names[0], eid, sid, r, c, riv, npxl, w, len, M, S, A);
                if (!tf){
                    std::cout << "Error while inserting the streamline for flow scenario [ "
                              << names[0] << " ] with Eid: " << IDS[0][i] << ", Sid: " << IDS[1][i] << std::endl;
                    return false;
                }
                // Debug
                //if (i == 1000){
                //    std::map<std::string ,WellList>::iterator flowit;
                //    flowit = FlowScenarios.find(names[0]);
                //    std::map<int, Well>::iterator wellit;
                //    wellit = flowit->second.Wells.find(174);
                //    wellit->second.log();
                //}
            }
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "Read URFS in " << elapsed.count() << std::endl;
            return true;
        }
#endif
        std::ifstream ifile1, ifile2;
        std::string filename1 = filename + "_data.mnts";
        std::string filename2 = filename + "_msa.mnts";
        ifile1.open(filename1);
        bool bfilesGood = true;
        // Test the first file
        if (!ifile1.is_open()){
            std::cout << "Cant open file: " << filename1 << std::endl;
            bfilesGood = false;
        }
        if (!bfilesGood){
            return false;
        }
        // Test the second file
        ifile2.open(filename2);
        if (!ifile2.is_open()){
            std::cout << "Cant open file: " << filename2 << std::endl;
            bfilesGood = false;
        }
        if (!bfilesGood){
            return false;
        }

        std::cout << "Reading " << filename1 << " and " << filename2 << std::endl;
        std::string line1, line2;
        int Nurfs, Npor;
        URFTYPE urftype;
        std::string name, type;
        {// Read the header from the first file
            getline(ifile1, line1);
            std::istringstream inp(line1.c_str());
            inp >> Nurfs;
            inp >> name;
            inp >> Npor;
        }
        {// Read porosity data from the first file
            bool tf;
            std::vector<int> por;
            getline(ifile1, line1);
            std::istringstream inp(line1.c_str());
            int p;
            for (int i = 0; i < Npor; ++i){
                inp >> p;
                por.push_back(p);
            }
            tf = addPorosityScenarios(name, por);
            if (!tf){
                std::cout << "Error while inserting the porosity indices for flow scenario [ "
                          << name << " ]" << std::endl;
                return false;
            }
        }

        {// Read the data from both files
            int eid, sid, r, c, riv, npxl;
            double m, s, w, age, len;
            bool tf;
            for (int i = 0; i < Nurfs; ++i){
                getline(ifile1, line1);
                getline(ifile2, line2);
                std::istringstream inp1(line1.c_str());
                inp1 >> eid;
                inp1 >> sid;
                inp1 >> r;
                inp1 >> c;
                inp1 >> riv;
                inp1 >> w;
                inp1 >> len;
                inp1 >> npxl;
                std::vector<double> M, S, A;
                std::istringstream inp2(line2.c_str());
                for (int j = 0; j < Npor; ++j){
                    inp2 >> m;
                    inp2 >> s;
                    inp2 >> age;
                    M.push_back(m);
                    S.push_back(s);
                    A.push_back(age);
                }


                tf = addStreamline(name, eid, sid, r, c, riv, npxl, w, len, M, S, A);
                if (!tf){
                    std::cout << "Error while inserting the streamline for flow scenario [ "
                              << name << " ] with Eid: " << eid << ", Sid: " << sid << std::endl;
                    return false;
                }
                // Debug
                //if (i == 1000){
                //    std::map<std::string ,WellList>::iterator flowit;
                //    flowit = FlowScenarios.find(name);
                //    std::map<int, Well>::iterator wellit;
                //    wellit = flowit->second.Wells.find(174);
                //    wellit->second.log();
                //}
            }
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read URFS in " << elapsed.count() << std::endl;
        return true;
    }
}

#endif //MANTISSERVER_WELLS_H
