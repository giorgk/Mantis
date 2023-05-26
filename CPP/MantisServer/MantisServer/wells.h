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
    /**
	 * @brief Stores data for each streamline.
	 *
	 * Although this is a class, it is used more like struct container.
	 *
	 */
    class Streamline{
    public:
        Streamline(){};
        /*! \brief streamlineClass constructor expects the parameters that define a streamline
		\param row_ind is the index in GNLM loading where this streamline starts from near the land surface.
		\param col_ind is the index in SWAT loading where this streamline starts from near the land surface.
		\param w_in is the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
		\param rch_in is the groundwater recharge rate in m/day according to the flow model
		\param type_in is the type of the unit response function
		\param paramA this is either the mean value or the streamline length.
		\param paramB this is either the standard deviation or the velocity.
		\param paramC if the type is both this is mean, while A and B are length and velocity
		\param paramD if the type is both this is standard deviation
		*/
        void setParameters(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
                        double paramA, double paramB, double paramC = 0, double paramD = 0);
        //void addSourceAreaCell(int lin_ind, int row, int col);
        void addSourceAreaCell(cell c);
        void clearSourceArea();

        void print();

        //! the row number of the pixel where this streamline starts from near the land surface.
        int row;

        //! the column number of the pixel where this streamline starts from near the land surface.
        int col;

        //! the mean value of the fitted unit response function.
        double mu;

        //! the standard deviation of the fitted unit response function.
        double std;

        //!the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
        double w;

        //! The streamline length
        double sl;

        //! the mean velocity along the streamline
        double vel;

        //! This is true for the streamlines that start from rivers
        bool inRiver;

        //! The type of streamline
        URFTYPE type;

        int Npxl;

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

    void Streamline::print() {
        std::cout << "r: " << row << ", c: " << col << ", m: " << mu << ", s: " << std << ", w: " << w << ", riv: " << inRiver << ", N: " << Npxl << std::endl;
    }

    void Streamline::setParameters(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
                                     double paramA, double paramB, double paramC, double paramD) {
        row = row_ind;
        col = col_ind;
        w = w_in;
        Npxl = npxl_in;
        type = type_in;
        inRiver = Riv == 1;
        switch (type)
        {
            case URFTYPE::LGNRM:
                mu = paramA;
                std = paramB;
                sl = paramC;
                vel = paramD;
                break;
            case URFTYPE::ADE:
                sl = paramA;
                vel = paramB;
                mu = paramC;
                std = paramD;
                break;
            case URFTYPE::BOTH:
                sl = paramA;
                vel = paramB;
                mu = paramC;
                std = paramD;
                break;
            default:
                break;
        }
    }

    class Well{
    public:
        void addStreamline(int Sid, int row_ind, int col_ind, double w, int npxl, URFTYPE type, int riv,
                           double paramA, double paramB, double paramC = 0, double paramD = 0);
        void setAdditionalData(double x, double y, double d, double s, double q, double r, double a);

        void calculateSourceArea(BackgroundRaster &braster, RechargeScenario &rch, int doCalc, bool debug = false);
        void calculateSourceAreaType0(BackgroundRaster &braster, bool debug = false);
        void calculateSourceAreaType1(BackgroundRaster &braster, RechargeScenario &rch, bool debug = false);
        void calculateSourceAreaType2(BackgroundRaster &braster, RechargeScenario &rch, bool debug = false);
        void calculateWeights();
        bool bsimulateThis(Scenario &scenario);

        std::map<int, Streamline> streamlines;
        double xcoord;
        double ycoord;
        double depth;
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

    void Well::addStreamline(int Sid, int row_ind, int col_ind, double w, int npxl, URFTYPE type, int riv,
                                  double paramA, double paramB, double paramC, double paramD) {
        Streamline s;
        s.setParameters(row_ind, col_ind, w, npxl, type, riv, paramA, paramB, paramC, paramD);
        streamlines.insert(std::pair<int, Streamline>(Sid, s));
    }

    void Well::setAdditionalData(double x, double y, double d, double s, double q, double r, double a){
        xcoord = x;
        ycoord = y;
        depth = d;
        screenLength = s;
        pumpingRate = q;
        ratio = r;
        angle = a;
    }

    void Well::calculateSourceArea(BackgroundRaster &braster, RechargeScenario &rch, int doCalc, bool debug){
        if (doCalc == 0){
            calculateSourceAreaType0(braster, debug);
        }
        else if (doCalc == 1){
            calculateSourceAreaType1(braster, rch, debug);
        }
        else if (doCalc == 2){
            calculateSourceAreaType2(braster, rch, debug);
        }
    }

    void Well::calculateSourceAreaType0(BackgroundRaster &braster, bool debug){
        std::map<int, Streamline>::iterator itstrml;
        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            int rasterValue = braster.IJ(itstrml->second.row, itstrml->second.col);
            if (rasterValue == -1){
                // The source area is outside the active area therefore the loading will be zero
                itstrml->second.mu = 0.0;
                itstrml->second.std = 0.0;
                continue;
            }
            if (itstrml->second.inRiver){
                itstrml->second.mu = 0.0;
                itstrml->second.std = 0.0;
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
        }

    }

    void Well::calculateSourceAreaType1(BackgroundRaster &braster, RechargeScenario &rch, bool debug) {
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
                itstrml->second.mu = 0.0;
                itstrml->second.std = 0.0;
                continue;
            }

            if (itstrml->second.inRiver){
                itstrml->second.mu = 0.0;
                itstrml->second.std = 0.0;
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
                niter++;
                itstrml->second.clearSourceArea();
            }

        }
    }

    void Well::calculateSourceAreaType2(BackgroundRaster &braster, RechargeScenario &rch, bool debug){
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
        std::vector<boost_point> srcPnts;
        boost_poly srcPoly;


        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            if (itstrml->second.inRiver){
                itstrml->second.mu = 0.0;
                itstrml->second.std = 0.0;
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
            if (Npxl == 1){
                Npxl = 0;
            }
            while (count_iter < 20){
                double lngth = cellSize + cellSize * Npxl;
                double wdth = std::max(lngth/ratio/1.95, cellSize/1.95);
                double SAxmin = -wdth;
                double SAxmax =  wdth;
                double SAymin = -lngth/1.95 - cellSize/1.95;
                double SAymax =  lngth/1.95 + cellSize/1.95;
                srcX.clear(); srcY.clear(); srcXrt.clear(), srcYrt.clear(); srcPnts.clear();
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
                    srcPnts.push_back(boost_point(srcXrt[i], srcYrt[i]));

                }
                boost::geometry::assign_points(srcPoly, srcPnts);
                boost::geometry::correct(srcPoly);

                //std::cout << "srcX = [" << srcXrt[0] << "," << srcXrt[1] << "," << srcXrt[2] << "," << srcXrt[3] << "]" << std::endl;
                //std::cout << "scrY = [" << srcYrt[0] << "," << srcYrt[1] << "," << srcYrt[2] << "," << srcYrt[3] << "]" << std::endl;

                int dJ = static_cast<int>(std::ceil((Xcntr - minSrcXrt)/cellSize));
                int dI = static_cast<int>(std::ceil((Ycntr - minSrcYrt)/cellSize));
                int JgridStart = std::max(1, col - dJ);
                int JgridEnd = std::min(Ncol, col + dJ);
                int IgridStart = std::max(1, row - dI);
                int IgridEnd = std::min(Ncol, row + dI);

                //Debug only
                int cnt_cells = 0;
                bool done_that = false;

                for (int i = IgridStart; i <= IgridEnd; ++i){
                    for (int j = JgridStart; j <= JgridEnd; ++j){
                        braster.cellCoords(i, j, xcell, ycell);
                        //std::cout << xcell << "," << ycell << std::endl;
                        bool tf = boost::geometry::within(boost_point(xcell, ycell), srcPoly);

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

                                    // Remove these lines
                                    cnt_cells++;
                                    if (cnt_cells > 3){
                                        done_that = true;
                                    }
                                }
                            }
                        }
                        // Remove these lines
                        if (done_that)
                            break;
                    }
                    // Remove these lines
                    if (done_that)
                        break;
                }
                if (Qtmp >= Qtarget){
                    break;
                }
                Npxl = Npxl + 1.0;
                count_iter++;
            }
            std::sort(itstrml->second.SourceArea.begin(), itstrml->second.SourceArea.end(), compareCellByDistance);
        }
    }

    bool Well::bsimulateThis(Scenario &scenario) {
        if (scenario.bNarrowSelection == true){
            if (scenario.useRadSelect){
                if (!scenario.RadSelect.isPointIn(xcoord, ycoord)){
                    return false;
                }
                if (scenario.useRectSelect){
                    if (!scenario.RectSelect.isPointIn(xcoord, ycoord)){
                        return false;
                    }
                }
                if (scenario.useDepthRange){
                    if (!scenario.DepthRange.isInRange(depth)){
                        return false;
                    }
                }
                if (scenario.useScreenLenghtRange){
                    if (!scenario.ScreenLengthRange.isInRange(screenLength)){
                        return false;
                    }
                }
            }
        }
        return true;
    }

    class WellList{
    public:
        WellList(){}
        void addWell(int Eid, Well w);
        bool addstreamline(int Eid,int Sid, int row_ind, int col_ind, double w,
                           int npxl, URFTYPE type, int riv,
                           double paramA, double paramB,
                           double paramC = 0, double paramD = 0);
        int calcSourceArea = false;
        std::string rch_map;
        std::map<int, Well> Wells;
    };

    bool WellList::addstreamline(int Eid,int Sid, int row_ind, int col_ind, double w,
                                 int npxl, URFTYPE type, int riv,
                                 double paramA, double paramB,
                                 double paramC, double paramD){
        std::map<int, Well>::iterator eidit;
        eidit = Wells.find(Eid);
        if (eidit == Wells.end()){
            std::cout << "I can't find a well with id [ " << Eid << " ]" << "in the list" << std::endl;
            return false;
        }
        eidit->second.addStreamline(Sid,row_ind,col_ind,w,npxl,type,riv,paramA,paramB,paramC,paramD);
        return true;
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
        bool hasFlowScenario(std::string &flowScen);
    private:
        bool readWells(std::string filename, BMapCollection &Bmaps);
        bool readURFs(std::string filename);
        bool addStreamline(std::string flowScenName, int wellid,
                           int Sid, int row_ind, int col_ind, double w,
                           int npxl, URFTYPE type, int riv,
                           double paramA, double paramB,
                           double paramC = 0, double paramD = 0);
    };

    bool FlowWellCollection::hasFlowScenario(std::string &flowScen) {
        std::map<std::string ,WellList>::iterator flowit;
        flowit = FlowScenarios.find(flowScen);
        if (flowit != FlowScenarios.end()){
            return true;
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
            for (wellit = flowit->second.Wells.begin(); wellit != flowit->second.Wells.end(); ++ wellit){
                //std::cout << wellit->first << std::endl;
                //if (wellit->first == 2893){
                //    dbg = true;
                //    std::cout << "Stop here" << std::endl;
                //}
                wellit->second.calculateSourceArea(braster, rchit->second, flowit->second.calcSourceArea, dbg);
                //dbg = false;
                count++;
                if (count % 1000 == 0){
                    std::cout << "----" << count << "----" << std::endl;
                }
            }
        }


        //flowit->second.calcSourceArea
    }

    bool FlowWellCollection::addStreamline(std::string flowScenName, int wellid,
                                           int Sid, int row_ind, int col_ind, double w,
                                           int npxl, URFTYPE type, int riv,
                                           double paramA, double paramB,
                                           double paramC, double paramD){
        std::map<std::string ,WellList>::iterator flowit;
        flowit = FlowScenarios.find(flowScenName);
        if (flowit == FlowScenarios.end()){
            std::cout << "I can't find wells under the [ " << flowScenName << " ] flow scenario" << std::endl;
            return false;
        }
        flowit->second.addstreamline(wellid,Sid,row_ind,col_ind,w,npxl,type,riv,paramA,paramB,paramC,paramD);
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
            int Nwells, Eid, calcSource;
            std::string setName, line, rch_scen;
            WellList wList;
            {
                getline(Welldatafile, line);
                std::istringstream inp(line.c_str());
                inp >> Nwells;
                inp >> setName;
                inp >> rch_scen;
                inp >> calcSource;
                wList.calcSourceArea = calcSource;
                wList.rch_map = rch_scen;
            }

            double xw, yw, D, SL, Q, ratio, angle;
            std::string regionCode;
            bool tf;
            for (int i = 0; i < Nwells; ++i){
                getline(Welldatafile, line);
                std::istringstream inp(line.c_str());
                inp >> Eid;
                //std::cout << Eid << std::endl;
                inp >> xw;
                inp >> yw;
                inp >> D;
                inp >> SL;
                inp >> Q;
                inp >> ratio;
                inp >> angle;
                Well w;
                w.setAdditionalData(xw,yw,D,SL,Q,ratio,angle);
                wList.addWell(Eid, w);
                for (int j = 0; j < Bmaps.Nbmaps(); ++j){
                    inp >> regionCode;
                    tf = Bmaps.addwell(regionCode,setName,Eid);
                    if (!tf){
                        std::cout << "I cannot find a region with name [ " << regionCode << " ]" << std::endl;
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
            const std::string FloatNameSet("MSW");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
            HighFive::DataSet datasetInts = HDFfile.getDataSet(IntsNameSet);
            HighFive::DataSet datasetFloat = HDFfile.getDataSet(FloatNameSet);
            std::vector<std::string> names;
            std::vector<std::vector<int>> IDS;
            std::vector<std::vector<double>> DATA;
            datasetNames.read(names);
            datasetInts.read(IDS);
            datasetFloat.read(DATA);
            if (names.size() != 2){
                std::cout << "2 names are needed for the URFset. " << names.size() << " provided" << std::endl;
                return false;
            }
            if (IDS[0].size() != DATA[0].size()){
                std::cout << "The rows of integer and float data do not match" << std::endl;
                return false;
            }
            if (IDS.size() != 6 || DATA.size() != 3){
                std::cout << "Incorrect number of columns" << std::endl;
                std::cout << "The size of integers must be 6 and for the floats 3" << std::endl;
                return false;
            }

            URFTYPE urftype;
            if (names[1].compare("LGNRM") == 0)
                urftype = URFTYPE::LGNRM;
            else if (names[1].compare("ADE") == 0)
                urftype = URFTYPE::ADE;
            else if (names[1].compare("BOTH") == 0)
                urftype = URFTYPE::BOTH;
            else{
                std::cout << "The URF type " << names[1] << " is not valid" << std::endl;
                return false;
            }

            int eid, sid, r, c, riv, npxl;
            double m, s, w;
            bool tf;
            for (unsigned int i = 0; i < IDS[0].size(); ++i){
                tf = addStreamline(names[0], IDS[0][i], IDS[1][i], IDS[2][i], IDS[3][i],
                              DATA[2][i], IDS[5][i], urftype, IDS[4][i],
                              DATA[0][i], DATA[1][i]);
                if (!tf){
                    std::cout << "Error while inserting the streamline for flow scenario [ "
                              << names[0] << " ] with Eid: " << IDS[0][i] << ", Sid: " << IDS[1][i] << std::endl;
                    return false;
                }
            }
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "Read URFS in " << elapsed.count() << std::endl;
            return true;
        }
#endif
        std::ifstream ifile;
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            int Nurfs;
            URFTYPE urftype;
            std::string name, type;
            {// Read the header
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> Nurfs;
                inp >> name;
                inp >> type;

                if (type.compare("LGNRM") == 0)
                    urftype = URFTYPE::LGNRM;
                else if (type.compare("ADE") == 0)
                    urftype = URFTYPE::ADE;
                else if (type.compare("BOTH") == 0)
                    urftype = URFTYPE::BOTH;
                else{
                    std::cout << "The URF type " << type << " is not valid" << std::endl;
                    return false;
                }
            }

            {// Read the data
                int eid, sid, r, c, riv, npxl;
                double m, s, w;
                bool tf;
                for (int i = 0; i < Nurfs; ++i){
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> eid;
                    inp >> sid;
                    inp >> r;
                    inp >> c;
                    inp >> riv;
                    inp >> m;
                    inp >> s;
                    inp >> w;
                    inp >> npxl;
                    tf = addStreamline(name, eid, sid, r, c, w, npxl, urftype, riv, m, s);
                    if (!tf){
                        std::cout << "Error while inserting the streamline for flow scenario [ "
                                  << name << " ] with Eid: " << eid << ", Sid: " << sid << std::endl;
                        return false;
                    }
                }
            }
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read URFS in " << elapsed.count() << std::endl;
        return true;
    }
}

#endif //MANTISSERVER_WELLS_H
