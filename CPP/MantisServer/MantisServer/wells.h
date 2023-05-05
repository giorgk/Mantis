//
// Created by giorg on 2/25/2023.
//

#ifndef MANTISSERVER_WELLS_H
#define MANTISSERVER_WELLS_H

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

        std::vector<cell> SourceArea;
    };

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

        void calculateSourceArea(BackroundRaster &braster, RechargeScenario &rch, bool doCalc);
        void calculateWeights();

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

    void Well::calculateSourceArea(BackroundRaster &braster, RechargeScenario &rch, bool doCalc) {
        bool tf;
        int lin_ind, rs, cs, rn, cn;
        double Xorig, Yorig, cellSize, xs, ys, dNr, side1, side2;
        double SAxmin, SAxmax, SAymin, SAymax;
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


        for (itstrml = streamlines.begin(); itstrml != streamlines.end(); ++itstrml){
            lin_ind = braster.IJ(itstrml->second.row, itstrml->second.col);
            if (!doCalc){
                if (lin_ind != -1){
                    itstrml->second.SourceArea.push_back(cell(itstrml->second.row, itstrml->second.col));
                    continue;
                }
            }
            else{
                rs = itstrml->second.row;
                cs = itstrml->second.col;
                braster.cellCoords(rs, cs, xs, ys);
                //xs = Xorig + cellSize/2 + cellSize*(cs);
                // For the Y the row numbers start from the top
                //ys =  Yorig + cellSize*dNr - cellSize/2 - cellSize*(rs);
                side1 = cellSize + cellSize * itstrml->second.Npxl;
                side2 = std::max(2.0*cellSize, side1/ratio/2.0);
                SAxmin = -side2 - 25;
                SAxmax =  side2 + 25;
                SAymin = -side1/2 - 25;
                SAymax =  side1/2 + 25;

                x1 = (cosd * SAxmin + sind * SAymin) + xs;
                y1 = (-sind * SAxmin + cosd * SAymin) + ys;

                x2 = (cosd * SAxmin + sind * SAymax) + xs;
                y2 = (-sind * SAxmin + cosd * SAymax) + ys;

                x3 = (cosd * SAxmax + sind * SAymax) + xs;
                y3 = (-sind * SAxmax + cosd * SAymax) + ys;

                x4 = (cosd * SAxmax + sind * SAymin) + xs;
                y4 = (-sind * SAxmax + cosd * SAymin) + ys;

                std::cout << "pp=[" << x1 << " " << y1 << "; " << x2 << " " << y2 << "; "
                                    << x3 << " " << y3 << "; " << x4 << " " << y4 << "]; " << std::endl;

                std::map< int, cell> for_test;
                std::map< int, cell> tested;
                std::map< int, cell> next_round;

                lin_ind = Nrow * cs + rs;
                for_test.insert(std::pair<int, cell>(lin_ind, cell(rs,cs)));
                if (braster.IJ(rs,cs) == -1){
                    // The first pixel is outside the study area.
                    // We will assume zero loading
                    itstrml->second.mu = 0.0;
                    itstrml->second.std = 0.0;
                    continue;
                }

                Qtmp = 0.0;
                Qtarget = pumpingRate * itstrml->second.w;

                bool sourceFound = false;
                while (true){
                    for (itcell1 = for_test.begin(); itcell1 != for_test.end(); ++itcell1){
                        tested.insert(std::pair<int,cell>(itcell1->first, itcell1->second));
                        braster.cellCoords(itcell1->second.row, itcell1->second.col, xtmp, ytmp);
                        std::cout << "plot(" << xtmp << "," << ytmp << ",'.k');" << std::endl;
                        //xtmp = -223275.0 + 50.0*(it4->second.col);
                        //ytmp =  298525.0 - 50.0*(it4->second.row);
                        // Check if the point is in the source area by testing the barycentric coordinates
                        tf = isInTriangle(x1, y1, x2, y2, x3, y3, xtmp, ytmp);
                        if (!tf){
                            tf = isInTriangle(x1, y1, x3, y3, x4, y4, xtmp, ytmp);
                            if (!tf){
                                continue;
                            }
                        }
                        lin_ind = braster.IJ(itcell1->second.row, itcell1->second.col);
                        if (lin_ind != -1){
                            tf = rch.getValues(lin_ind, rch_v, dummy);
                            if (tf){
                                if (rch_v > 10.0){
                                    itstrml->second.SourceArea.push_back(itcell1->second);
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
                                    if (braster.IJ(rn,cn) == -1)
                                        continue;
                                    lin_ind = Nrow * cn + rn;
                                    itcell2 = tested.find(lin_ind);
                                    if (itcell2 == tested.end()){
                                        for_test.insert(std::pair<int,cell>(lin_ind, cell(rn,cn)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    class WellList{
    public:
        WellList(){}
        void addWell(int Eid, Well w);
        bool addstreamline(int Eid,int Sid, int row_ind, int col_ind, double w,
                           int npxl, URFTYPE type, int riv,
                           double paramA, double paramB,
                           double paramC = 0, double paramD = 0);
        bool calcSourceArea = false;
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
        void calcWellSourceArea(BackroundRaster &braster, RechargeScenarioList &rchList);
        std::map<std::string ,WellList> FlowScenarios;
    private:
        bool readWells(std::string filename, BMapCollection &Bmaps);
        bool readURFs(std::string filename);
        bool addStreamline(std::string flowScenName, int wellid,
                           int Sid, int row_ind, int col_ind, double w,
                           int npxl, URFTYPE type, int riv,
                           double paramA, double paramB,
                           double paramC = 0, double paramD = 0);
    };

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

    void FlowWellCollection::calcWellSourceArea(BackroundRaster &braster, RechargeScenarioList &rchList) {
        std::map<std::string ,WellList>::iterator flowit;
        std::map<int, Well>::iterator wellit;
        std::map<std::string, RechargeScenario>::iterator rchit;
        for (flowit = FlowScenarios.begin(); flowit != FlowScenarios.end(); ++flowit){
            rchit = rchList.RechargeList.find(flowit->second.rch_map);
            if (rchit == rchList.RechargeList.end()){
                std::cout << "I can't find a recharge scenario with name " << flowit->second.rch_map << std::endl;
                continue;
            }
            for (wellit = flowit->second.Wells.begin(); wellit != flowit->second.Wells.end(); ++ wellit){
                wellit->second.calculateSourceArea(braster, rchit->second, flowit->second.calcSourceArea);
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
                inp >> calcSource;
                inp >> rch_scen;
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
