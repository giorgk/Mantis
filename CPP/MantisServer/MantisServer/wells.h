//
// Created by giorg on 2/25/2023.
//

#ifndef MANTISSERVER_WELLS_H
#define MANTISSERVER_WELLS_H

# include "MShelper.h"
#include "BMaps.h"

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

        std::map<int, Streamline> streamlines;
        double xcoord;
        double ycoord;
        double depth;
        double screenLength;
        double pumpingRate;
        double ratio;
        double angle;
    };

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
