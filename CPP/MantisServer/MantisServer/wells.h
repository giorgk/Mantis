//
// Created by giorg on 2/25/2023.
//

#ifndef MANTISSERVER_WELLS_H
#define MANTISSERVER_WELLS_H

# include "MShelper.h"

namespace mantisServer{
    /**
	 * @brief Stores data for each streamline.
	 *
	 * Although this is a class, it is used more like struct container.
	 *
	 */
    class Streamline{
    public:
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
        Streamline(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
                        double paramA, double paramB, double paramC = 0, double paramD = 0);

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

    Streamline::Streamline(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
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
        streamlines.insert(std::pair<int, Streamline>(Sid, Streamline( row_ind, col_ind, w, npxl, type, riv, paramA, paramB, paramC, paramD)));
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

        std::map<int, Well> Wells;
    };

    void WellList::addWell(int Eid, Well w) {
        Wells.insert(std::pair<int, Well>(Eid, w));
    }

    class FlowWellCollection{
    public:
        FlowWellCollection(){}

        bool readMainfile(std::string path, std::string filename, bool bReadWells);

        std::map<std::string ,WellList> FlowScenarios;
    private:
        bool readWells(std::string filename);
        bool readURFs(std::string filename);
    };

    bool FlowWellCollection::readMainfile(std::string path, std::string filename, bool bReadWells) {
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
                    tf = readWells(filename1);
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

    bool FlowWellCollection::readWells(std::string filename) {
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
            getline(Welldatafile, line);
            std::istringstream inp(line.c_str());
        }



        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Wells in " << elapsed.count() << std::endl;
        return true;
    }

    bool FlowWellCollection::readURFs(std::string filename) {
        auto start = std::chrono::high_resolution_clock::now();


        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read URFS in " << elapsed.count() << std::endl;
        return false;
    }
}

#endif //MANTISSERVER_WELLS_H
