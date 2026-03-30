//
// Created by giorg on 10/22/2024.
//

#ifndef MANTISSA_MS_UNSAT_H
#define MANTISSA_MS_UNSAT_H
namespace MS{

    class UNSAT{
    public:
        UNSAT()
        : wc(0.0), minD(0.0), minR(0.0), Ncells(0)
        {}

        bool readdata(const std::string& dpth_file, const std::string& dpth_name, const std::string& rch_file,
                      double wc_in, double d_in, double r_in, boost::mpi::communicator& world);

        int traveltime(int IJ) const;
        double getDepth(int IJ) const;
        double getRch(int IJ) const;
        double getSurfPerc(int IJ) const;
    private:
        double wc;
        double minD;
        double minR;
        std::vector<double> Depth;
        std::vector<double> Rch;
        std::vector<double> SurfPerc;
        int Ncells;

        bool read_depth(const std::string& dpth_file, const std::string& dpth_name,
            boost::mpi::communicator& world);
        bool read_rch(const std::string& rch_file, boost::mpi::communicator& world);
    };

    inline bool UNSAT::readdata(const std::string& dpth_file, const std::string& dpth_name, const std::string& rch_file,
                                double wc_in, double d_in, double r_in,
                                boost::mpi::communicator& world) {
        // Reset state first, so a failed reread does not leave stale data behind
        wc = wc_in;
        minD = d_in;
        minR = r_in;
        Ncells = 0;
        Depth.clear();
        Rch.clear();
        SurfPerc.clear();

        // -----------------------------
        // Read Depth
        // -----------------------------
        bool tf = read_depth(dpth_file, dpth_name, world);
        if (!tf)
        {
            return false;
        }


        // -----------------------------
        // Read Recharge + SurfPerc
        // -----------------------------
        tf = read_rch(rch_file, world);
        if (!tf)
        {
            return false;
        }

        return true;
    }

    inline int UNSAT::traveltime(int IJ) const {
        if (IJ < 0 || IJ >= Ncells || Ncells <= 0){
            return 0;
        }
        double d = Depth[IJ];
        if (d < minD){
            d = minD;
        }

        double r = Rch[IJ];
        if (r < minR){
            r = minR;
        }

        // d : meters
        // r : mm/year
        // r/1000 : m/year
        // tau : years
        const double tau = wc * d / (r/1000.0);

        if (tau <= 0.0)
        {
            return 0;
        }

        if (tau > static_cast<double>(std::numeric_limits<int>::max()))
        {
            return std::numeric_limits<int>::max();
        }

        return static_cast<int>(std::ceil(tau));
    }

    inline double UNSAT::getDepth(int IJ) const {
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return Depth[IJ];
    }

    inline double UNSAT::getRch(int IJ) const {
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return Rch[IJ];
    }

    inline double UNSAT::getSurfPerc(int IJ) const{
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return SurfPerc[IJ];
    }

    inline bool UNSAT::read_depth(const std::string& dpth_file, const std::string& dpth_name,
                                    boost::mpi::communicator& world) {
        int readSuccess = 0;
        int localNcells = 0;

        Depth.clear();

        if (world.rank() == 0) {
            std::vector<std::string> names;
            std::vector<std::vector<double>> data_NxNames;

            try {
                const std::string ext = getExtension(dpth_file);

                if (ext == "h5" || ext == "H5") {
#if _USEHF>0
                    std::cout << "Reading HDF5 file " << dpth_file << std::endl;

                    const std::string NamesNameSet("Names");
                    const std::string DataNameSet("Data");

                    HighFive::File HDFfile(dpth_file, HighFive::File::ReadOnly);
                    HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
                    HighFive::DataSet datasetData  = HDFfile.getDataSet(DataNameSet);

                    std::vector<std::vector<double>> tmp;
                    datasetNames.read(names);
                    datasetData.read(tmp);

                    if (names.empty())
                    {
                        std::cout << "Error: dataset 'Names' is empty in " << dpth_file << std::endl;
                    }
                    else if (tmp.empty())
                    {
                        std::cout << "Error: dataset 'Data' is empty in " << dpth_file << std::endl;
                    }
                    else if (tmp.size() == names.size())
                    {
                        // names x N  -> transpose to N x names
                        const std::size_t nNames = names.size();
                        const std::size_t N = tmp[0].size();

                        data_NxNames.assign(N, std::vector<double>(nNames, 0.0));
                        for (std::size_t j = 0; j < nNames; ++j)
                        {
                            for (std::size_t i = 0; i < N; ++i)
                            {
                                data_NxNames[i][j] = tmp[j][i];
                            }
                        }
                    }
                    else if (tmp[0].size() == names.size())
                    {
                        // already N x names
                        data_NxNames = tmp;
                    }
                    else
                    {
                        std::cout << "Error: cannot interpret dimensions in " << dpth_file
                                  << ". Names.size() = " << names.size()
                                  << ", Data = (" << tmp.size()
                                  << " x " << (tmp.empty() ? 0 : tmp[0].size()) << ")"
                                  << std::endl;
                    }
#else
                    std::cout << "Error: input file " << dpth_file
                          << " is HDF5, but this executable was built without HDF5 support."
                          << std::endl;
#endif

                }
                else {
                    std::vector<std::vector<double>> tmp;
                    const bool tf = readLinearData(dpth_file, names, data_NxNames, 5000000);

                    if (!tf)
                    {
                        std::cout << "Error: failed to read depth file " << dpth_file << std::endl;
                    }
                }

                if (!data_NxNames.empty()) {
                    if (data_NxNames[0].size() != names.size())
                    {
                        std::cout << "Error: normalized depth matrix in " << dpth_file
                                  << " does not have names.size() columns." << std::endl;
                    }
                    else
                    {
                        int idx = -1;
                        for (int i = 0; i < static_cast<int>(names.size()); ++i)
                        {
                            if (names[i] == dpth_name)
                            {
                                idx = i;
                                break;
                            }
                        }

                        if (idx < 0)
                        {
                            std::cout << "Error: depth scenario '" << dpth_name
                                      << "' was not found in " << dpth_file << std::endl;
                        }
                        else
                        {
                            localNcells = static_cast<int>(data_NxNames.size());
                            Depth.resize(localNcells, 0.0);

                            for (int i = 0; i < localNcells; ++i)
                            {
                                Depth[i] = data_NxNames[i][idx];
                            }

                            readSuccess = 1;
                        }
                    }
                }
            }
            catch (const std::exception& e)
            {
                std::cout << "Error reading depth file " << dpth_file
                          << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
        }

        world.barrier();
        sendScalarFromRoot2AllProc<int>(readSuccess, world);

        if (readSuccess == 0){
            return false;
        }

        sendScalarFromRoot2AllProc<int>(localNcells, world);
        Ncells = localNcells;

        if (world.rank() > 0)
        {
            Depth.resize(Ncells, 0.0);
        }

        world.barrier();
        sendVectorFromRoot2AllProc<double>(Depth, world);

        if (PrintMatrices){
            printVectorForAllProc<double>(Depth,world, 0, 30);
        }

        return true;
    }

    inline bool UNSAT::read_rch(const std::string& rch_file,
                            boost::mpi::communicator& world)
    {
        int readSuccess = 0;

        Rch.clear();
        SurfPerc.clear();

        if (world.rank() == 0)
        {
            try
            {
                const std::string ext = getExtension(rch_file);

                if (ext == "h5" || ext == "H5")
                {
    #if _USEHF > 0
                    std::cout << "Reading HDF5 file " << rch_file << std::endl;

                    const std::string NamesNameSet("Names");
                    const std::string DataNameSet("Data");

                    HighFive::File HDFfile(rch_file, HighFive::File::ReadOnly);
                    HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
                    HighFive::DataSet datasetData  = HDFfile.getDataSet(DataNameSet);

                    std::vector<std::string> names;
                    std::vector<std::vector<double>> tmp;

                    datasetNames.read(names);
                    datasetData.read(tmp);

                    if (tmp.empty())
                    {
                        std::cout << "Error: recharge file " << rch_file
                                  << " contains empty dataset 'Data'." << std::endl;
                    }
                    else if (tmp.size() == 2)
                    {
                        // HDF5 case 1: 2 x N
                        Rch = tmp[0];
                        SurfPerc = tmp[1];

                        if (static_cast<int>(Rch.size()) != Ncells)
                        {
                            std::cout << "Error: size mismatch between Depth and Recharge in "
                                      << rch_file << ". Depth size = " << Ncells
                                      << ", Recharge size = " << Rch.size() << std::endl;
                        }
                        else if (static_cast<int>(SurfPerc.size()) != Ncells)
                        {
                            std::cout << "Error: size mismatch between Depth and SurfPerc in "
                                      << rch_file << ". Depth size = " << Ncells
                                      << ", SurfPerc size = " << SurfPerc.size() << std::endl;
                        }
                        else
                        {
                            readSuccess = 1;
                        }
                    }
                    else if (!tmp.empty() && tmp[0].size() == 2)
                    {
                        // HDF5 case 2: N x 2
                        const int N = static_cast<int>(tmp.size());

                        if (N != Ncells)
                        {
                            std::cout << "Error: size mismatch between Depth and Recharge in "
                                      << rch_file << ". Depth size = " << Ncells
                                      << ", Recharge rows = " << N << std::endl;
                        }
                        else
                        {
                            Rch.resize(Ncells, 0.0);
                            SurfPerc.resize(Ncells, 0.0);

                            for (int i = 0; i < Ncells; ++i)
                            {
                                Rch[i]      = tmp[i][0];
                                SurfPerc[i] = tmp[i][1];
                            }

                            readSuccess = 1;
                        }
                    }
                    else
                    {
                        std::cout << "Error: cannot interpret dimensions in " << rch_file
                                  << ". Expected 2xN or Nx2, got ("
                                  << tmp.size() << " x "
                                  << (tmp.empty() ? 0 : tmp[0].size()) << ")"
                                  << std::endl;
                    }
    #else
                    std::cout << "Error: input file " << rch_file
                              << " is HDF5, but this executable was built without HDF5 support."
                              << std::endl;
    #endif
                }
                else
                {
                    // ASCII is always N x 2
                    std::vector<std::string> names;
                    std::vector<std::vector<double>> data_Nx2;

                    const bool tf = readLinearData(rch_file, names, data_Nx2, 5000000);

                    if (!tf)
                    {
                        std::cout << "Error: failed to read recharge file " << rch_file << std::endl;
                    }
                    else if (data_Nx2.empty())
                    {
                        std::cout << "Error: recharge file " << rch_file
                                  << " contains no data." << std::endl;
                    }
                    else if (data_Nx2[0].size() != 2)
                    {
                        std::cout << "Error: recharge ASCII file " << rch_file
                                  << " must have exactly 2 columns, but has "
                                  << data_Nx2[0].size() << std::endl;
                    }
                    else if (static_cast<int>(data_Nx2.size()) != Ncells)
                    {
                        std::cout << "Error: size mismatch between Depth and Recharge in "
                                  << rch_file << ". Depth size = " << Ncells
                                  << ", Recharge rows = " << data_Nx2.size() << std::endl;
                    }
                    else
                    {
                        Rch.resize(Ncells, 0.0);
                        SurfPerc.resize(Ncells, 0.0);

                        for (int i = 0; i < Ncells; ++i)
                        {
                            Rch[i]      = data_Nx2[i][0];
                            SurfPerc[i] = data_Nx2[i][1];
                        }

                        readSuccess = 1;
                    }
                }
            }
            catch (const std::exception& e)
            {
                std::cout << "Error reading recharge file " << rch_file
                          << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
        }

        world.barrier();
        sendScalarFromRoot2AllProc<int>(readSuccess, world);

        if (readSuccess == 0)
        {
            return false;
        }

        if (world.rank() > 0)
        {
            Rch.resize(Ncells, 0.0);
            SurfPerc.resize(Ncells, 0.0);
        }

        world.barrier();
        sendVectorFromRoot2AllProc<double>(Rch, world);
        sendVectorFromRoot2AllProc<double>(SurfPerc, world);

        if (PrintMatrices)
        {
            printVectorForAllProc<double>(Rch, world, 3483753, 3483753 + 20);
            printVectorForAllProc<double>(SurfPerc, world, 3483753, 3483753 + 20);
        }

        return true;
    }

}

#endif //MANTISSA_MS_UNSAT_H
