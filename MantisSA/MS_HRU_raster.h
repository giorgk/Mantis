//
// Created by giorg on 11/11/2024.
//

#ifndef MANTISSA_MS_HRU_RASTER_H
#define MANTISSA_MS_HRU_RASTER_H

#include <string>


namespace MS{
    class HRU_Raster{
    public:
        HRU_Raster()
        : Ncells(0)
        {}

        bool read(const std::string& filename, boost::mpi::communicator &world);
        int getHRU(int IJ) const;

    private:
        std::vector<int> hrus;
        int Ncells;
    };

    bool HRU_Raster::read(const std::string& filename, boost::mpi::communicator &world) {
        // Reset state first
        hrus.clear();
        Ncells = 0;

        std::string ext = getExtension(filename);
        if (ext == "h5" || ext == "H5"){
#if _USEHF>0
            try {
                if (world.rank() == 0)
                {
                    std::cout << "Reading HDF5 HRU raster file " << filename << std::endl;
                }

                const std::string NameSet("HRUs");
                HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);

                dataset.read(hrus);
                Ncells = static_cast<int>(hrus.size());

                if (Ncells <= 0)
                {
                    if (world.rank() == 0)
                    {
                        std::cout << "Error: dataset 'HRUs' in file " << filename
                                  << " is empty." << std::endl;
                    }
                    return false;
                }

                return true;
            }
            catch (const std::exception& e) {
                if (world.rank() == 0)
                {
                    std::cout << "Error reading HDF5 HRU raster file " << filename
                              << ": " << e.what() << std::endl;
                }
                return false;
            }
#else
            if (world.rank() == 0)
            {
                std::cout << "Error: input file " << filename
                          << " is HDF5, but this executable was built without HDF5 support."
                          << std::endl;
            }
            return false;
#endif
        }
        else{
            std::vector<std::vector<int>> tmp;
            bool tf = RootReadsMatrixFileDistribFlat<int>(filename, tmp, 1, world, 5000000);

            if (!tf)
            {
                if (world.rank() == 0)
                {
                    std::cout << "Error: failed to read HRU raster file " << filename << std::endl;
                }
                return false;
            }

            if (tmp.empty())
            {
                if (world.rank() == 0)
                {
                    std::cout << "Error: HRU raster file " << filename
                              << " contains no data." << std::endl;
                }
                return false;
            }

            if (PrintMatrices)
            {
                printMatrixForAllProc(tmp, world, 345, 353, 0, 1);
            }

            Ncells = static_cast<int>(tmp.size());

            if (Ncells <= 0)
            {
                if (world.rank() == 0)
                {
                    std::cout << "Error: HRU raster file " << filename
                              << " contains an empty HRU vector." << std::endl;
                }
                return false;
            }

            hrus.resize(Ncells);

            for (int i = 0; i < Ncells; ++i)
            {
                if (tmp[i].size() != 1)
                {
                    if (world.rank() == 0)
                    {
                        std::cout << "Error: HRU raster row " << i
                                  << " does not have exactly 1 column." << std::endl;
                    }
                    return false;
                }
                hrus[i] = tmp[i][0];
            }

            return true;
        }

    }

    inline int HRU_Raster::getHRU(int IJ) const {
        if (IJ < 0 || IJ >= Ncells ){
            return -9;
        }
        return hrus[IJ];
    }
}

#endif //MANTISSA_MS_HRU_RASTER_H
