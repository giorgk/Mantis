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

        bool read(const std::string& filename, boost::mpi::communicator &world, int n_reserve);
        int getHRU(int IJ) const;

    private:
        bool read_HDF5(const std::string& filename, std::vector<int> &V, boost::mpi::communicator &world);
        std::vector<int> hrus;
        int Ncells;
    };

    bool HRU_Raster::read(const std::string& filename, boost::mpi::communicator &world, int n_reserve) {
        // Reset state first
        hrus.clear();
        Ncells = 0;

        std::string ext = getExtension(filename);
        std::vector<int> tmp;
        bool tf = false;
        if (ext == "h5" || ext == "H5"){
            tf = read_HDF5(filename, tmp, world);
        }
        else{
            if (world.rank() == 0) {
                std::cout << "Reading " << filename << std::endl;
            }
            tf = RootReadsVectorFileDistrib<int>(filename, tmp, world, 5000000,n_reserve);
        }
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
            printVectorForAllProc(tmp, world, 345, 353);
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
        hrus = tmp;

        return true;

    }

    inline int HRU_Raster::getHRU(int IJ) const {
        if (IJ < 0 || IJ >= Ncells ){
            return -9;
        }
        return hrus[IJ];
    }

    bool HRU_Raster::read_HDF5(const std::string& filename, std::vector<int> &V,
        boost::mpi::communicator &world) {
#if _USEHF > 0
        int readSuccess = 1;
        if (world.rank() == 0) {
            std::cout << "Reading " << filename << std::endl;
            const std::string NameSet("HRUs");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);

            // Read as 2D to preserve shape
            dataset.read(V);

            if (V.empty()) {
                std::cout << "Error: empty dataset in " << filename << std::endl;
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


#endif
        return false;
    }

}

#endif //MANTISSA_MS_HRU_RASTER_H
