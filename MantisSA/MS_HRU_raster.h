//
// Created by giorg on 11/11/2024.
//

#ifndef MANTISSA_MS_HRU_RASTER_H
#define MANTISSA_MS_HRU_RASTER_H

#include <string>


namespace MS{
    class HRU_Raster{
    public:
        HRU_Raster(){}
        bool read(std::string filename, boost::mpi::communicator &world);
        int getHRU(int IJ);
    private:
        std::vector<int> hrus;
        int Ncells;
    };

    bool HRU_Raster::read(std::string filename, boost::mpi::communicator &world) {
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
#if _USEHF>0
            const std::string NameSet("HRUs");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
            dataset.read(hrus);
            Ncells = hrus.size();
            return true;
#endif
        }
        else{
            std::vector<std::vector<int>> tmp;
            bool tf = RootReadsMatrixFileDistrib<int>(filename, tmp, 1, true, world, 5000000);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(tmp,world,0,1,345,353);}
            hrus = tmp[0];
            Ncells = hrus.size();
            return true;
        }

    }

    int HRU_Raster::getHRU(int IJ) {
        if (IJ < 0 || IJ >= Ncells ){
            return -9;
        }
        return hrus[IJ];
    }
}

#endif //MANTISSA_MS_HRU_RASTER_H
