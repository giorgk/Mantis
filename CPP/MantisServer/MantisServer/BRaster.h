//
// Created by giorg on 2/24/2023.
//

#ifndef MANTISSERVER_BRASTER_H
#define MANTISSERVER_BRASTER_H

namespace mantisServer{

    class BackroundRaster{
    public:
        BackroundRaster(){}
        bool readData(std::string filename, int Nr, int Nc, int Ncells);
        int IJ(int row, int col);
        void getSurroundingPixels(int row, int col, int Ndepth, std::vector<int>& lin_inds);
        int Nr(){return Nrows;}
        int Nc(){return Ncols;}
        int Ncell(){return Ncells;}
    private:
        int Nrows;
        int Ncols;
        int Ncells;
        std::vector<std::vector<int>> raster;
        bool bHFread;
    };

    bool BackroundRaster::readData(std::string filename, int Nr, int Nc, int Ncell) {
        auto start = std::chrono::high_resolution_clock::now();
        Nrows = Nr;
        Ncols = Nc;
        Ncells = Ncell;
#if _USEHF>0
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
            const std::string NameSet("Raster");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
            dataset.read(raster);
            if (raster.size() != Ncols){
                std::cout << "The number of columns of the raster (" << raster.size() << ") is not equal to Ncols: " << Ncols << std::endl;
            }
            if (raster[0].size() != Nrows){
                std::cout << "The number of rows of the raster (" << raster[0].size() << ") is not equal to Nrows: " << Nrows << std::endl;
            }
            bHFread = true;
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
            return true;
        }
#endif

        bHFread = false;
        std::ifstream ifile;
        raster.clear();
        raster.resize(Nrows,std::vector<int>(Ncols,-1));
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            int r, c;
            for (int i = 0; i < Ncells; ++i){
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> r;
                inp >> c;
                if (r < Nrows && c < Ncols)
                    raster[r][c] = i;
                else{
                    std::cout << "I can't assign pixel (" << r << "," << c << ") in raster map" << std::endl;
                }
            }
            ifile.close();
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
        return true;
    }

    int BackroundRaster::IJ(int row, int col) {
        if (row >=0 && col >=0 && row < Nrows && col < Ncols) {
            if (bHFread)
                return raster[col][row];
            else
                return raster[row][col];
        }
        else
            return -1;
    }

    void BackroundRaster::getSurroundingPixels(int row, int col, int Ndepth, std::vector<int> &lin_inds) {
        if (Ndepth == 0){
            int lid = IJ(row,col);
            if ( lid >= 0)
                lin_inds.push_back(lid);
        }
        else{
            for (int i = row - Ndepth; i <= row + Ndepth; ++i){
                int lid = IJ(row,col);
                if (lid == -1)
                    continue;
                for (int j = col - Ndepth; j < col + Ndepth; ++j){
                    lid = IJ(row,col);
                    if (lid == -1)
                        continue;
                    lin_inds.push_back(lid);
                }
            }
        }
    }
}

#endif //MANTISSERVER_BRASTER_H
