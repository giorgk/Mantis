//
// Created by giorg on 2/24/2023.
//

#ifndef MANTISSERVER_BRASTER_H
#define MANTISSERVER_BRASTER_H

namespace mantisServer{

    class BackgroundRaster{
    public:
        BackgroundRaster(){}
        bool readData(std::string filename, int Nr, int Nc, int Ncells);
        int IJ(int row, int col) const;
        int linear_index(int row, int col);
        bool locateCell(double x, double y, int& row, int& col, int& activeCellIndex) const;
        void getSurroundingPixels(int row, int col, int Ndepth, std::vector<int>& lin_inds);
        int Nr(){return Nrows;}
        int Nc(){return Ncols;}
        int Ncell(){return Ncells;}
        void setGridLocation(double x, double y, double cs);
        void getGridLocation(double &x, double &y, double &cs);
        void cellCoords(int r, int c, double &x, double &y);
    private:
        int Nrows;
        int Ncols;
        int Ncells;
        std::vector<std::vector<int>> raster;
        bool bHFread;
        double Xorig;
        double Yorig;
        double cellSize;
    };

    void BackgroundRaster::setGridLocation(double x, double y, double cs) {
        Xorig = x;
        Yorig = y;
        cellSize = cs;
    }

    void BackgroundRaster::getGridLocation(double &x, double &y, double &cs) {
        x = Xorig;
        y = Yorig;
        cs = cellSize;
    }

    void BackgroundRaster::cellCoords(int r, int c, double &x, double &y) {
        x = Xorig + cellSize/2 + cellSize*(c);
        // For the Y the row numbers start from the top
        y =  Yorig + cellSize*Nrows - cellSize/2 - cellSize*(r);
    }

    bool BackgroundRaster::readData(std::string filename, int Nr, int Nc, int Ncell) {
        auto start = std::chrono::high_resolution_clock::now();
        Nrows = Nr;
        Ncols = Nc;
        Ncells = Ncell;
        bHFread = false;
        raster.clear();
        if (Nrows <= 0 || Ncols <= 0 || Ncells <= 0) {
            std::cout << "Raster dimensions and active-cell count must be positive" << std::endl;
            return false;
        }
#if _USEHF>0
        if (isHdf5Filename(filename)){
            try {
                std::cout << "Reading " << filename << std::endl;
                const std::string NameSet("Raster");
                HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
                std::vector<std::vector<int>> tmp_raster;
                dataset.read(tmp_raster);

                bool isFullRaster = tmp_raster.size() == static_cast<std::size_t>(Ncols);
                if (isFullRaster) {
                    for (std::size_t col = 0; col < tmp_raster.size(); ++col) {
                        if (tmp_raster[col].size() != static_cast<std::size_t>(Nrows)) {
                            isFullRaster = false;
                            break;
                        }
                    }
                }

                bool isCoordinateList = tmp_raster.size() == static_cast<std::size_t>(Ncells);
                if (isCoordinateList) {
                    for (std::size_t i = 0; i < tmp_raster.size(); ++i) {
                        if (tmp_raster[i].size() != 2) {
                            isCoordinateList = false;
                            break;
                        }
                    }
                }

                if (isFullRaster) {
                    std::vector<bool> seen(static_cast<std::size_t>(Ncells), false);
                    int activeCells = 0;
                    for (int col = 0; col < Ncols; ++col) {
                        for (int row = 0; row < Nrows; ++row) {
                            const int value = tmp_raster[col][row];
                            if (value == -1) continue;
                            if (value < 0 || value >= Ncells || seen[value]) {
                                std::cout << "Invalid or duplicate active-cell index " << value
                                          << " in " << filename << std::endl;
                                return false;
                            }
                            seen[value] = true;
                            ++activeCells;
                        }
                    }
                    if (activeCells != Ncells) {
                        std::cout << "Raster in " << filename << " contains " << activeCells
                                  << " active indices but " << Ncells << " were expected" << std::endl;
                        return false;
                    }
                    raster = tmp_raster;
                    bHFread = true;
                }
                else if (isCoordinateList) {
                    raster.assign(Nrows, std::vector<int>(Ncols,-1));
                    for (int i = 0; i < Ncells; ++i){
                        const int r = tmp_raster[i][0];
                        const int c = tmp_raster[i][1];
                        if (r < 0 || c < 0 || r >= Nrows || c >= Ncols || raster[r][c] != -1) {
                            std::cout << "Invalid or duplicate raster pixel (" << r << "," << c
                                      << ") in " << filename << std::endl;
                            return false;
                        }
                        raster[r][c] = i;
                    }
                }
                else{
                    std::cout << "The raster dimensions in " << filename
                              << " do not agree with the Raster options" << std::endl;
                    return false;
                }

                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
                return true;
            }
            catch (const std::exception& e) {
                std::cout << "HDF5 error while reading background raster " << filename
                          << ": " << e.what() << std::endl;
                return false;
            }
            catch (...) {
                std::cout << "Unknown HDF5 error while reading background raster "
                          << filename << std::endl;
                return false;
            }
        }
#else
        if (isHdf5Filename(filename)) {
            std::cout << "This MantisServer build does not include HDF5 support" << std::endl;
            return false;
        }
#endif

        std::ifstream ifile;
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
                if (!getline(ifile, line)) {
                    std::cout << "Not enough coordinate rows in " << filename << std::endl;
                    return false;
                }
                std::istringstream inp(line.c_str());
                if (!(inp >> r >> c)) {
                    std::cout << "Invalid coordinate row in " << filename << std::endl;
                    return false;
                }
                if (r >= 0 && c >= 0 && r < Nrows && c < Ncols && raster[r][c] == -1)
                    raster[r][c] = i;
                else{
                    std::cout << "Invalid or duplicate raster pixel (" << r << "," << c
                              << ") in " << filename << std::endl;
                    return false;
                }
            }
            ifile.close();
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
        return true;
    }

    int BackgroundRaster::IJ(int row, int col) const {
        if (row >=0 && col >=0 && row < Nrows && col < Ncols) {
            if (bHFread)
                return raster[col][row];
            else
                return raster[row][col];
        }
        else
            return -1;
    }

    bool BackgroundRaster::locateCell(double x, double y, int& row, int& col,
                                      int& activeCellIndex) const {
        row = -1;
        col = -1;
        activeCellIndex = -1;
        if (cellSize <= 0.0 || x < Xorig || y < Yorig) {
            return false;
        }

        col = static_cast<int>(std::floor((x - Xorig) / cellSize));
        row = Nrows - 1 -
              static_cast<int>(std::floor((y - Yorig) / cellSize));
        if (row < 0 || col < 0 || row >= Nrows || col >= Ncols) {
            row = -1;
            col = -1;
            return false;
        }

        activeCellIndex = IJ(row, col);
        return activeCellIndex >= 0;
    }

    int BackgroundRaster::linear_index(int row, int col) {
        return Nrows * col + row;
    }

    void BackgroundRaster::getSurroundingPixels(int row, int col, int Ndepth, std::vector<int> &lin_inds) {
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
