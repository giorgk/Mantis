//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_RASTER_H
#define MANTISSA_MS_RASTER_H

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

namespace MS{
    class BackgroundRaster{
    public:
        BackgroundRaster(){}
        bool readData(std::string filename, int Nr, int Nc, int Ncells);
        int IJ(int row, int col);
        int linear_index(int row, int col);
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

    bool BackgroundRaster::readData(std::string filename, int Nr, int Nc, int Ncell){
        Nrows = Nr;
        Ncols = Nc;
        Ncells = Ncell;
#if _USEHF>0

        const std::string NameSet("Raster");
        HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
        HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
        dataset.read(raster);

        return true;

#endif
    }

    int BackgroundRaster::IJ(int row, int col) {
        if (row >=0 && col >=0 && row < Nrows && col < Ncols) {
                return raster[col][row];
        }
        else
            return -1;
    }

    int BackgroundRaster::linear_index(int row, int col) {
        return Nrows * col + row;
    }
}

#endif //MANTISSA_MS_RASTER_H