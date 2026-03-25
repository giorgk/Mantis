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

#include <unordered_map>
#include <cstdint>
#include <utility>

#include "MS_mpi_utils.h"

namespace MS{
    class BackgroundRaster{
    public:
        BackgroundRaster();
        bool readData(const std::string& filename, int Nr, int Nc, int Ncells,
            boost::mpi::communicator &world);
        int IJ(int row, int col) const;
        bool cellFromIndex(int idx, int& row, int& col) const;
        int Nr() const{return Nrows;}
        int Nc() const{return Ncols;}
        int Ncell() const{return Ncells;}
        void setGridLocation(double x, double y, double cs);
        void getGridLocation(double &x, double &y, double &cs) const;
        void cellCoords(int r, int c, double &x, double &y) const;

    private:
        static std::uint64_t packRC(int row, int col) {
            return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(row)) << 32) |
                   static_cast<std::uint32_t>(col);
        }
        int Nrows;
        int Ncols;
        int Ncells;
        bool bHFread;
        double Xorig;
        double Yorig;
        double cellSize;
        //std::vector<std::vector<int>> raster;
        std::unordered_map<std::uint64_t, int> rc_to_idx;
        std::vector<std::pair<int,int>> idx_to_rc;

        bool readData_h5(const std::string& filename, std::vector<std::vector<int>>& rasterCol,
                 boost::mpi::communicator& world);
    };

    inline BackgroundRaster::BackgroundRaster()
        : Nrows(0), Ncols(0), Ncells(0),
          bHFread(false), Xorig(0.0), Yorig(0.0), cellSize(0.0)
    {}

    inline void BackgroundRaster::setGridLocation(double x, double y, double cs) {
        Xorig = x;
        Yorig = y;
        cellSize = cs;
    }

    inline void BackgroundRaster::getGridLocation(double &x, double &y, double &cs) const {
        x = Xorig;
        y = Yorig;
        cs = cellSize;
    }

    inline void BackgroundRaster::cellCoords(int r, int c, double &x, double &y) const {
        x = Xorig + cellSize/2 + cellSize*(c);
        // For the Y the row numbers start from the top
        y =  Yorig + cellSize*Nrows - cellSize/2 - cellSize*(r);
    }

    inline bool BackgroundRaster::readData(const std::string& filename, int Nr, int Nc, int Ncell,
                                           boost::mpi::communicator &world){
        Nrows = Nr;
        Ncols = Nc;
        Ncells = Ncell;
        std::vector<std::vector<int>> rasterCol;
        if (world.rank() == 0) {
            std::cout << "Reading raster file " << filename << std::endl;
        }
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
            const bool tf = readData_h5(filename, rasterCol, world);
            if (!tf)
            {
                return false;
            }
        }
        else{
            const bool tf = RootReadsMatrixFileDistribFlat(filename, rasterCol,2, world, 5000000);
            world.barrier();
            if (!tf){return false;}
        }

        if (PrintMatrices){
            printMatrixForAllProc<int>(rasterCol, world, 0, 10, 0, 2);
        }

        if (rasterCol.empty()) {
            std::cout << "Error: raster file " << filename
                      << " contains no valid rows" << std::endl;
            return false;
        }

        world.barrier();
        std::cout << "Build raster indexing ..." << std::endl;

        rc_to_idx.clear();
        rc_to_idx.reserve(rasterCol.size());
        idx_to_rc.clear();
        idx_to_rc.reserve(rasterCol.size());

        for (std::size_t i = 0; i < rasterCol.size(); ++i) {
            if (rasterCol[i].size() != 2) {
                std::cout << "Error: raster row " << i << " does not have 2 columns" << std::endl;
                return false;
            }

            const int r = rasterCol[i][0];
            const int c = rasterCol[i][1];

            if (r < 0 || r >= Nrows || c < 0 || c >= Ncols) {
                std::cout << "I can't assign pixel (" << r << "," << c
                          << ") in raster map" << std::endl;
                return false;
            }

            const int idx = static_cast<int>(i);
            const std::uint64_t key = packRC(r, c);
            auto ret = rc_to_idx.emplace(key, idx);
            if (!ret.second) {
                std::cout << "Error: duplicate raster cell (" << r << "," << c << ")"
                          << std::endl;
                return false;
            }

            idx_to_rc.emplace_back(r, c);
        }

        if (Ncells > 0 && static_cast<int>(idx_to_rc.size()) != Ncells) {
            std::cout << "Warning: Ncells = " << Ncells
                      << " but raster file contains " << idx_to_rc.size()
                      << " active cells" << std::endl;
        }

        return true;
    }

    inline int BackgroundRaster::IJ(int row, int col) const {
        if (row < 0 || row >= Nrows || col < 0 || col >= Ncols) {
            return -1;
        }
        const auto it = rc_to_idx.find(packRC(row, col));
        if (it == rc_to_idx.end()) {
            return -1;
        }
        return it->second;
    }

    inline bool BackgroundRaster::cellFromIndex(int idx, int& row, int& col) const
    {
        if (idx < 0 || idx >= static_cast<int>(idx_to_rc.size())) {
            return false;
        }

        row = idx_to_rc[idx].first;
        col = idx_to_rc[idx].second;
        return true;
    }

    inline bool BackgroundRaster::readData_h5(const std::string& filename,
                                          std::vector<std::vector<int>>& rasterCol,
                                          boost::mpi::communicator& world) {
        int readSuccess = 0;
        rasterCol.clear();

        if (world.rank() == 0) {
#if _USEHF > 0
            try {
                const std::string NameSet("Raster");
                HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);

                std::vector<std::vector<int>> tmp;
                dataset.read(tmp);

                if (tmp.empty())
                {
                    std::cout << "Error: dataset 'Raster' is empty in " << filename << std::endl;
                }
                else if (tmp.size() == 2)
                {
                    // 2 x N  -> transpose to N x 2
                    const std::size_t N = tmp[0].size();
                    rasterCol.assign(N, std::vector<int>(2, 0));

                    for (std::size_t i = 0; i < N; ++i)
                    {
                        rasterCol[i][0] = tmp[0][i];
                        rasterCol[i][1] = tmp[1][i];
                    }

                    readSuccess = 1;
                }
                else if (tmp[0].size() == 2)
                {
                    // already N x 2
                    rasterCol = tmp;
                    readSuccess = 1;
                }
                else
                {
                    std::cout << "Error: cannot interpret Raster dimensions in " << filename
                              << ". Expected Nx2 or 2xN, got ("
                              << tmp.size() << " x "
                              << (tmp.empty() ? 0 : tmp[0].size()) << ")" << std::endl;
                }
            }
            catch (const std::exception& e) {
                std::cout << "Error reading HDF5 raster file " << filename
                      << ": " << e.what() << std::endl;
            }
#else
            std::cout << "Error: HDF5 support is disabled but file " << filename
                  << " has .h5 extension" << std::endl;
#endif
        }
        world.barrier();
        sendScalarFromRoot2AllProc<int>(readSuccess, world);

        if (readSuccess == 0)
        {
            return false;
        }

        sendFlatMatrixFromRoot2AllProc<int>(rasterCol, world);

        return true;
    }

}

#endif //MANTISSA_MS_RASTER_H
