//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_NPSAT_DATA_H
#define MANTISSA_NPSAT_DATA_H

#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <zlib.h>

#include "MS_structures.h"
#include "MS_urf_calc.h"
#include "MS_unsat.h"
#include "MS_HRU_raster.h"
#include "MS_mpi_utils.h"
#include "MSdebug.h"

namespace MS{

    bool has_gz_extension(const std::string& filename) {
        return filename.size() >= 3 &&
               filename.compare(filename.size() - 3, 3, ".gz") == 0;
    }

    inline bool readNPSATdata_h5_from_root(const std::string& filename,
                                       Matrix<int>& ints,
                                       Matrix<double>& dbls,
                                       Matrix<double>& msas,
                                       std::vector<int>& version,
                                       std::vector<double>& porosities,
                                       boost::mpi::communicator& world) {
        int readSuccess = 0;

        ints.clear();
        dbls.clear();
        msas.clear();
        version.clear();
        porosities.clear();

        if (world.rank() == 0) {
#if _USEHF > 0
            try {
                std::cout << "Reading " << filename << std::endl;
                const std::string INT_NameSet("INT");
                const std::string DBL_NameSet("DBL");
                const std::string MSA_NameSet("MSA");
                const std::string VER_NameSet("Version");
                const std::string POR_NameSet("Por");

                HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

                // --- INT
                {
                    std::vector<std::vector<int>> tmp;
                    HDFNfile.getDataSet(INT_NameSet).read(tmp);

                    const size_t nRows = tmp.size();
                    const size_t nCols = (nRows > 0 ? tmp[0].size() : 0);

                    ints.allocate(nRows, nCols);
                    for (size_t i = 0; i < nRows; ++i)
                        for (size_t j = 0; j < nCols; ++j)
                            ints(i, j) = tmp[i][j];
                }

                // --- DBL
                {
                    std::vector<std::vector<double>> tmp;
                    HDFNfile.getDataSet(DBL_NameSet).read(tmp);

                    const size_t nRows = tmp.size();
                    const size_t nCols = (nRows > 0 ? tmp[0].size() : 0);

                    dbls.allocate(nRows, nCols);
                    for (size_t i = 0; i < nRows; ++i)
                        for (size_t j = 0; j < nCols; ++j)
                            dbls(i, j) = tmp[i][j];
                }

                // --- MSA
                {
                    std::vector<std::vector<double>> tmp;
                    HDFNfile.getDataSet(MSA_NameSet).read(tmp);

                    const size_t nRows = tmp.size();
                    const size_t nCols = (nRows > 0 ? tmp[0].size() : 0);

                    msas.allocate(nRows, nCols);
                    for (size_t i = 0; i < nRows; ++i)
                        for (size_t j = 0; j < nCols; ++j)
                            msas(i, j) = tmp[i][j];
                }

                // --- 1D datasets (unchanged)
                HDFNfile.getDataSet(VER_NameSet).read(version);
                HDFNfile.getDataSet(POR_NameSet).read(porosities);

                readSuccess = 1;
            }
            catch (const std::exception& e)
            {
                std::cout << "Error reading HDF5 NPSAT file " << filename
                          << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
#else
            std::cout << "Error: HDF5 support is disabled but file "
                      << filename << " has .h5 extension" << std::endl;
            readSuccess = 0;
#endif
        }

        world.barrier();
        sendScalarFromRoot2AllProc<int>(readSuccess, world);

        if (readSuccess == 0)
        {
            return false;
        }

        sendMatrixFromRoot2AllProc<int>(ints, world);
        sendMatrixFromRoot2AllProc<double>(dbls, world);
        sendMatrixFromRoot2AllProc<double>(msas, world);
        sendVectorFromRoot2AllProc<int>(version, world);
        sendVectorFromRoot2AllProc<double>(porosities, world);

        return true;
    }


    inline bool readNPSATdata(const std::string& filename, WELLS &wells, int por, int Nyears,
                              BackgroundRaster& braster, RiverOptions &rivopt,
                              boost::mpi::communicator &world){

        const std::string ext = getExtension(filename);

        Matrix<int> ints;
        Matrix<double> dbls;
        Matrix<double> msas;
        std::vector<int> version;
        std::vector<double> porosities;

        bool tf = false;

        if (ext == "h5" || ext == "H5") {
            tf = readNPSATdata_h5_from_root(filename, ints, dbls, msas, version, porosities, world);
            if (!tf)
            {
                return false;
            }
        }
        else {
            std::vector<int> info;
            bool tf = RootReadsNPSATinfo(filename + "info.dat", info, world);
            if (!tf){return false;}

            // Expected flattened format:
            // [version, nPor, por1, por2, ..., porN, intRows, intCols, dblRows, dblCols, msaRows, msaCols]
            if (info.size() < 2){
                std::cout << filename + "info.dat should contain version, porosity and the int, dbl, msa matrices size" << std::endl;
                return false;
            }

            const int nPor = info[1];
            if (nPor <= 0) {
                std::cout << filename + "info.dat has invalid number of porosities" << std::endl;
                return false;
            }

            const int idx = nPor + 2;
            if (static_cast<int>(info.size()) < idx + 6) {
                std::cout << filename + "info.dat should contain version, porosity values, "
                          << "and INT/DBL/MSA matrix sizes" << std::endl;
                return false;
            }

            version.push_back(info[0]);

            for (int i = 0; i < nPor; ++i) {
                porosities.push_back(static_cast<double>(info[i + 2]));
            }

            const int intRows = info[idx + 0];
            const int intCols = info[idx + 1];
            const int dblRows = info[idx + 2];
            const int dblCols = info[idx + 3];
            const int msaRows = info[idx + 4];
            const int msaCols = info[idx + 5];

            if (world.rank() == 0) {std::cout << "Reading " << filename + "INT.dat" << std::endl;}
            tf = RootReadsMatrixFileDistrib<int>(filename + "INT.dat", ints, intCols, world, 500000, intRows);
            if (!tf){return false;}
            if (world.rank() == 0) {std::cout << "Reading " << filename + "DBL.dat" << std::endl;}
            tf = RootReadsMatrixFileDistrib<double>(filename + "DBL.dat", dbls, dblCols, world, 500000, dblRows);
            if (!tf){return false;}
            if (world.rank() == 0) {std::cout << "Reading " << filename + "MSA.dat" << std::endl;}
            tf = RootReadsMatrixFileDistrib<double>(filename + "MSA.dat", msas, msaCols, world, 500000, msaRows);
            if (!tf){return false;}

            if (static_cast<int>(ints.size()) != intRows) {
                std::cout << "Error: INT.dat row count does not match info.dat" << std::endl;
                return false;
            }
            if (static_cast<int>(dbls.size()) != dblRows) {
                std::cout << "Error: DBL.dat row count does not match info.dat" << std::endl;
                return false;
            }
            if (static_cast<int>(msas.size()) != msaRows) {
                std::cout << "Error: MSA.dat row count does not match info.dat" << std::endl;
                return false;
            }
        }

        if (PrintMatrices){
            printMatrixForAllProc<int>(ints, world, 0, 10, 0, ints.num_cols());
            printMatrixForAllProc<double>(dbls, world, 0, 10, 0, dbls.num_cols());
            printMatrixForAllProc<double>(msas, world, 0, 10, 0, msas.num_cols());
        }

        // ----------------------------
        // Global consistency checks
        // ----------------------------
        if (version.empty()) {
            std::cout << "Error: missing version information in " << filename << std::endl;
            return false;
        }

        if (porosities.empty()) {
            std::cout << "Error: missing porosity list in " << filename << std::endl;
            return false;
        }

        if (ints.empty() || dbls.empty() || msas.empty()) {
            std::cout << "Error: one or more NPSAT matrices are empty in " << filename << std::endl;
            return false;
        }

        if (ints.num_rows() != dbls.num_rows() || ints.num_rows() != msas.num_rows()) {
            std::cout << "Error: INT, DBL, and MSA row counts do not match in "
                      << filename << std::endl;
            return false;
        }

        const int ver = version[0];
        if (ver != 0 && ver != 1) {
            std::cout << "Error: unsupported NPSAT version " << ver << std::endl;
            return false;
        }

        int por_pos = -1;
        for (int i = 0; i < static_cast<int>(porosities.size()); ++i) {
            if (static_cast<int>(porosities[i]) == por) {
                por_pos = i;
                break;
            }
        }

        if (por_pos < 0) {
            std::cout << "Error: requested porosity " << por
                      << " was not found in file " << filename << std::endl;
            return false;
        }

        const int por_idx = 3 * por_pos;

        if (world.rank() == 0){
            std::cout << "Assembling NPSAT urf data ..." << std::endl;
        }

        wells.version = version[0];

        tmpWELLS tw;
        tmpWELLS::iterator itwtmp = tw.end();
        std::pair<tmpWELLS::iterator,bool> tmpret;

        //std::map<int, WELL>::iterator itw = wells.wells.end();
        //std::pair<std::map<int, WELL>::iterator,bool> ret;

        for (std::size_t i = 0; i < ints.size(); ++i){
            const int wellId = ints(i,0);

            if (itwtmp == tw.end() || itwtmp->first != wellId) {
                itwtmp = tw.find(wellId);
                if (itwtmp == tw.end()) {
                    tmpret = tw.insert(std::pair<int, NPSATTMP>(wellId, NPSATTMP()));
                    if (!tmpret.second) {
                        std::cout << "Error: failed to insert temporary well " << wellId << std::endl;
                        return false;
                    }
                    itwtmp = tmpret.first;
                }
            }

            itwtmp->second.Sid.push_back(ints(i,1));
            itwtmp->second.urfI.push_back(ints(i,2));
            itwtmp->second.urfJ.push_back(ints(i,3));
            if (version[0] == 0){
                itwtmp->second.inRiv.push_back(ints(i,4) == 1);
            }

            itwtmp->second.W.push_back(dbls(i,0));
            itwtmp->second.Len.push_back(dbls(i,1));

            if (version[0] == 1){
                itwtmp->second.rivDist.push_back(dbls(i,2));
            }

            itwtmp->second.m.push_back(msas(i,por_idx + 0));
            itwtmp->second.s.push_back(msas(i,por_idx + 1));
            itwtmp->second.a.push_back(msas(i,por_idx + 2));
            itwtmp->second.sumW += dbls(i,0);
        }

        //Normalize streamlines and split them according to the processors
        int count = 0;
        for (itwtmp = tw.begin(); itwtmp != tw.end(); ++itwtmp){
            if (count % world.size() != world.rank()){
                ++count;
                continue;
            }

            WELL w;

            for (std::size_t i = 0; i < itwtmp->second.urfJ.size(); ++i){
                //if (itwtmp->first == 21539){
                //    bool stop = true;
                //}

                STRML s;
                s.Sid = itwtmp->second.Sid[i];
                s.urfI = itwtmp->second.urfI[i]-1;
                s.urfJ = itwtmp->second.urfJ[i]-1;
                // Its ok if the source area is negative. that means loading is zero
                s.IJ = braster.IJ(s.urfI, s.urfJ);

                if (ver == 0){
                    s.inRiv = itwtmp->second.inRiv[i];
                }
                if (ver == 1){
                    s.rivInfl = calcRiverInfluence(itwtmp->second.rivDist[i], rivopt);
                    s.inRiv = false;
                }

                if (itwtmp->second.sumW == 0.0) {
                    std::cout << "Error: zero total weight for well " << itwtmp->first << std::endl;
                    return false;
                }

                //s.hru_idx = itwtmp->second.hru_idx[i];
                s.W = itwtmp->second.W[i]/itwtmp->second.sumW;
                s.Len = itwtmp->second.Len[i];
                s.m = itwtmp->second.m[i];
                s.s = itwtmp->second.s[i];
                s.a = itwtmp->second.a[i];
                calcURFs(Nyears,s.m, s.s, s.a, s.Len, s.urf);
                w.strml.push_back(s);
            }
            w.m_rmv.resize(Nyears,0.0);

            wells.wells.insert(std::pair<int,WELL>(itwtmp->first, w));
            count = count + 1;
            //std::cout << count << std::endl;
        }
        return true;
    }

    inline bool readInitSaltConc(const std::string &filename, WELLS &wells, boost::mpi::communicator &world){
        Matrix<double> eid_conc;
        if (world.rank() == 0) {
            std::cout << "Reading initial salt concentration from " << filename << std::endl;
        }
        const bool tf = RootReadsMatrixFileDistrib<double>(filename, eid_conc, 2, world,1000000, 250000);
        if (tf) {
            std::map<int, WELL>::iterator itw;
            for (std::size_t i = 0; i < eid_conc.size(); ++i){
                int eid = static_cast<int>(eid_conc(i,0));
                itw = wells.wells.find(eid);
                if (itw != wells.wells.end()){
                    itw->second.initConc = eid_conc(i,1);
                }
            }
        }
        return tf;
    }

    void linearizeWellBTCs(WELLS &W, std::vector<double> &v){
        v.clear();
        v.push_back(static_cast<double>(W.wells.size()));

        std::map<int, WELL>::iterator itw;
        for (itw = W.wells.begin(); itw != W.wells.end(); ++itw){
            v.push_back(static_cast<double>(itw->first));

            const std::vector<double> &btc = itw->second.wellBtc;
            for (unsigned int i = 0; i < btc.size(); ++i){
                v.push_back(btc[i]);
            }
        }
    }

    void linearizeBTCs(WELLS &W, std::vector<double> &v, std::vector<int> ids, SWAT_data &swat, HRU_Raster &hru_raster, int Nyears){
        v.clear();
        std::map<int, WELL>::iterator itw;
        int countWells = 0;
        for (unsigned int j = 0; j < ids.size(); ++j){
            itw = W.wells.find(ids[j]);
            if (itw != W.wells.end()){
                countWells = countWells + 1;
                v.push_back(static_cast<double>(itw->first));// Eid
                int countS = 0;
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru_idx = swat.hru_index(hru_raster.getHRU(itw->second.strml[i].IJ));
                    if (hru_idx >= 0){
                        countS = countS + 1;
                    }
                }
                v.push_back(static_cast<double>(countS));// Number of streamlines
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru = hru_raster.getHRU(itw->second.strml[i].IJ);
                    int hru_idx = swat.hru_index(hru);
                    if (hru_idx >= 0){
                        v.push_back(static_cast<double>(itw->second.strml[i].Sid)); // Sid
                        v.push_back(static_cast<double>(itw->second.strml[i].W)); // W
                        for (int k = 0; k < Nyears; ++k){
                            v.push_back(itw->second.strml[i].btc[k]);
                        }
                    }
                }
            }
        }
        v.insert(v.begin(), countWells);
    }

    void linearizeURFS(WELLS &W, std::vector<double> &v, std::vector<int> ids, SWAT_data &swat, HRU_Raster &hru_raster, int Nyears){
        v.clear();
        std::map<int, WELL>::iterator itw;
        int countWells = 0;
        for (unsigned int j = 0; j < ids.size(); ++j){
            itw = W.wells.find(ids[j]);
            if (itw != W.wells.end()){
                countWells = countWells + 1;
                v.push_back(static_cast<double>(itw->first));// Eid
                int countS = 0;
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru_idx = swat.hru_index(hru_raster.getHRU(itw->second.strml[i].IJ));
                    if (hru_idx >= 0){
                        countS = countS + 1;
                    }
                }
                v.push_back(static_cast<double>(countS));// Number of streamlines
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru = hru_raster.getHRU(itw->second.strml[i].IJ);
                    int hru_idx = swat.hru_index(hru);
                    if (hru_idx >= 0){
                        v.push_back(static_cast<double>(itw->second.strml[i].Sid)); // Sid
                        v.push_back(static_cast<double>(itw->second.strml[i].m)); // m
                        v.push_back(static_cast<double>(itw->second.strml[i].s)); // s
                        v.push_back(static_cast<double>(itw->second.strml[i].a)); // a
                        v.push_back(static_cast<double>(itw->second.strml[i].Len)); // Len
                        v.push_back(static_cast<double>(itw->second.strml[i].W)); // W
                        for (int k = 0; k < Nyears; ++k){
                            v.push_back(itw->second.strml[i].urf[k]);
                        }
                    }
                }

            }
        }
        v.insert(v.begin(), countWells);
    }

    void BundleDetailData(WELLS &W, std::vector<double> &v, std::vector<double> &m,
                          std::vector<int> ids,
                          UNSAT &UN, HRU_Raster &hru_raster, SWAT_data &swat, int Nyears){
        v.clear();

        std::map<int, WELL>::iterator itw;
        int countWells = 0;
        for (unsigned int j = 0; j < ids.size(); ++j){
            itw = W.wells.find(ids[j]);
            //std::cout << itw->first << std::endl;
            if (itw != W.wells.end()){
                countWells = countWells + 1;
                v.push_back(static_cast<double>(itw->first));// Eid
                m.push_back(static_cast<double>(itw->first));// Eid m
                int countS = 0;
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru_idx = swat.hru_index(hru_raster.getHRU(itw->second.strml[i].IJ));
                    if (hru_idx >= 0){
                        countS = countS + 1;
                    }
                }
                v.push_back(static_cast<double>(countS));// Number of streamlines
                m.push_back(static_cast<double>(countS));// Number of streamlines m
                for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                    int hru = hru_raster.getHRU(itw->second.strml[i].IJ);
                    int hru_idx = swat.hru_index(hru);
                    if (hru_idx >= 0){
                        v.push_back(static_cast<double>(itw->second.strml[i].Sid)); // Sid
                        m.push_back(static_cast<double>(itw->second.strml[i].Sid)); // Sid m
                        m.push_back(static_cast<double>(itw->second.strml[i].wellSourceId)); // Well id that this well receives its applied water from
                        v.push_back(static_cast<double>(itw->second.strml[i].urfI)); // UrfI
                        v.push_back(static_cast<double>(itw->second.strml[i].urfJ)); // UrfJ
                        v.push_back(static_cast<double>(hru)); // hru
                        v.push_back(static_cast<double>(hru_idx)); // hru_idx
                        v.push_back(UN.getDepth(itw->second.strml[i].IJ)); // Depth
                        v.push_back(UN.getRch(itw->second.strml[i].IJ)); // Recharge
                        v.push_back(UN.getSurfPerc(itw->second.strml[i].IJ)); // Surface percentage
                        v.push_back(static_cast<double>(itw->second.strml[i].rivInfl)); // RiverInfluence
                        for (int k = 0; k < Nyears; ++k) {
                            v.push_back(itw->second.strml[i].lf_conc[k]);
                            m.push_back(itw->second.strml[i].gw_conc[k]);
                        }
                    }
                }
            }
        }
        v.insert(v.begin(), countWells);
        m.insert(m.begin(), countWells);
    }

    void printWELLSfromAllProc(std::vector<std::vector<double>> &AllProcBTC, std::string filename,
                               int Nyears, bool compress = false){

        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // Open output stream
        if (use_compression){
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else{
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        std::ostringstream line;
        line << "Eid";
        for (int i = 1; i <= Nyears; ++i)
            line << ", btc" << i;
        line << "\n";
        write_line(line.str());

        for (const auto& procData : AllProcBTC){
            int idx = 0;
            int Nbtc = static_cast<int>(procData[idx++]);
            for (int j = 0; j < Nbtc; ++j){
                int eid = static_cast<int>(procData[idx++]);
                std::ostringstream row;
                row << std::scientific << std::setprecision(6);
                row << eid;
                for (int iyr = 0; iyr < Nyears; ++iyr){
                    row << ", " << procData[idx++];
                }
                row << "\n";
                write_line(row.str());
            }
        }

        // Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    void printMfeedFromAllProc(std::vector<std::vector<double>> &AllProcData, std::string filename,
                               int Nyears, bool compress = false){
        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // Open output stream
        if (use_compression){
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else{
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        // Set common formatting
        std::ostringstream line;
        line << std::scientific << std::setprecision(6);
        line << "Eid, Sid, Aid";
        for (int i = 1; i <= Nyears; ++i)
            line << ", MF" << i;
        line << "\n";
        write_line(line.str());

        for (const auto& procData : AllProcData){
            int idx = 0;
            int Nw = static_cast<int>(procData[idx++]);
            for (int j = 0; j < Nw; ++j){
                int eid = static_cast<int>(procData[idx++]);
                int Ns  = static_cast<int>(procData[idx++]);
                for (int k = 0; k < Ns; ++k){
                    int sid     = static_cast<int>(procData[idx++]);
                    int aid    = static_cast<int>(procData[idx++]);

                    std::ostringstream row;
                    row << std::scientific << std::setprecision(6);
                    row << eid << ", " << sid << ", " << aid;

                    for (int iyr = 0; iyr < Nyears; ++iyr) {
                        row << ", " << procData[idx++];
                    }
                    row << "\n";
                    write_line(row.str());
                }
            }
        }

        // Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    void printBTCssFromAllProc(std::vector<std::vector<double>> &AllProcData, std::string filename,
                               int Nyears, bool compress = false){

        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // Open output stream
        if (use_compression){
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else{
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        // Set common formatting
        std::ostringstream line;
        line << std::scientific << std::setprecision(6);
        line << "Eid, Sid, W";
        for (int i = 1; i <= Nyears; ++i)
            line << ", btc" << i;
        line << "\n";
        write_line(line.str());

        for (const auto& procData : AllProcData){
            int idx = 0;
            int Nw = static_cast<int>(procData[idx++]);
            for (int j = 0; j < Nw; ++j){
                int eid = static_cast<int>(procData[idx++]);
                int Ns  = static_cast<int>(procData[idx++]);
                for (int k = 0; k < Ns; ++k){
                    int sid     = static_cast<int>(procData[idx++]);
                    double w        = procData[idx++];

                    std::ostringstream row;
                    row << std::scientific << std::setprecision(6);
                    row << eid << ", " << sid << ", " << w;

                    for (int iyr = 0; iyr < Nyears; ++iyr) {
                        row << ", " << procData[idx++];
                    }
                    row << "\n";
                    write_line(row.str());
                }
            }
        }
// Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    void printURFsFromAllProc(std::vector<std::vector<double>> &AllProcData, std::string filename,
                              int Nyears, bool compress = false){

        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // Open output stream
        if (use_compression){
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else{
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        // Set common formatting
        std::ostringstream line;
        line << std::scientific << std::setprecision(6);
        line << "Eid, Sid, m, s, age, len, W";
        for (int i = 1; i <= Nyears; ++i)
            line << ", urf" << i;
        line << "\n";
        write_line(line.str());

        for (const auto& procData : AllProcData){
            int idx = 0;
            int Nw = static_cast<int>(procData[idx++]);
            for (int j = 0; j < Nw; ++j){
                int eid = static_cast<int>(procData[idx++]);
                int Ns  = static_cast<int>(procData[idx++]);
                for (int k = 0; k < Ns; ++k){
                    int sid     = static_cast<int>(procData[idx++]);
                    double m        = procData[idx++];
                    double s        = procData[idx++];
                    double age      = procData[idx++];
                    double len      = procData[idx++];
                    double w        = procData[idx++];

                    std::ostringstream row;
                    row << std::scientific << std::setprecision(6);
                    row << eid << ", " << sid << ", " << m << ", " << s << ", "
                        << age << ", " << len << ", " << w;

                    for (int iyr = 0; iyr < Nyears; ++iyr) {
                        row << ", " << procData[idx++];
                    }
                    row << "\n";
                    write_line(row.str());
                }
            }
        }
        // Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    void printDetailOutputFromAllProc(const std::vector<std::vector<double>>& AllProcData, const std::string& filename,
                                      int Nyears, bool compress = false){

        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // Open output stream
        if (use_compression){
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else{
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        // Set common formatting
        std::ostringstream line;
        line << std::scientific << std::setprecision(6);
        line << "Eid, Sid, UrfI, UrfJ, hru, hru_idx, UnsatD, UnsatR, SrfPrc, RivInfl";
        for (int i = 1; i <= Nyears; ++i)
            line << ", LF" << i;
        line << "\n";
        write_line(line.str());

        for (const auto& procData : AllProcData){
            int idx = 0;
            int Nw = static_cast<int>(procData[idx++]);
            for (int j = 0; j < Nw; ++j){
                int eid = static_cast<int>(procData[idx++]);
                int Ns  = static_cast<int>(procData[idx++]);
                for (int k = 0; k < Ns; ++k){
                    int sid     = static_cast<int>(procData[idx++]);
                    int urfI    = static_cast<int>(procData[idx++]);
                    int urfJ    = static_cast<int>(procData[idx++]);
                    int hru     = static_cast<int>(procData[idx++]);
                    int hru_idx = static_cast<int>(procData[idx++]);

                    double d        = procData[idx++];
                    double r        = procData[idx++];
                    double sp       = procData[idx++];
                    double rivinfl  = procData[idx++];

                    std::ostringstream row;
                    row << std::scientific << std::setprecision(6);
                    row << eid << ", " << sid << ", " << urfI << ", " << urfJ << ", "
                        << hru << ", " << hru_idx << ", "
                        << d << ", " << r << ", " << sp << ", " << rivinfl;

                    for (int iyr = 0; iyr < Nyears; ++iyr) {
                        row << ", " << procData[idx++];
                    }
                    row << "\n";
                    write_line(row.str());
                }
            }
        }
        // Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    void printHRUMassRemoved(const std::map<int,int> &hru_idx_map, const Matrix<double> &mass_removed,
                            const std::string &filename, bool compress = false) {

        const bool use_compression = compress || has_gz_extension(filename);

        std::function<void(const std::string&)> write_line;
#ifdef _USEZLIB
        gzFile gz_out = nullptr;
#endif
        std::ofstream out_file;

        // --- Open output stream
        if (use_compression) {
#ifdef _USEZLIB
            gz_out = gzopen(filename.c_str(), "wb");
            if (!gz_out) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                return;
            }
            write_line = [&gz_out](const std::string& line){
                gzputs(gz_out, line.c_str());
            };
#else
            std::cerr << "Compression requested, but zlib support is not compiled in." << std::endl;
            return;
#endif
        }
        else {
            out_file.open(filename);
            if (!out_file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return;
            }
            write_line = [&out_file](const std::string& line) {
                out_file << line;
            };
        }

        const int Nyears = mass_removed.num_cols();
        // --- Header
        std::ostringstream header;
        header << "hru";
        for (int i = 1; i <= Nyears; ++i)
            header << ", mrmv" << i;
        header << "\n";
        write_line(header.str());
        // --- Data rows
        for (const auto &kv : hru_idx_map){
            const int hru     = kv.first;
            const int hru_idx = kv.second;

            std::ostringstream row;
            row << std::scientific << std::setprecision(6);
            row << hru;

            for (int iyr = 0; iyr < Nyears; ++iyr){
                row << ", " << mass_removed(hru_idx, iyr);
            }
            row << "\n";
            write_line(row.str());
        }

        // --- Close output
#ifdef _USEZLIB
        if (gz_out) gzclose(gz_out);
#endif
        if (out_file.is_open()) out_file.close();
    }

    inline bool readSelectedWells(const std::string &filename, const std::string &groupFilename,
                                  /*std::vector<int> &idsVI, std::vector<int> &idsVD,*/
                                  std::map<int,MS::SelectedWellsGroup> &SWGmap,
                                  boost::mpi::communicator &world){

        bool tf = readSelectedWellsGroupInfo(groupFilename, SWGmap, world);
        if (!tf){
            return false;
        }
        Matrix<int> M;
        tf = RootReadsMatrixFileDistrib<int>(filename, M,3, world,100000, 20000);
        if (!tf){
            return false;
        }
        if (PrintMatrices){
            printMatrixForAllProc<int>(M, world, 910, 917, 0, 3);
        }
        //bool tf = readMatrix<int>(filename,T,nSelectedWells,2);

        std::map<int,MS::SelectedWellsGroup>::iterator it;
        for (std::size_t i = 0; i < M.num_rows(); ++i){
            it = SWGmap.find(M(i,2));
            if (it == SWGmap.end()) {
                if (world.rank() == 0) {
                    std::cout << "The group id " << M(i,2)
                              << " of the well id " << M(i,1)
                              << " of type " << M(i,0)
                              << " is not listed in the Groups file" << std::endl;
                }
                return false;
            }

            if (M(i,0) == 1){
                //idsVI.push_back(T[i][1]);
                it->second.idVI.push_back(M(i,1));
            }
            else if (M(i,0) == 2){
                //idsVD.push_back(T[i][1]);
                it->second.idVD.push_back(M(i,1));
            }
            else {
                if (world.rank() == 0) {
                    std::cout << "Invalid well type " << M(i,0)
                              << " for well id " << M(i,1)
                              << " in file " << filename << std::endl;
                }
                return false;
            }
        }
        return true;
    }
}

#endif //MANTISSA_NPSAT_DATA_H
