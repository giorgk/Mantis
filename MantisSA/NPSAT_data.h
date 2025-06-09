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

#include "MS_urf_calc.h"
#include "MS_unsat.h"
#include "MS_HRU_raster.h"
#include "MSdebug.h"

namespace MS{

    bool has_gz_extension(const std::string& filename) {
        return filename.size() >= 3 &&
               filename.compare(filename.size() - 3, 3, ".gz") == 0;
    }


    bool readNPSATdata(std::string filename, WELLS &wells, int por, int Nyears,
                       BackgroundRaster& braster, RiverOptions &rivopt,
                       boost::mpi::communicator &world){
        std::string ext = getExtension(filename);
        std::vector<std::vector<int>> ints;
        std::vector<std::vector<double>> dbls;
        std::vector<std::vector<double>> msas;
        std::vector<int> version;
        std::vector<double> porosities;

        if (ext.compare("h5") == 0) {
#if _USEHF > 0
            if (world.rank() == 0) {
                std::cout << "Reading " << filename << std::endl;
            }
            const std::string INT_NameSet("INT");
            const std::string DBL_NameSet("DBL");
            const std::string MSA_NameSet("MSA");
            const std::string VER_NameSet("Version");
            const std::string POR_NameSet("Por");

            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

            HighFive::DataSet dataset_INT = HDFNfile.getDataSet(INT_NameSet);
            HighFive::DataSet dataset_DBL = HDFNfile.getDataSet(DBL_NameSet);
            HighFive::DataSet dataset_MSA = HDFNfile.getDataSet(MSA_NameSet);
            HighFive::DataSet dataset_VER = HDFNfile.getDataSet(VER_NameSet);
            HighFive::DataSet dataset_POR = HDFNfile.getDataSet(POR_NameSet);

            dataset_INT.read(ints);
            dataset_DBL.read(dbls);
            dataset_MSA.read(msas);
            dataset_VER.read(version);
            dataset_POR.read(porosities);
#endif
        }
        else {
            std::vector<int> info;
            bool tf = RootReadsNPSATinfo(filename + "info.dat", info, world);
            if (!tf){return false;}
            if (info.size() < 2){
                std::cout << filename + "info.dat should contain version, porosity and the int, dbl, msa matrices size" << std::endl;
                return false;
            }
            int idx = info[1] + 2;
            if (info.size() < idx+5){
                std::cout << filename + "info.dat should contain version, porosity and the int, dbl, msa matrices size" << std::endl;
                return false;
            }
            version.push_back(info[0]);
            for (int i = 0; i < info[1]; ++i){
                porosities.push_back(info[i + 2]);
            }
            std::vector<int> int_size;
            std::vector<int> dbl_size;
            std::vector<int> msa_size;
            int_size.push_back(info[idx]);
            int_size.push_back(info[idx+1]);
            dbl_size.push_back(info[idx+2]);
            dbl_size.push_back(info[idx+3]);
            msa_size.push_back(info[idx+4]);
            msa_size.push_back(info[idx+5]);

            tf = RootReadsMatrixFileDistrib<int>(filename + "INT.dat", ints, int_size[1], false, world, 500000);
            if (!tf){return false;}
            tf = RootReadsMatrixFileDistrib<double>(filename + "DBL.dat", dbls, dbl_size[1], false, world, 500000);
            if (!tf){return false;}
            tf = RootReadsMatrixFileDistrib<double>(filename + "MSA.dat", msas, msa_size[1], false, world, 500000);
            if (!tf){return false;}

            if (PrintMatrices){
                printMatrixForAllProc<int>(ints, world, 0, int_size[1], 0, 10);
                printMatrixForAllProc<double>(dbls, world, 0, dbl_size[1], 0, 10);
                printMatrixForAllProc<double>(msas, world, 0, msa_size[1], 0, 10);
            }
            //return false;

        }

        if (world.rank() == 0){
            std::cout << "Assembling NPSAT urf data ..." << std::endl;
        }
        wells.version = version[0];
        int por_idx = 0;
        {// find porosity index
            for (int i = 0; i < porosities.size(); ++i){
                if (porosities[i] == por){
                    por_idx = i;
                    break;
                }
            }
        }

        por_idx = 3*(por_idx-1);
        std::map<int, WELL>::iterator itw = wells.wells.end();
        std::pair<std::map<int, WELL>::iterator,bool> ret;

        tmpWELLS tw;
        tmpWELLS ::iterator itwtmp = tw.end();
        std::pair<tmpWELLS::iterator,bool> tmpret;
        //int count_wells = 0;
        //int wells_per_proc = Nwells/world.size();
        //int count_proc = 0;
        for (int i = 0; i < ints.size(); ++i){
            if (itwtmp == tw.end()){
                itwtmp = tw.find(ints[i][0]);
                if (itwtmp == tw.end()){
                    tmpret = tw.insert(std::pair<int,NPSATTMP>(ints[i][0], NPSATTMP()));
                    if (tmpret.second){
                        itwtmp = tmpret.first;
                    }
                }
            }
            else{
                if (itwtmp->first != ints[i][0]){
                    itwtmp = tw.find(ints[i][0]);
                    if (itwtmp == tw.end()){

                        tmpret = tw.insert(std::pair<int,NPSATTMP>(ints[i][0], NPSATTMP()));
                        if (tmpret.second){
                            itwtmp = tmpret.first;
                        }
                    }
                }
            }
            itwtmp->second.Sid.push_back(ints[i][1]);
            itwtmp->second.urfI.push_back(ints[i][2]);
            itwtmp->second.urfJ.push_back(ints[i][3]);
            if (version[0] == 0){
                itwtmp->second.inRiv.push_back(ints[i][4] == 1);
            }
            //itwtmp->second.hru_idx.push_back(ints[i][4]);
            itwtmp->second.W.push_back(dbls[i][0]);
            itwtmp->second.Len.push_back(dbls[i][1]);
            if (version[0] == 1){
                itwtmp->second.rivDist.push_back(dbls[i][2]);
            }
            itwtmp->second.m.push_back(msas[i][por_idx]);
            itwtmp->second.s.push_back(msas[i][por_idx+1]);
            itwtmp->second.a.push_back(msas[i][por_idx+2]);
            itwtmp->second.sumW = itwtmp->second.sumW + dbls[i][0];
        }

        //Normalize streamlines and split them according to the processors
        int count = 0;
        for (itwtmp = tw.begin(); itwtmp != tw.end(); ++itwtmp){
            if (count % world.size() != world.rank()){
                count = count + 1;
                continue;
            }
            WELL w;
            for (unsigned int i = 0; i < itwtmp->second.urfJ.size(); ++i){
                //if (itwtmp->first == 21539){
                //    bool stop = true;
                //}

                STRML s;
                s.Sid = itwtmp->second.Sid[i];
                s.urfI = itwtmp->second.urfI[i]-1;
                s.urfJ = itwtmp->second.urfJ[i]-1;
                s.IJ = braster.IJ(s.urfI, s.urfJ);
                if (version[0] == 0){
                    s.inRiv = itwtmp->second.inRiv[i];
                }
                if (version[0] == 1){
                    s.rivInfl = calcRiverInfluence(itwtmp->second.rivDist[i], rivopt);
                    s.inRiv = false;
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

            wells.wells.insert(std::pair<int,WELL>(itwtmp->first, w));
            count = count + 1;
            //std::cout << count << std::endl;
        }


        return true;

    }

    bool readInitSaltConc(std::string filename, WELLS &wells, int rank){
        std::ifstream datafile;
        datafile.open(filename);
        if (!datafile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            if (rank == 0){
                std::cout << "Reading " << filename << std::endl;
            }
            std::map<int, WELL>::iterator itw;
            std::string line;
            int eid;
            double conc;
            while (getline(datafile, line)){
                std::istringstream inp(line.c_str());
                inp >> eid;
                inp >> conc;
                itw = wells.wells.find(eid);
                if (itw != wells.wells.end()){
                    itw->second.initConc = conc;
                }
            }
            return true;
        }

    }

    bool readDistribPumping(std::string filename,  WELL_CELLS &well_cells){
        std::string ext = getExtension(filename);
        std::vector<int> cellWell;
        WELL_CELLS::iterator it;
        if (ext.compare("h5") == 0){
#if _USEHF > 0
            const std::string NameSet("WellID");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet dataset = HDFNfile.getDataSet(NameSet);
            dataset.read(cellWell);
#endif
        }
        else{
            bool tf = readVector<int>(filename, cellWell);
            if (!tf){
                return false;
            }
        }

        int step = cellWell.size()/10;
        int countSteps = 1;
        for (int i = 0; i < cellWell.size(); ++i){
            if (i > countSteps*step){
                std::cout << countSteps*10 << "% " << std::flush;
                countSteps++;
            }

            it = well_cells.find(cellWell[i]);
            if (it == well_cells.end()){
                std::vector<int> tmp(1,i);
                well_cells.insert(std::pair<int, std::vector<int>>(cellWell[i], tmp));
            }
            else{
                it->second.push_back(i);
            }
        }
        std::cout << std::endl;
        return true;
    }

    void linearizeWellBTCs(WELLS &W, std::vector<double> &v){
        v.clear();
        v.push_back(static_cast<double>(W.wells.size()));
        std::map<int, WELL>::iterator itw;
        for (itw = W.wells.begin(); itw != W.wells.end(); ++itw){
            v.push_back(static_cast<double>(itw->first));
            for (unsigned int i = 0; i < itw->second.wellBtc.size(); ++i){
                v.push_back(itw->second.wellBtc[i]);
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

    bool readSelectedWells(std::string filename, std::string groupFilename,
                           /*std::vector<int> &idsVI, std::vector<int> &idsVD,*/
                           std::map<int,MS::SelectedWellsGroup> &SWGmap,
                           boost::mpi::communicator &world){

        bool tf = readSelectedWellsGroupInfo(groupFilename, SWGmap);
        if (!tf){
            return false;
        }
        std::vector<std::vector<int>> T;
        tf = RootReadsMatrixFileDistrib<int>(filename,T,3, false,world);
        if (!tf){
            return false;
        }
        if (PrintMatrices){
            printMatrixForAllProc<int>(T, world, 910, 917, 0, 3);
        }
        //bool tf = readMatrix<int>(filename,T,nSelectedWells,2);

        std::map<int,MS::SelectedWellsGroup>::iterator it;
        for (unsigned int i = 0; i < T.size(); ++i){
            it = SWGmap.find(T[i][2]);
            if (it != SWGmap.end()){
                if (T[i][0] == 1){
                    //idsVI.push_back(T[i][1]);
                    it->second.idVI.push_back(T[i][1]);
                }
                else if (T[i][0] == 2){
                    //idsVD.push_back(T[i][1]);
                    it->second.idVD.push_back(T[i][1]);
                }
            }
            else{
                std::cout << "The group id " << T[i][2] << "of the well id " << T[i][1] << " of type "  << T[i][0] << " is not listed in the Groups file" << std::endl;
                return false;
            }
        }
        return true;
    }
}

#endif //MANTISSA_NPSAT_DATA_H
