//
// Created by giorg on 10/29/2024.
//

#ifndef MANTISSA_MSDEBUG_H
#define MANTISSA_MSDEBUG_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "MS_structures.h"
#include "MS_urf_calc.h"

namespace MS{

    // CHANGE (v0.1.15): all debug levels use one fixed set of columns. Empty
    // columns do not apply to that record. This makes one CSV easy to filter by
    // level, year, streamline, cell, HRU, or any calculation input/output.
    enum DebugColumn {
        DBG_LEVEL = 0, DBG_IYR, DBG_YEAR, DBG_ISWAT, DBG_EID, DBG_SID,
        DBG_BUFFER_CELL, DBG_CELL_I, DBG_CELL_J, DBG_IJ,
        DBG_HRU, DBG_HRUIDX, DBG_REGIME, DBG_STRM_W, DBG_INIT_CONC,
        DBG_TREATED_CONC_FROM_PUMP, DBG_SHIFT, DBG_SHIFTED_IYR,
        DBG_BLEND_U, DBG_SURF_PERC, DBG_RIV_INFL, DBG_IN_RIVER,
        DBG_AGE, DBG_LENGTH, DBG_M, DBG_S, DBG_WELL_SOURCE_ID, DBG_IRRGW_MM,
        DBG_IRRSW_MM, DBG_IRRSALTSW_KGHA, DBG_FERTSALT_KGHA,
        DBG_DSSL_KGHA, DBG_QSALT_KGHA, DBG_UPTK_KGHA, DBG_PGW,
        DBG_M_NPSAT, DBG_M_TOTAL_RAW, DBG_M_TOTAL, DBG_PERC_MM,
        DBG_SALT_PERC_PPM, DBG_WATER_WEIGHT, DBG_WEIGHTED_C, DBG_C_HIST,
        DBG_C_FUTURE, DBG_C_BEFORE_RIVER, DBG_C_AFTER_RIVER,
        DBG_C_AFTER_SURFACE, DBG_C_CELL,
        DBG_C_GW, DBG_C_TOT, DBG_BTC, DBG_PREBTC, DBG_CONVOLVED,
        DBG_STRM_CONTRIBUTION, DBG_SUMW, DBG_WELL_CONC, DBG_COLUMN_COUNT
    };

    struct DebugCsvRow {
        std::vector<std::string> value;

        DebugCsvRow() : value(DBG_COLUMN_COUNT, "") {}

        void set(int column, const std::string &v) {
            value[column] = v;
        }

        void set(int column, int v) {
            std::ostringstream text;
            text << v;
            value[column] = text.str();
        }

        void set(int column, double v) {
            if (std::isnan(v)) {
                return;
            }
            std::ostringstream text;
            text << std::setprecision(17) << v;
            value[column] = text.str();
        }

        void write(std::ofstream &file) const {
            for (int i = 0; i < DBG_COLUMN_COUNT; ++i) {
                if (i > 0) {
                    file << ',';
                }
                file << value[i];
            }
            file << '\n';
        }
    };

    inline void writeDebugHeader(std::ofstream &file) {
        file << "level,iyr,year,iswat,Eid,Sid,buffer_cell,cell_i,cell_j,"
                "IJ,hru,hruidx,regime,strm_W,init_conc,treated_ConcFromPump,shift,"
                "shifted_iyr,blend_u,surf_perc,riv_infl,in_river,age,length,m,s,"
                "well_source_id,irrGW_mm,irrSW_mm,"
                "irrsaltSW_kgha,fertsalt_kgha,dssl_kgha,Qsalt_kgha,uptk_kgha,pGW,"
                "m_npsat,m_total_raw,m_total,perc_mm,Salt_perc_ppm,water_weight,"
                "weighted_c,c_hist,c_future,c_before_river,c_after_river,"
                "c_after_surface,c_cell,c_gw,c_tot,btc,prebtc,convolved,"
                "strm_contribution,sumW,well_conc\n";
    }

    inline void writeDebugStaticRows(std::ofstream &file, int eid,
                                     const WELL &well) {
        for (std::size_t i = 0; i < well.strml.size(); ++i) {
            const STRML &strm = well.strml[i];
            // CHANGE (v0.1.17): print one metadata row per streamline. The
            // expanded URF values and indices are intentionally not printed.
            DebugCsvRow row;
            row.set(DBG_LEVEL, "STREAMLINE_STATIC");
            row.set(DBG_EID, eid);
            row.set(DBG_SID, strm.Sid);
            row.set(DBG_IJ, strm.IJ);
            row.set(DBG_STRM_W, strm.W);
            row.set(DBG_INIT_CONC, well.initConc);
            row.set(DBG_RIV_INFL, strm.rivInfl);
            row.set(DBG_IN_RIVER, strm.inRiv ? 1 : 0);
            row.set(DBG_AGE, strm.a);
            row.set(DBG_LENGTH, strm.Len);
            row.set(DBG_M, strm.m);
            row.set(DBG_S, strm.s);
            row.set(DBG_WELL_SOURCE_ID, strm.wellSourceId);
            row.write(file);

            // OLD (pre-v0.1.17): one STREAMLINE_STATIC row was written for
            // every URF index and value.
            // for (std::size_t k = 0; k < strm.urf.size(); ++k) { ... }
        }
    }

    inline void writeDebugCellRow(std::ofstream &file, int iyr, int year,
                                  int iswat, int eid, int sid, int buffer_cell,
                                  int cell_i, int cell_j, int IJ, double strm_W,
                                  double init_conc, double riv_infl,
                                  const CellDebugData &data) {
        DebugCsvRow row;
        row.set(DBG_LEVEL, "CELL_YEAR");
        row.set(DBG_IYR, iyr);
        row.set(DBG_YEAR, year);
        row.set(DBG_ISWAT, iswat);
        row.set(DBG_EID, eid);
        row.set(DBG_SID, sid);
        row.set(DBG_BUFFER_CELL, buffer_cell);
        row.set(DBG_CELL_I, cell_i);
        row.set(DBG_CELL_J, cell_j);
        row.set(DBG_IJ, IJ);
        row.set(DBG_HRU, data.hru);
        row.set(DBG_HRUIDX, data.hruidx);
        row.set(DBG_REGIME, data.regime);
        row.set(DBG_STRM_W, strm_W);
        row.set(DBG_INIT_CONC, init_conc);
        row.set(DBG_TREATED_CONC_FROM_PUMP, data.treated_conc_from_pump);
        row.set(DBG_SHIFT, data.shift);
        row.set(DBG_SHIFTED_IYR, data.shifted_iyr);
        row.set(DBG_BLEND_U, data.blend_u);
        row.set(DBG_SURF_PERC, data.surf_perc);
        row.set(DBG_RIV_INFL, riv_infl);
        row.set(DBG_IRRGW_MM, data.irrGW_mm);
        row.set(DBG_IRRSW_MM, data.irrSW_mm);
        row.set(DBG_IRRSALTSW_KGHA, data.irrsaltSW_kgha);
        row.set(DBG_FERTSALT_KGHA, data.fertsalt_kgha);
        row.set(DBG_DSSL_KGHA, data.dssl_kgha);
        row.set(DBG_QSALT_KGHA, data.Qsalt_kgha);
        row.set(DBG_UPTK_KGHA, data.uptk_kgha);
        row.set(DBG_PGW, data.pGW);
        row.set(DBG_M_NPSAT, data.m_npsat);
        row.set(DBG_M_TOTAL_RAW, data.m_total_raw);
        row.set(DBG_M_TOTAL, data.m_total);
        row.set(DBG_PERC_MM, data.perc_mm);
        row.set(DBG_SALT_PERC_PPM, data.Salt_perc_ppm);
        row.set(DBG_WATER_WEIGHT, data.water_weight);
        row.set(DBG_WEIGHTED_C, data.c_cell * data.water_weight);
        row.set(DBG_C_HIST, data.c_hist);
        row.set(DBG_C_FUTURE, data.c_future);
        row.set(DBG_C_BEFORE_RIVER, data.c_before_river);
        row.set(DBG_C_AFTER_RIVER, data.c_after_river);
        row.set(DBG_C_AFTER_SURFACE, data.c_after_surface);
        row.set(DBG_C_CELL, data.c_cell);
        row.set(DBG_C_GW, data.c_gw);
        row.write(file);
    }

    inline void writeDebugStreamlineRow(std::ofstream &file, int iyr, int year,
                                        int iswat, int eid, const STRML &strm,
                                        double init_conc, double c_gw,
                                        double c_tot, double btc, double prebtc,
                                        double sumW, double max_age) {
        const double convolved = btc + prebtc;
        DebugCsvRow row;
        row.set(DBG_LEVEL, "STREAMLINE_YEAR");
        row.set(DBG_IYR, iyr);
        row.set(DBG_YEAR, year);
        row.set(DBG_ISWAT, iswat);
        row.set(DBG_EID, eid);
        row.set(DBG_SID, strm.Sid);
        row.set(DBG_IJ, strm.IJ);
        row.set(DBG_STRM_W, strm.W);
        row.set(DBG_INIT_CONC, init_conc);
        row.set(DBG_RIV_INFL, strm.rivInfl);
        row.set(DBG_IN_RIVER, strm.inRiv ? 1 : 0);
        row.set(DBG_AGE, strm.a);
        row.set(DBG_LENGTH, strm.Len);
        row.set(DBG_M, strm.m);
        row.set(DBG_S, strm.s);
        row.set(DBG_WELL_SOURCE_ID, strm.wellSourceId);
        if (strm.inRiv) {
            row.set(DBG_REGIME, "river");
        }
        else if (strm.a > max_age) {
            row.set(DBG_REGIME, "max_age");
        }
        else {
            row.set(DBG_REGIME, "cells");
        }
        row.set(DBG_C_GW, c_gw);
        row.set(DBG_C_TOT, c_tot);
        row.set(DBG_BTC, btc);
        row.set(DBG_PREBTC, prebtc);
        row.set(DBG_CONVOLVED, convolved);
        row.set(DBG_STRM_CONTRIBUTION, strm.W * convolved);
        row.set(DBG_SUMW, sumW);
        row.write(file);
    }

    inline void writeDebugWellRow(std::ofstream &file, int iyr, int year,
                                  int iswat, int eid, double init_conc,
                                  double sumW, double well_conc) {
        DebugCsvRow row;
        row.set(DBG_LEVEL, "WELL_YEAR");
        row.set(DBG_IYR, iyr);
        row.set(DBG_YEAR, year);
        row.set(DBG_ISWAT, iswat);
        row.set(DBG_EID, eid);
        row.set(DBG_INIT_CONC, init_conc);
        row.set(DBG_SUMW, sumW);
        row.set(DBG_WELL_CONC, well_conc);
        row.write(file);
    }

    void testConvolution(){
        std::vector<double> urf, lf, btc, prebtc;
        MS::calcURFs(200, 4.782, 0.4213, 140.32, 1074.7, urf);
        double x = 1;
        for (int i = 0; i < 100; ++i){
            lf.push_back(x);
            x = x+1.0;
        }
        x = 100;
        for (int i = 0; i < 100; ++i){
            lf.push_back(x);
            x = x - 1;
        }
        btc.resize(lf.size(),0);
        prebtc.resize(lf.size(),0);
        MS::convolute(urf,lf,btc,prebtc,100.0);
        for (int i = 0; i < urf.size(); ++i){
            std::cout << urf[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < lf.size(); ++i){
            std::cout << lf[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < btc.size(); ++i){
            std::cout << btc[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < prebtc.size(); ++i){
            std::cout << prebtc[i] << " ";
        }
        std::cout << std::endl;
    }

    void testBroadcast(boost::mpi::communicator& world){
        int N=30000000;
        std::vector<double> tmp(N,0);
        if (world.rank() == 0) {
            for (int i = 0; i < tmp.size(); ++i)
                tmp[i] = static_cast<double>(i)*2.36547;
        }
        //boost::mpi::broadcast(world, &tmp[0],N, 0);
        MS::sendVectorFromRoot2AllProc<double>(tmp, world);

        for (int k = 0; k < world.size(); ++k){
            for (int i = N-20; i < N; ++i){
                std::cout << tmp[i] << " ";
            }
            std::cout << std::endl;
        }
    }

    void printConcFromPump(std::string &fname, std::vector<double> &conc){
        std::ofstream out_file;
        out_file.open(fname.c_str());
        for (unsigned int i = 0; i < conc.size(); ++i){
            out_file << conc[i] << std::endl;
        }
        out_file.close();

    }

    template<typename T>
    void testMatrixBroadcast(boost::mpi::communicator& world){
        // Create the matrix
        std::vector<std::vector<T>> vMat;
        if (world.rank() == 0){
            int idx = 1;
            for (int i = 0; i < 7; ++i){
                std::vector<T> vVec;
                for (int j = 0; j < 10; ++j ){
                    vVec.push_back(static_cast<T>(idx));
                    idx = idx + 1;
                }
                vMat.push_back(vVec);
            }
        }
        else{
            std::vector<T> vVec;
            vVec.resize(10,0);
            vMat.resize(7,vVec);
        }
        world.barrier();

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                for (unsigned int i = 0; i < vMat.size(); ++i){
                    for (unsigned int j = 0; j < vMat[i].size(); ++j){
                        std::cout << vMat[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            world.barrier();
        }


        sentMatrixFromRoot2AllProc(vMat, world);

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                for (unsigned int i = 0; i < vMat.size(); ++i){
                    for (unsigned int j = 0; j < vMat[i].size(); ++j){
                        std::cout << vMat[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            world.barrier();
        }
    }

    void testScalarBroadcast(boost::mpi::communicator& world){
        int ia = 0;
        double ib= 0.0;
        if (world.rank() == 0){
            ia = 100;
            ib = 99.9999;
        }
        world.barrier();
        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " " << ia << " " << ib << std::endl;
            }
            world.barrier();
        }
        MS::sendScalarFromRoot2AllProc<int>(ia, world);
        MS::sendScalarFromRoot2AllProc<double>(ib, world);

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " " << ia << " " << ib << std::endl;
            }
            world.barrier();
        }

    }


}
#endif //MANTISSA_MSDEBUG_H
