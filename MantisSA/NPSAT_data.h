//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_NPSAT_DATA_H
#define MANTISSA_NPSAT_DATA_H

#include <vector>
#include <map>
#include <string>

#include "MS_urf_calc.h"

namespace MS{
    struct STRML{
        int urfI;
        int urfJ;
        int IJ;
        int hru_idx;
        double W;
        double Len;
        double m;
        double s;
        double a;
        std::vector<double> urf;
        std::vector<double> lf;

    };

    struct WELL{
        double initConc = 0.0;
        std::vector<STRML> strml;
    };

    struct NPSATTMP{
        double sumW = 0;
        std::vector<int> urfI;
        std::vector<int> urfJ;
        std::vector<int> hru_idx;
        std::vector<double> W;
        std::vector<double> Len;
        std::vector<double> m;
        std::vector<double> s;
        std::vector<double> a;
    };

    typedef std::map<int,NPSATTMP> tmpWELLS;

    typedef std::map<int, WELL> WELLS;

    typedef std::map<int, std::vector<int>> WELL_CELLS;

    bool readNPSATdata(std::string filename, WELLS &wells, int por, int Nyears, BackgroundRaster& braster,
                       boost::mpi::communicator &world){
#if _USEHF > 0
        {
            std::vector<std::vector<int>> ints;
            std::vector<std::vector<double>> dbls;
            std::vector<std::vector<double>> msas;
            const std::string INT_NameSet("INT");
            const std::string DBL_NameSet("DBL");
            const std::string MSA_NameSet("MSA");

            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

            HighFive::DataSet dataset_INT = HDFNfile.getDataSet(INT_NameSet);
            HighFive::DataSet dataset_DBL = HDFNfile.getDataSet(DBL_NameSet);
            HighFive::DataSet dataset_MSA = HDFNfile.getDataSet(MSA_NameSet);

            dataset_INT.read(ints);
            dataset_DBL.read(dbls);
            dataset_MSA.read(msas);

            int por_idx = 3*(por-1);
            WELLS::iterator itw = wells.end();
            std::pair<WELLS::iterator,bool> ret;

            tmpWELLS tw;
            tmpWELLS ::iterator itwtmp = tw.end();
            std::pair<tmpWELLS::iterator,bool> tmpret;
            //int count_wells = 0;
            //int wells_per_proc = Nwells/world.size();
            //int count_proc = 0;
            for (int i = 0; i < ints[0].size(); ++i){
                if (itwtmp == tw.end()){
                    itwtmp = tw.find(ints[0][i]);
                    if (itwtmp == tw.end()){
                        tmpret = tw.insert(std::pair<int,NPSATTMP>(ints[0][i], NPSATTMP()));
                            if (tmpret.second){
                                itwtmp = tmpret.first;
                            }
                    }
                }
                else{
                    if (itwtmp->first != ints[0][i]){
                        itwtmp = tw.find(ints[0][i]);
                        if (itwtmp == tw.end()){

                            tmpret = tw.insert(std::pair<int,NPSATTMP>(ints[0][i], NPSATTMP()));
                            if (tmpret.second){
                                itwtmp = tmpret.first;
                            }
                        }
                    }
                }
                itwtmp->second.urfI.push_back(ints[1][i]);
                itwtmp->second.urfJ.push_back(ints[2][i]);
                itwtmp->second.hru_idx.push_back(ints[3][i]);
                itwtmp->second.W.push_back(dbls[0][i]);
                itwtmp->second.Len.push_back(dbls[1][i]);
                itwtmp->second.m.push_back(msas[por_idx][i]);
                itwtmp->second.s.push_back(msas[por_idx+1][i]);
                itwtmp->second.a.push_back(msas[por_idx+2][i]);
                itwtmp->second.sumW = itwtmp->second.sumW + dbls[0][i];
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
                    STRML s;
                    s.urfI = itwtmp->second.urfI[i]-1;
                    s.urfJ = itwtmp->second.urfJ[i]-1;
                    s.IJ = braster.IJ(s.urfI, s.urfJ);
                    s.hru_idx = itwtmp->second.hru_idx[i];
                    s.W = itwtmp->second.W[i]/itwtmp->second.sumW;
                    s.Len = itwtmp->second.Len[i];
                    s.m = itwtmp->second.m[i];
                    s.s = itwtmp->second.s[i];
                    s.a = itwtmp->second.a[i];
                    calcURFs(Nyears,s.m, s.s, s.a, s.Len, s.urf);
                    w.strml.push_back(s);
                }

                wells.insert(std::pair<int,WELL>(itwtmp->first, w));
                count = count + 1;
            }
        }


        return true;
#endif
    }

    bool readInitSaltConc(std::string filename, WELLS &wells){
        std::ifstream datafile;
        datafile.open(filename);
        if (!datafile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            WELLS::iterator itw;
            std::string line;
            int eid;
            double conc;
            while (getline(datafile, line)){
                std::istringstream inp(line.c_str());
                inp >> eid;
                inp >> conc;
                itw = wells.find(eid);
                if (itw != wells.end()){
                    itw->second.initConc = conc;
                }
            }
            return true;
        }

    }

    bool readDistribPumping(std::string filename,  WELL_CELLS &well_cells){
#if _USEHF > 0
        const std::string NameSet("WellID");
        HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
        HighFive::DataSet dataset = HDFNfile.getDataSet(NameSet);
        std::vector<int> cellWell;
        dataset.read(cellWell);
#endif
        WELL_CELLS::iterator it;
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
}

#endif //MANTISSA_NPSAT_DATA_H
