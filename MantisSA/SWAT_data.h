//
// Created by giorg on 10/16/2024.
//

#ifndef MANTISSA_SWAT_DATA_H
#define MANTISSA_SWAT_DATA_H

#include <fstream>
#include <vector>

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

#include "MS_mpi_utils.h"

namespace MS {

    /*struct SWAT_row{
        int hru_num = 0;
        int year = 0;
        double irrtotal_mm = 0.0;
        double irrSW_mm = 0.0;
        double irrGW_mm = 0.0;
        double irrsaltSW_kgha = 0.0;
        double irrsaltGW_kgha = 0.0;
        double totpercsalt_kgha = 0.0;
        double perc_mm = 0.0;
    };*/

    class SWAT_data{
    public:
        SWAT_data(){}
        bool read(const std::string filename, int NSwatYears, int swat_version, boost::mpi::communicator &world);
        bool read_v0(const std::string filename, int NSwatYears, boost::mpi::communicator &world);
        bool read_v1(const std::string filename, int NSwatYears, boost::mpi::communicator &world);
        bool read_HRU_idx_Map(std::string filename, boost::mpi::communicator &world);
        int hru_index(int HRU);


        //std::vector<SWAT_row> SWAT_TAB;
        //std::vector<std::vector<int>> HRUS;
        std::vector<std::vector<double>> irrtotal_mm;
        std::vector<std::vector<double>> irrSW_mm;
        std::vector<std::vector<double>> irrGW_mm;
        std::vector<std::vector<double>> irrsaltSW_kgha;
        std::vector<std::vector<double>> irrsaltGW_Kgha;
        std::vector<std::vector<double>> fertsalt_kgha;
        std::vector<std::vector<double>> dssl_kgha;
        std::vector<std::vector<double>> uptk_kgha;
        std::vector<std::vector<double>> pGW;
        std::vector<std::vector<double>> totpercsalt_kgha;
        std::vector<std::vector<double>> perc_mm;
        std::vector<std::vector<double>> Salt_perc_ppm;
        std::map<int,int> hru_idx_map;
    private:
        bool readASCIIset(std::string filename, std::vector<std::vector<double>> &data, int nHRU, int Nyrears);
        bool readASCIIHRUS(std::string filename, std::vector<std::vector<int>> &data);
    };


    bool SWAT_data::read_v0(const std::string filename, int NSwatYears, boost::mpi::communicator &world) {

        std::string ext = getExtension(filename);

        if (ext.compare("h5") == 0){
#if _USEHF > 0
                const std::string HRUSNameSet("HRUS");
                const std::string irrtotal_mm_NameSet("irrtotal_mm");
                const std::string irrSW_mm_NameSet("irrSW_mm");
                const std::string irrGW_mm_NameSet("irrGW_mm");
                const std::string irrsaltSW_kgha_NameSet("irrsaltSW_kgha");
                const std::string irrsaltGW_Kgha_NameSet("irrsaltGW_Kgha");
                const std::string totpercsalt_kgha_NameSet("totpercsalt_kgha");
                const std::string perc_mm_NameSet("perc_mm");
                const std::string Salt_perc_ppm_NameSet("Salt_perc_ppm");

                HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

                //HighFive::DataSet dataset_HRUS = HDFNfile.getDataSet(HRUSNameSet);
                HighFive::DataSet dataset_irrtotal_mm = HDFNfile.getDataSet(irrtotal_mm_NameSet);
                HighFive::DataSet dataset_irrSW_mm = HDFNfile.getDataSet(irrSW_mm_NameSet);
                HighFive::DataSet dataset_irrGW_mm = HDFNfile.getDataSet(irrGW_mm_NameSet);
                HighFive::DataSet dataset_irrsaltSW_kgha = HDFNfile.getDataSet(irrsaltSW_kgha_NameSet);
                HighFive::DataSet dataset_irrsaltGW_Kgha = HDFNfile.getDataSet(irrsaltGW_Kgha_NameSet);
                HighFive::DataSet dataset_totpercsalt_kgha = HDFNfile.getDataSet(totpercsalt_kgha_NameSet);
                HighFive::DataSet dataset_perc_mm = HDFNfile.getDataSet(perc_mm_NameSet);
                HighFive::DataSet dataset_Salt_perc_ppm = HDFNfile.getDataSet(Salt_perc_ppm_NameSet);

                //dataset_HRUS.read(HRUS);
                dataset_irrtotal_mm.read(irrtotal_mm);
                dataset_irrSW_mm.read(irrSW_mm);
                dataset_irrGW_mm.read(irrGW_mm);
                dataset_irrsaltSW_kgha.read(irrsaltSW_kgha);
                dataset_irrsaltGW_Kgha.read(irrsaltGW_Kgha);
                dataset_totpercsalt_kgha.read(totpercsalt_kgha);
                dataset_perc_mm.read(perc_mm);
                dataset_Salt_perc_ppm.read(Salt_perc_ppm);

                return true;
#endif
        }
        else{
            //bool tf = RootReadsMatrixFileDistrib<int>(filename + "hrus.dat", HRUS,1,false,world);
            //if (!tf){ return false;}
            //int nHRUs = HRUS[0].size();

            bool tf = RootReadsMatrixFileDistrib<double>(filename + "irrtotal_mm.dat", irrtotal_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrtotal_mm,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "irrSW_mm.dat", irrSW_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrSW_mm,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "irrGW_mm.dat", irrGW_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrGW_mm,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltSW_kgha.dat", irrsaltSW_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrsaltSW_kgha,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltGW_Kgha.dat", irrsaltGW_Kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrsaltGW_Kgha,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "totpercsalt_kgha.dat", totpercsalt_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(totpercsalt_kgha,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "perc_mm.dat", perc_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(perc_mm,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "Salt_perc_ppm.dat", Salt_perc_ppm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(Salt_perc_ppm,world,0,10,34,40);}

            return true;
        }
    }

    bool SWAT_data::read_v1(const std::string filename, int NSwatYears, boost::mpi::communicator &world) {
        std::string ext = getExtension(filename);

        if (ext.compare("h5") == 0){
#if _USEHF > 0
            //const std::string HRUSNameSet("HRUS");
            //const std::string irrtotal_mm_NameSet("irrtotal_mm");
            const std::string irrSW_mm_NameSet("irrSW_mm");
            const std::string irrGW_mm_NameSet("irrGW_mm");
            const std::string irrsaltSW_kgha_NameSet("irrsaltSW_kgha");
            const std::string irrsaltGW_Kgha_NameSet("irrsaltGW_Kgha");
            const std::string fertsalt_kgha_NameSet("fertsalt_kgha");
            const std::string dssl_kgha_NameSet("dssl_kgha");
            const std::string uptk_kgha_NameSet("uptk_kgha");
            const std::string pGW_NameSet("pGW");
            //const std::string totpercsalt_kgha_NameSet("totpercsalt_kgha");
            const std::string perc_mm_NameSet("perc_mm");
            const std::string Salt_perc_ppm_NameSet("Salt_perc_ppm");

            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

            //HighFive::DataSet dataset_HRUS = HDFNfile.getDataSet(HRUSNameSet);
            //HighFive::DataSet dataset_irrtotal_mm = HDFNfile.getDataSet(irrtotal_mm_NameSet);
            HighFive::DataSet dataset_irrSW_mm = HDFNfile.getDataSet(irrSW_mm_NameSet);
            HighFive::DataSet dataset_irrGW_mm = HDFNfile.getDataSet(irrGW_mm_NameSet);
            HighFive::DataSet dataset_irrsaltSW_kgha = HDFNfile.getDataSet(irrsaltSW_kgha_NameSet);
            HighFive::DataSet dataset_irrsaltGW_Kgha = HDFNfile.getDataSet(irrsaltGW_Kgha_NameSet);
            HighFive::DataSet dataset_fertsalt_kgha = HDFNfile.getDataSet(fertsalt_kgha_NameSet);
            HighFive::DataSet dataset_dssl_kgha = HDFNfile.getDataSet(dssl_kgha_NameSet);
            HighFive::DataSet dataset_uptk_kgha = HDFNfile.getDataSet(uptk_kgha_NameSet);
            HighFive::DataSet dataset_pGW = HDFNfile.getDataSet(pGW_NameSet);
            //HighFive::DataSet dataset_totpercsalt_kgha = HDFNfile.getDataSet(totpercsalt_kgha_NameSet);
            HighFive::DataSet dataset_perc_mm = HDFNfile.getDataSet(perc_mm_NameSet);
            HighFive::DataSet dataset_Salt_perc_ppm = HDFNfile.getDataSet(Salt_perc_ppm_NameSet);

            //dataset_HRUS.read(HRUS);
            //dataset_irrtotal_mm.read(irrtotal_mm);
            dataset_irrSW_mm.read(irrSW_mm);
            dataset_irrGW_mm.read(irrGW_mm);
            dataset_irrsaltSW_kgha.read(irrsaltSW_kgha);
            dataset_irrsaltGW_Kgha.read(irrsaltGW_Kgha);
            dataset_fertsalt_kgha.read(fertsalt_kgha);
            dataset_dssl_kgha.read(dssl_kgha);
            dataset_uptk_kgha.read(uptk_kgha);
            dataset_pGW.read(pGW);
            //dataset_totpercsalt_kgha.read(totpercsalt_kgha);
            dataset_perc_mm.read(perc_mm);
            dataset_Salt_perc_ppm.read(Salt_perc_ppm);

            return true;

#endif
        }
        else{
            bool tf = false;
            //bool tf = RootReadsMatrixFileDistrib<int>(filename + "hrus.dat", HRUS,1,false,world);
            //if (!tf){ return false;}
            //int nHRUs = HRUS[0].size();

            //bool tf = RootReadsMatrixFileDistrib<double>(filename + "irrtotal_mm.dat", irrtotal_mm, NSwatYears, true, world);
            //if (!tf){ return false;}
            //if (PrintMatrices){printMatrixForAllProc(irrtotal_mm,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "irrSW_mm.dat", irrSW_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrSW_mm,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "irrGW_mm.dat", irrGW_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrGW_mm,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltSW_kgha.dat", irrsaltSW_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrsaltSW_kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltGW_Kgha.dat", irrsaltGW_Kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(irrsaltGW_Kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "fertsalt_kgha.dat", fertsalt_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(fertsalt_kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "dssl_kgha.dat", dssl_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(dssl_kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "uptk_kgha.dat", uptk_kgha, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(uptk_kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "pGW.dat", pGW, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(pGW,world,0,10,34,40);}

            //tf = RootReadsMatrixFileDistrib<double>(filename + "totpercsalt_kgha.dat", totpercsalt_kgha, NSwatYears, true, world);
            //if (!tf){ return false;}
            //if (PrintMatrices){printMatrixForAllProc(totpercsalt_kgha,world,0,10,34,40);}

            tf = RootReadsMatrixFileDistrib<double>(filename + "perc_mm.dat", perc_mm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(perc_mm,world,0,10,34,40);}
            tf = RootReadsMatrixFileDistrib<double>(filename + "Salt_perc_ppm.dat", Salt_perc_ppm, NSwatYears, true, world);
            if (!tf){ return false;}
            if (PrintMatrices){printMatrixForAllProc(Salt_perc_ppm,world,0,10,34,40);}

            return true;
        }
    }

    bool SWAT_data::read(const std::string filename, int NSwatYears, int swat_version,
                         boost::mpi::communicator &world) {
        bool tf = false;
        if (swat_version == 0){
            tf = read_v0(filename, NSwatYears,world);
        }
        else if(swat_version == 1){
            tf = read_v1(filename, NSwatYears, world);
        }
        return tf;
    }

    bool SWAT_data::read_HRU_idx_Map(std::string filename, boost::mpi::communicator &world){
        std::vector<std::vector<int>> tmp;
        bool tf = RootReadsMatrixFileDistrib<int>(filename, tmp, 1, true, world);
        if (!tf){ return false;}
        if (PrintMatrices){printMatrixForAllProc(tmp,world,0,1,0,10);}

        for (int i = 0; i < tmp[0].size(); ++i){
            hru_idx_map.insert(std::pair<int,int>(tmp[0][i], i));
        }
        return true;
    }

    int SWAT_data::hru_index(int HRU){
        if (HRU < 0){
            return -9;
        }
        std::map<int,int>::iterator it = hru_idx_map.find(HRU);
        if (it != hru_idx_map.end()){
            return it->second;
        }
        else{
            return -9;
        }

    }

    bool SWAT_data::readASCIIset(std::string filename, std::vector<std::vector<double>> &data, int nHRU, int Nyrears){
        data.clear();
        data.resize(Nyrears,std::vector<double>(nHRU,0));

        std::ifstream ifile;
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            double val;
            for (int i = 0; i < nHRU; ++i){
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                for (int j = 0; j < Nyrears; ++j){
                    inp >> val;
                    data[j][i] = val;
                }
            }
            ifile.close();
            return true;
        }
    }

    bool SWAT_data::readASCIIHRUS(std::string filename, std::vector<std::vector<int>> &data){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            std::vector<int> tmp;
            int v;
            while (getline(datafile, line)){
                std::istringstream inp(line.c_str());
                inp >> v;
                tmp.push_back(v);
            }
            data.push_back(tmp);
            datafile.close();
        }
        return true;
    }
}

#endif //MANTISSA_SWAT_DATA_H
