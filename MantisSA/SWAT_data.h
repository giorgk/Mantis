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
        bool read(const std::string filename, int NSwatYears);


        //std::vector<SWAT_row> SWAT_TAB;
        std::vector<std::vector<int>> HRUS;
        std::vector<std::vector<double>> irrtotal_mm;
        std::vector<std::vector<double>> irrSW_mm;
        std::vector<std::vector<double>> irrGW_mm;
        std::vector<std::vector<double>> irrsaltSW_kgha;
        std::vector<std::vector<double>> irrsaltGW_Kgha;
        std::vector<std::vector<double>> totpercsalt_kgha;
        std::vector<std::vector<double>> perc_mm;
    private:
        bool readASCIIset(std::string filename, std::vector<std::vector<double>> &data, int nHRU, int Nyrears);
        bool readASCIIHRUS(std::string filename, std::vector<std::vector<int>> &data);
    };




    bool SWAT_data::read(const std::string filename, int NSwatYears) {

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

                HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);

                HighFive::DataSet dataset_HRUS = HDFNfile.getDataSet(HRUSNameSet);
                HighFive::DataSet dataset_irrtotal_mm = HDFNfile.getDataSet(irrtotal_mm_NameSet);
                HighFive::DataSet dataset_irrSW_mm = HDFNfile.getDataSet(irrSW_mm_NameSet);
                HighFive::DataSet dataset_irrGW_mm = HDFNfile.getDataSet(irrGW_mm_NameSet);
                HighFive::DataSet dataset_irrsaltSW_kgha = HDFNfile.getDataSet(irrsaltSW_kgha_NameSet);
                HighFive::DataSet dataset_irrsaltGW_Kgha = HDFNfile.getDataSet(irrsaltGW_Kgha_NameSet);
                HighFive::DataSet dataset_totpercsalt_kgha = HDFNfile.getDataSet(totpercsalt_kgha_NameSet);
                HighFive::DataSet dataset_perc_mm = HDFNfile.getDataSet(perc_mm_NameSet);

                dataset_HRUS.read(HRUS);
                dataset_irrtotal_mm.read(irrtotal_mm);
                dataset_irrSW_mm.read(irrSW_mm);
                dataset_irrGW_mm.read(irrGW_mm);
                dataset_irrsaltSW_kgha.read(irrsaltSW_kgha);
                dataset_irrsaltGW_Kgha.read(irrsaltGW_Kgha);
                dataset_totpercsalt_kgha.read(totpercsalt_kgha);
                dataset_perc_mm.read(perc_mm);

                return true;
#endif
        }
        else{
            int nHRUs;
            {// Read HRUS;
                bool tf = readASCIIHRUS(filename + "hrus.dat",HRUS);
                if (!tf){
                    return false;
                }
                nHRUs = HRUS[0].size();
            }
            { //Read irrtotal_mm
                bool tf = readASCIIset(filename + "irrtotal_mm.dat", irrtotal_mm,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read irrSW_mm
                bool tf = readASCIIset(filename + "irrSW_mm.dat", irrSW_mm,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read irrGW_mm
                bool tf = readASCIIset(filename + "irrGW_mm.dat", irrGW_mm,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read irrsaltSW_kgha
                bool tf = readASCIIset(filename + "irrsaltSW_kgha.dat", irrsaltSW_kgha,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read irrsaltGW_Kgha
                bool tf = readASCIIset(filename + "irrsaltGW_Kgha.dat", irrsaltGW_Kgha,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read totpercsalt_kgha
                bool tf = readASCIIset(filename + "totpercsalt_kgha.dat", totpercsalt_kgha,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            { //Read perc_mm
                bool tf = readASCIIset(filename + "perc_mm.dat", perc_mm,nHRUs,NSwatYears);
                if (!tf){
                    return false;
                }
            }
            return true;
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
