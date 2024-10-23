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

    struct SWAT_row{
        int hru_num = 0;
        int year = 0;
        double irrtotal_mm = 0.0;
        double irrSW_mm = 0.0;
        double irrGW_mm = 0.0;
        double irrsaltSW_kgha = 0.0;
        double irrsaltGW_kgha = 0.0;
        double totpercsalt_kgha = 0.0;
        double perc_mm = 0.0;
    };

    class SWAT_data{
    public:
        SWAT_data(){}
        bool read(const std::string& filename);


        std::vector<SWAT_row> SWAT_TAB;
        std::vector<std::vector<int>> HRUS;
        std::vector<std::vector<double>> irrtotal_mm;
        std::vector<std::vector<double>> irrSW_mm;
        std::vector<std::vector<double>> irrGW_mm;
        std::vector<std::vector<double>> irrsaltSW_kgha;
        std::vector<std::vector<double>> irrsaltGW_Kgha;
        std::vector<std::vector<double>> totpercsalt_kgha;
        std::vector<std::vector<double>> perc_mm;
    };




    bool SWAT_data::read(const std::string& filename) {
        bool ishdf = true;
        if (ishdf){
#if _USEHF > 0
            {

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
            }



#endif
        }


        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line, substr;
            double dTmp;
            int iTmp;
            std::string sTmp;
            int line_cnt = 0;
            while (getline(datafile, line)){
                std::istringstream inp(line);
                line_cnt++;
                if (line_cnt == 1)
                    continue;

                SWAT_row sr;
                getline(inp, substr,','); // Watershed
                getline(inp, substr,','); // gis_id
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.hru_num;
                }
                getline(inp, substr,','); // name
                getline(inp, substr,','); // yr
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.year;
                }
                getline(inp, substr,','); // Landuse
                getline(inp, substr,','); // Napplied_kgha
                getline(inp, substr,','); // irrsaltSW_kgha
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.irrsaltSW_kgha;
                }
                getline(inp, substr,','); // irrsaltGW_kgha
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.irrsaltGW_kgha;
                }
                getline(inp, substr,','); // irrsalt_outsidesource_kgha,
                getline(inp, substr,','); // fertsalt_kgha,
                getline(inp, substr,','); // percsalt_kgha,
                getline(inp, substr,','); // totpercsalt_kgha,
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.totpercsalt_kgha;
                }
                getline(inp, substr,','); // dssl_kgha,
                getline(inp, substr,','); // precip_mm,
                getline(inp, substr,','); // irrtotal_mm,
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.irrtotal_mm;
                }
                getline(inp, substr,','); // irrSW_mm,
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.irrSW_mm;
                }
                getline(inp, substr,','); // irrGW_mm,
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.irrGW_mm;
                }
                getline(inp, substr,','); // irrUNL_mm,
                getline(inp, substr,','); // perc_mm,
                {
                    std::istringstream inp1(substr);
                    inp1 >> sr.perc_mm;
                }
                //getline(inp, substr,','); // et_mm

                std::vector<int> ttt;
                ttt.push_back(7);

                SWAT_TAB.push_back(sr);
            }
        }
        datafile.close();
        return true;
    }

}

#endif //MANTISSA_SWAT_DATA_H
