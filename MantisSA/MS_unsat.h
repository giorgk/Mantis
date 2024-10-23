//
// Created by giorg on 10/22/2024.
//

#ifndef MANTISSA_MS_UNSAT_H
#define MANTISSA_MS_UNSAT_H
namespace MS{

    class UNSAT{
    public:
        UNSAT(){}
        bool readdata(std::string dpth_file, std::string dpth_name, std::string rch_file,
                      double wc_in, double d_in, double r_in);

        int traveltime(int IJ);
    private:
        double wc;
        double minD;
        double minR;
        std::vector<double> Depth;
        std::vector<double> Rch;
        int Ncells;
    };

    bool UNSAT::readdata(std::string dpth_file, std::string dpth_name, std::string rch_file,
                         double wc_in, double d_in, double r_in) {
        wc = wc_in;
        minD = d_in;
        minR = r_in;


#if _USEHF>0
        {   // Read the Depth file
            const std::string NamesNameSet("Names");
            const std::string DataNameSet("Data");
            HighFive::File HDFfile(dpth_file, HighFive::File::ReadOnly);
            HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
            HighFive::DataSet datasetData = HDFfile.getDataSet(DataNameSet);

            std::vector<std::string> names;
            std::vector<std::vector<double>> tmp;
            datasetNames.read(names);
            datasetData.read(tmp);
            int idx = 0;
            bool index_found = false;
            for (int i = 0; i < static_cast<int>(names.size()); ++i){
                if (dpth_name == names[i]){
                    idx = i;
                    index_found = true;
                    continue;
                }
            }
            if (index_found){
                Depth = tmp[idx];
                Ncells = static_cast<int>(Depth.size());
            }
            else{
                std::cout << " The depth scenario << " << dpth_name << " was not found in the list of scenarios in the " << dpth_file << std::endl;
                return false;
            }
        }

        {   // Read Recharge
            const std::string NamesNameSet("Names");
            const std::string DataNameSet("Data");
            HighFive::File HDFfile(rch_file, HighFive::File::ReadOnly);
            HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
            HighFive::DataSet datasetData = HDFfile.getDataSet(DataNameSet);

            std::vector<std::string> names;
            std::vector<std::vector<double>> tmp;
            datasetNames.read(names);
            datasetData.read(tmp);

            Rch = tmp[0];
            if (Ncells != Rch.size()){
                std::cout << "The size of Depth is different than the size of Recharge" << std::endl;
                return false;
            }
        }
#endif


        return true;
    }

    int UNSAT::traveltime(int IJ) {
        if (IJ < 0 || IJ >= Ncells ){
            return 0;
        }
        double d = Depth[IJ];
        if (d < minD){
            d = minD;
        }

        double r = Rch[IJ];
        if (r < minR){
            r = minR;
        }

        double tau = wc * d / (r/1000);
        return static_cast<int>(std::ceil(tau));
    }

}

#endif //MANTISSA_MS_UNSAT_H
