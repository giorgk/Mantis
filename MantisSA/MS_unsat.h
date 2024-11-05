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
                      double wc_in, double d_in, double r_in, int nRasterCells);

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
                         double wc_in, double d_in, double r_in, int nRasterCells) {
        wc = wc_in;
        minD = d_in;
        minR = r_in;
        std::string ext = getExtension(dpth_file);
        if (ext.compare("h5") == 0) { // Read the Depth file
#if _USEHF>0
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
                    break;
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
#endif
        }
        else{
            std::ifstream ifile;
            ifile.open(dpth_file);
            if (!ifile.is_open()){
                std::cout << "Cant open file: " << dpth_file << std::endl;
                return false;
            }
            else{
                std::cout << "Reading " << dpth_file << std::endl;
                int nData;
                std::string line;
                { //Get the number of Unsaturated scenarios
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> nData;
                    Depth.clear();
                }

                int idx = 0;
                {// Get the scenario names
                    std::string name;
                    bool index_found = false;
                    for (int i = 0; i < nData; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        inp >> name;
                        if (dpth_name == name){
                            idx = i;
                            index_found = true;
                            break;
                        }
                    }
                }
                {// Read the values
                    Depth.resize(nRasterCells,0.0);
                    for (int i = 0; i < nRasterCells; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        double v;
                        for (int j = 0; j < nData; ++j){
                            inp >> v;
                            if (j == idx){
                                Depth[i] = v;
                                break;
                            }
                        }
                    }
                }
                ifile.close();
            }
        }

        std::string ext1 = getExtension(rch_file);
        if (ext1.compare("h5") == 0)
        {   // Read Recharge
#if _USEHF>0
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
#endif
        }
        else{
            std::ifstream ifile;
            ifile.open(dpth_file);
            if (!ifile.is_open()){
                std::cout << "Cant open file: " << rch_file << std::endl;
                return false;
            }
            else{
                std::cout << "Reading " << rch_file << std::endl;
                int nData;
                std::string line;
                { //Get the number of Unsaturated scenarios
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> nData;
                    Rch.clear();
                }
                {// Get the scenario names
                    std::string name;
                    for (int i = 0; i < nData; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        inp >> name;
                    }
                }
                {// Read the values
                    Rch.resize(nRasterCells,0.0);
                    for (int i = 0; i < nRasterCells; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        double v;
                        for (int j = 0; j < nData; ++j){
                            inp >> v;
                            if (j == 0){
                                Rch[i] = v;
                                break;
                            }
                        }
                    }
                }
                ifile.close();
            }
        }


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
