//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_STRUCTURES_H
#define MANTISSA_MS_STRUCTURES_H
#include <string>
#include <fstream>
namespace MS{

    struct SelectedWellsGroup{
        std::string groupName;
        std::vector<int> idVI;
        std::vector<int> idVD;
    };


    struct RasterOptions {
        int Ncells;
        int Nrows;
        int Ncols;
        double Xorig;
        double Yorig;
        double CellSize;
        std::string File;
    };

    struct SwatOptions{
        int Nyears;
        int StartYear;
        int version;
        std::string HRU_raster_file;
        std::string HRU_index_file;
        std::string Data_file;
    };

    struct HistoricOptions{
        int StartYear;
        int EndYear;
        int Interval;
        int BlendStart;
        int BlendEnd;
        std::string filename;
        std::string ext;

    };

    struct NpsatOptions{
        bool bUseInitConcVI;
        bool bUseInitConcVD;
        int version;
        std::string VIdataFile;
        std::string VDdataFile;
        std::string InitSaltVIFile;
        std::string InitSaltVDFile;
        std::string DistribuPumpFile;
    };

    struct RiverOptions{
        double StartDist;
        double EndDist;
        double ConcValue;
    };

    struct SimOptions{
        bool EnableFeedback;
        bool OutofAreaUseInitConc;
        int StartYear;
        int Nyears;
        int Porosity;
        int nBuffer;
        double MaxConc;
        double MaxAge;
        double SurfConcValue;
        double OutofAreaConc;
    };

    struct UnsatOptions{
        double wc;
        double minDepth;
        double minRch;
        std::string Depth_file;
        std::string Depth_name;
        std::string Rch_file;
    };

    struct OutputOptions{
        bool printLoad;
        bool printURFs;
        bool printBTCs;
        bool printSelectedWells;
        std::string OutFile;
        std::string SelectedWells;
        std::string SelectedWellGroups;
    };

    struct STRML{
        int Sid;
        int urfI;
        int urfJ;
        int IJ;
        bool inRiv;
        double rivRist;
        //int hru_idx;
        double W;
        double Len;
        double m;
        double s;
        double a;
        int wellSourceId;
        std::vector<double> urf;
        std::vector<double> lf_conc;
        //std::vector<double> lf_mass;
        std::vector<double> gw_conc;
        //std::vector<double> gw_conc;
        std::vector<double> btc;

    };

    struct WELL{
        double initConc = 0.0;
        std::vector<STRML> strml;
        std::vector<double> wellBtc;
    };

    struct NPSATTMP{
        double sumW = 0;
        std::vector<int> Sid;
        std::vector<int> urfI;
        std::vector<int> urfJ;
        std::vector<bool> inRiv;
        std::vector<double> rivDist;
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

    std::string getExtension(std::string filename){
        // https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
        // Solution 1.3: stepping away from the iterators
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(filename);
        while (std::getline(tokenStream, token, '.')){
            tokens.push_back(token);
        }
        return tokens.back();
    }

    template<typename T>
    bool readVector(std::string filename, std::vector<T> & data){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std:: cout << "Reading " << filename << std::endl;
            std::string line;
            T v;
            while (getline(datafile, line)){
                std::istringstream inp(line.c_str());
                inp >> v;
                data.push_back(v);
            }
            datafile.close();
        }
        return true;
    }

    bool readSelectedWellsGroupInfo(std::string filename, std::map<int,SelectedWellsGroup>& swg_map){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            swg_map.clear();
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            while (getline(datafile, line)){
                if (line.size() > 0){
                    std::istringstream inp(line.c_str());
                    std::string groupName;
                    int groupId;
                    inp >> groupId;
                    inp >> groupName;
                    std::map<int,SelectedWellsGroup>::iterator it;
                    it = swg_map.find(groupId);
                    if (it != swg_map.end()){
                        std::cout << "Selected Wells Group with id " << groupId << " was found more than once" <<std::endl;
                        return false;
                    }
                    else{
                        SelectedWellsGroup swg;
                        swg.groupName = groupName;
                        swg_map.insert(std::pair<int, SelectedWellsGroup>(groupId, swg));
                    }
                }
            }
            datafile.close();
        }
        return true;
    }

    template<typename T>
    bool readMatrix(std::string filename, std::vector<std::vector<T>> & data, int nCols, int freq = 500000){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            T v;
            int countLines = freq;
            int ir = 0;
            //for (int ir = 0; ir < nRows; ++ir){
            while (getline(datafile, line)){
                if (line.size() > 0){
                    if (ir > countLines){
                        std::cout << ir << std::endl;
                        countLines = countLines + freq;
                    }
                    //std::cout << line << std::endl;
                    std::istringstream inp(line.c_str());
                    std::vector<T> d;
                    for (int i = 0; i < nCols; ++i){
                        inp >> v;
                        d.push_back(v);
                    }
                    data.push_back(d);
                    ir = ir + 1;
                }
            }
            datafile.close();
        }
        return true;
    }

    template<typename T>
    void transposeMatrix(std::vector<std::vector<T>> & data_in, std::vector<std::vector<T>> & data_out){
        data_out.clear();
        data_out.resize(data_in[0].size(),std::vector<T>(data_in.size(),0));
        for (unsigned int i = 0; i < data_in.size(); ++i){
            for (unsigned int j = 0; j < data_in[i].size(); ++j){
                data_out[j][i] = data_in[i][j];
            }
        }
    }


    bool readLinearData(std::string filename, std::vector<std::string> &names,
                        std::vector<std::vector<double>> &data, int freq = 500000){
        std::ifstream ifile;
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            int nData;
            std::string line;
            { //Get the number of Unsaturated scenarios
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> nData;
                data.clear();
                data.resize(nData,std::vector<double>());
            }

            {// Get the scenario names
                std::string name;
                for (int i = 0; i < nData; ++i){
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> name;
                    names.push_back(name);
                }
            }
            {// Read the values
                int ir = 0;
                int countLines = freq;
                while (getline(ifile, line)) {
                    if (line.size() > 0){
                        if (ir > countLines){
                            std::cout << ir << std::endl;
                            countLines = countLines + freq;
                        }

                        std::istringstream inp(line.c_str());
                        double v;
                        for (int j = 0; j < nData; ++j){
                            inp >> v;
                            data[j].push_back(v);
                        }
                        ir = ir + 1;
                    }
                }
            }
            ifile.close();
        }
        return true;
    }

    double calcRiverInfluence(const double &dist, const RiverOptions &ropt){
        double u = 1.0;
        if (dist < ropt.StartDist){
            u = 0.0;
        }
        else if (dist > ropt.EndDist){
            u = 1.0;
        }
        else{
            u = (dist - ropt.StartDist)/(ropt.EndDist - ropt.EndDist);
        }
        return u;
    }

    bool readNPSATinfo(std::string filename, std::vector<int> &M){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()){
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line, propname;
            int count = 0;
            int version;
            std::vector<int> por_values, int_size, dbl_size, msa_size;
            while (getline(datafile, line)){
                std::istringstream inp(line);
                inp >> propname;
                if (propname.compare("Version") == 0){
                    inp >> version;
                    continue;
                }
                if (propname.compare("Por") == 0){
                    int value;
                    while (inp >> value) {
                        por_values.push_back(value);
                    }
                    continue;
                }
                if (propname.compare("INT") == 0){
                    int value;
                    while (inp >> value) {
                        int_size.push_back(value);
                    }
                    continue;
                }
                if (propname.compare("DBL") == 0){
                    int value;
                    while (inp >> value) {
                        dbl_size.push_back(value);
                    }
                    continue;
                }
                if (propname.compare("MSA") == 0){
                    int value;
                    while (inp >> value) {
                        msa_size.push_back(value);
                    }
                    continue;
                }

                count++;
                if (count > 40){
                    return false;
                }
            }
            datafile.close();
            M.push_back(version);
            M.push_back(por_values.size());
            M.insert(M.end(), por_values.begin(), por_values.end());
            M.insert(M.end(), int_size.begin(), int_size.end());
            M.insert(M.end(), dbl_size.begin(), dbl_size.end());
            M.insert(M.end(), msa_size.begin(), msa_size.end());

            return true;
        }
    }

}

#endif //MANTISSA_MS_STRUCTURES_H
