//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_STRUCTURES_H
#define MANTISSA_MS_STRUCTURES_H
#include <string>
#include <fstream>
namespace MS{




    struct RasterOptions {
        int Ncells;
        int Nrows;
        int Ncols;
        double Xorig;
        double Yorig;
        double CellSize;
        std::string File;
    };

    struct STRML{
        int Sid;
        int urfI;
        int urfJ;
        int IJ;
        //int hru_idx;
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
        std::vector<double> btc;
    };

    struct NPSATTMP{
        double sumW = 0;
        std::vector<int> Sid;
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
                if (line.size() > 1){
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
                    if (line.size() > 1){
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

}

#endif //MANTISSA_MS_STRUCTURES_H
