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
    bool readMatrix(std::string filename, std::vector<std::vector<T>> & data, int nCols){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            T v;
            while (getline(datafile, line)){
                std::istringstream inp(line.c_str());
                std::vector<T> d;
                for (int i = 0; i < nCols; ++i){
                    inp >> v;
                    d.push_back(v);
                }
                data.push_back(d);
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


}

#endif //MANTISSA_MS_STRUCTURES_H
