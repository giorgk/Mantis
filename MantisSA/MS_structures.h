//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_STRUCTURES_H
#define MANTISSA_MS_STRUCTURES_H
#include <string>
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
}

#endif //MANTISSA_MS_STRUCTURES_H
