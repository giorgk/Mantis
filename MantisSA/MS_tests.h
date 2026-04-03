//
// Created by giorg on 6/9/2025.
//

#include <iostream>
#include <boost/mpi.hpp>

#include "MS_mpi_utils.h"

#ifndef MANTISSA_MS_TESTS_H
#define MANTISSA_MS_TESTS_H

namespace MS{
    namespace TESTS{
        void sendVec2Root(boost::mpi::communicator& world){
            std::vector<int> tmp;
            for (int i = 0; i < world.rank()*2+5; ++i){
                tmp.push_back(i + +world.rank()*10);
            }
            for (int irank = 0; irank < world.size(); ++irank){
                if (irank == world.rank()){
                    std::cout << "Rank " << world.rank() << " : ";
                    for (int i = 0; i < tmp.size(); ++i){
                        std::cout << " " << tmp[i];
                    }
                    std::cout << std::endl;
                }
            }
            std::vector<std::vector<int>> out;
            MS::sendVec2Root<int>(tmp, out,world);
            MS::printVectorForAllProc<int>(tmp, world);
            std::cout << world.rank() << " : " << out.size() << std::endl;
            if (world.rank() == 0){
                for (int i = 0; i < out.size(); ++i){
                    std::cout << i << " |";
                    for (int j = 0; j < out[i].size(); ++j){
                        std::cout << " " << out[i][j];
                    }
                    std::cout << "|" << std::endl;
                }
            }
        }

        void sendScalarFromRoot2AllProc(boost::mpi::communicator& world){
            std::string message = "___";
            int ii = 0;
            double v =0.0;

            if (world.rank() == 0) {
                message = "Hello from root!";
                ii = 111;
                v = 123.123;
            }
            std::cout << "BEFORE: Rank " << world.rank() << " received: " << message << " " << ii << " " << v << std::endl;
            MS::sendScalarFromRoot2AllProc(message, world);
            MS::sendScalarFromRoot2AllProc(ii, world);
            MS::sendScalarFromRoot2AllProc(v, world);
            std::cout << "AFTER: Rank " << world.rank() << " received: " << message << " " << ii << " " << v << std::endl;

        }

        void sendVectorFromRoot2AllProc(boost::mpi::communicator& world){
            std::vector<int> tmp_i;
            std::vector<double> tmp_d;
            for (int i = 0; i < 20; i = i+2){
                tmp_i.push_back(i);
                tmp_d.push_back(static_cast<double>(i)*2.3);
            }
            MS::sendVectorFromRoot2AllProc(tmp_i, world);
            MS::sendVectorFromRoot2AllProc(tmp_d, world);

            MS::printVectorForAllProc(tmp_i,world);
            MS::printVectorForAllProc(tmp_d,world);
        }

        void sentMatrixFromRoot2AllProc(boost::mpi::communicator& world){
            Matrix<int> m_i;
            Matrix<double> m_d;

            int nr;
            int nc;

            if (world.rank() == 0){
                nr = 10;
                nc = 5;

                m_i.allocate(nr, nc);
                m_d.allocate(nr, nc);

                for (int i = 0; i < nr; ++i) {
                    for (int j = 0; j < nc; ++j) {
                        m_i(i, j) = j + i * 10;
                        m_d(i, j) = static_cast<double>(j) + static_cast<double>(i) * 10.343;
                    }
                }
            }
            sendMatrixFromRoot2AllProc(m_i, world);
            sendMatrixFromRoot2AllProc(m_d, world);

            printMatrixForAllProc(m_i, world);
            printMatrixForAllProc(m_d, world);
        }
    }
}

#endif //MANTISSA_MS_TESTS_H
