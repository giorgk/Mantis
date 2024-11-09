//
// Created by giorg on 10/29/2024.
//

#ifndef MANTISSA_MSDEBUG_H
#define MANTISSA_MSDEBUG_H

namespace MS{
    void testConvolution(){
        std::vector<double> urf, lf, btc, prebtc;
        MS::calcURFs(200, 4.782, 0.4213, 140.32, 1074.7, urf);
        double x = 1;
        for (int i = 0; i < 100; ++i){
            lf.push_back(x);
            x = x+1.0;
        }
        x = 100;
        for (int i = 0; i < 100; ++i){
            lf.push_back(x);
            x = x - 1;
        }
        btc.resize(lf.size(),0);
        prebtc.resize(lf.size(),0);
        MS::convolute(urf,lf,btc,prebtc,100.0);
        for (int i = 0; i < urf.size(); ++i){
            std::cout << urf[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < lf.size(); ++i){
            std::cout << lf[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < btc.size(); ++i){
            std::cout << btc[i] << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < prebtc.size(); ++i){
            std::cout << prebtc[i] << " ";
        }
        std::cout << std::endl;
    }

    void testBroadcast(boost::mpi::communicator& world){
        int N=30000000;
        std::vector<double> tmp(N,0);
        if (world.rank() == 0) {
            for (int i = 0; i < tmp.size(); ++i)
                tmp[i] = static_cast<double>(i)*2.36547;
        }
        //boost::mpi::broadcast(world, &tmp[0],N, 0);
        MS::sendVectorFromRoot2AllProc<double>(tmp, world);

        for (int k = 0; k < world.size(); ++k){
            for (int i = N-20; i < N; ++i){
                std::cout << tmp[i] << " ";
            }
            std::cout << std::endl;
        }
    }

    void printConcFromPump(std::string &fname, std::vector<double> &conc){
        std::ofstream out_file;
        out_file.open(fname.c_str());
        for (unsigned int i = 0; i < conc.size(); ++i){
            out_file << conc[i] << std::endl;
        }
        out_file.close();

    }

    template<typename T>
    void testMatrixBroadcast(boost::mpi::communicator& world){
        // Create the matrix
        std::vector<std::vector<T>> vMat;
        if (world.rank() == 0){
            int idx = 1;
            for (int i = 0; i < 7; ++i){
                std::vector<T> vVec;
                for (int j = 0; j < 10; ++j ){
                    vVec.push_back(static_cast<T>(idx));
                    idx = idx + 1;
                }
                vMat.push_back(vVec);
            }
        }
        else{
            std::vector<T> vVec;
            vVec.resize(10,0);
            vMat.resize(7,vVec);
        }
        world.barrier();

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                for (unsigned int i = 0; i < vMat.size(); ++i){
                    for (unsigned int j = 0; j < vMat[i].size(); ++j){
                        std::cout << vMat[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            world.barrier();
        }


        sentMatrixFromRoot2AllProc(vMat, world);

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                for (unsigned int i = 0; i < vMat.size(); ++i){
                    for (unsigned int j = 0; j < vMat[i].size(); ++j){
                        std::cout << vMat[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            world.barrier();
        }
    }

    void testScalarBroadcast(boost::mpi::communicator& world){
        int ia = 0;
        double ib= 0.0;
        if (world.rank() == 0){
            ia = 100;
            ib = 99.9999;
        }
        world.barrier();
        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " " << ia << " " << ib << std::endl;
            }
            world.barrier();
        }
        MS::sendScalarFromRoot2AllProc<int>(ia, world);
        MS::sendScalarFromRoot2AllProc<double>(ib, world);

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " " << ia << " " << ib << std::endl;
            }
            world.barrier();
        }

    }


}
#endif //MANTISSA_MSDEBUG_H
