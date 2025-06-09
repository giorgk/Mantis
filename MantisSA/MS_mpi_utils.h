//
// Created by giorg on 10/18/2024.
//

#ifndef MANTISSA_MS_MPI_UTILS_H
#define MANTISSA_MS_MPI_UTILS_H
namespace MS{
    template<typename T>
    void sendVec2Root(std::vector<T> &v, std::vector<std::vector<T>> &out, boost::mpi::communicator &world){
        std::vector<int> vSize;
        std::vector<int> displ;
        int totSize = 0;

        // Step 1: Gather sizes at root
        // Send the size to calculate the displacements
        int localSize = static_cast<int>(v.size());
        if (world.rank() == 0){
            boost::mpi::gather(world, localSize, vSize, 0);
            assert(static_cast<int>(vSize.size()) == world.size() && "Size mismatch in gathered sizes");

            // Step 2: Compute displacements
            displ.reserve(world.size());
            displ.push_back(0);
            for (int i = 1; i < vSize.size(); ++i){
                //if (i!=0){
                    displ.push_back(displ[i-1] + vSize[i-1]);
                //}
                //totSize = totSize + vSize[i];
                //std::cout << "Size from rank:" << i << " " <<  vSize[i] << " displ " <<  displ[i] << std::endl;
            }
            totSize = displ.back() + vSize.back();
        }
        else{
            boost::mpi::gather(world, static_cast<int>(v.size()),0);
        }
        world.barrier();

        // Step 3: Gather data from all processes
        //Send the data
        std::vector<T> allValues;
        if (world.rank() == 0){
            allValues.resize(totSize);
            boost::mpi::gatherv(world, v, allValues.data(), vSize, displ,0);
        }
        else{
            boost::mpi::gatherv(world, v, 0);
        }
        world.barrier();

        // Step 4: Reconstruct vectors from gathered data
        // Root processor gathers the data
        if (world.rank() == 0){
            out.clear();
            out.resize(world.size());
            for (int i = 0; i < world.size(); ++i){
                int start = displ[i];
                int end = start + vSize[i];
                out[i].assign(allValues.begin() + start, allValues.begin() + end);
            }

        }
        //displ.push_back(totSize);
        //int idx = 0;
        //if (world.rank() == 0){
        //    for (int i = 0; i < world.size(); ++i){
        //        for (int j = displ[i]; j < displ[i+1]; ++j){
        //            out[i].push_back(allValues[idx]);
        //            idx = idx + 1;
        //        }
        //    }
            //for (int i = 0; i < world.size(); ++i){
            //    for (int j = 0; j < out[i].size(); ++j){
            //        std::cout << out[i][j] << " ";
            //    }
            //    std::cout << std::endl;
            //}
        //}
        world.barrier();
    }


    void sendWellEids(WELLS & wells, std::vector<std::vector<int>> & wellids, boost::mpi::communicator & world){
        //std::cout << "Rank " << world.rank() << " nWells: " << wells.size() << std::endl;
        std::vector<int> proc_eids;
        std::map<int, WELL>::iterator itw;
        int cnt = 0;
        for (itw = wells.wells.begin(); itw != wells.wells.end(); ++itw){
            proc_eids.push_back(itw->first);
            //if (cnt < 20){
            //    std::cout << world.rank() << ": " << itw->first <<std::endl;
            //    cnt++;
            //}
        }

        sendVec2Root<int>(proc_eids, wellids, world);


        //if (world.rank() == 0){
            //for (int i = 0; i < world.size(); ++i){
            //    std::cout << "Rank " << i << " wellids: " << wellids[i].size() << std::endl;
            //}

            //for (int i = 0; i < wellids.size(); ++i){
            //    for (int j = 0; j < 50 /*wellids[i].size()*/; ++j){
            //        std::cout << wellids[i][j] << " ";
            //    }
            //    std::cout << std::endl;
            //}
        //}
    }

    template<typename T>
    void sendScalarFromRoot2AllProc(T &a, boost::mpi::communicator &world){
        //boost::mpi::broadcast(world, &a, 1, 0);
        boost::mpi::broadcast(world, a,0);
    }

    template<typename T>
    void sendVectorFromRoot2AllProc(std::vector<T>& v, boost::mpi::communicator& world){
        // boost::mpi::broadcast(world, v, 0); we cannot use this with templates
        int size = static_cast<int>(v.size());
        boost::mpi::broadcast(world, size, 0);  // Broadcast the size of the vector

        if (world.rank() != 0) {
            v.resize(size);  // Ensure the vector has the right size before receiving data
        }

        if (size > 0){
            boost::mpi::broadcast(world, v.data(), size, 0);
        }
        //boost::mpi::broadcast(world, &v[0], size, 0);
    }

    template<typename T>
    void sentMatrixFromRoot2AllProc(std::vector<std::vector<T>> &v, boost::mpi::communicator &world){
        int nRows = 0;
        int nCols = 0;
        if (world.rank() == 0){
            nRows = static_cast<int>(v.size());
            nCols = nRows > 0 ? static_cast<int>(v[0].size()) : 0;
        }
        world.barrier();
        // Send the dimensions and initialize the matrix for the other processors
        sendScalarFromRoot2AllProc(nRows,world);
        sendScalarFromRoot2AllProc(nCols, world);
        if (world.rank() > 0){
            v.assign(nRows, std::vector<T>(nCols, T{}));
        }
        world.barrier();
        for (int i = 0; i < nRows; ++i){
            sendVectorFromRoot2AllProc<T>(v[i], world);
        }
    }

    template<typename T>
    void sendFlatMatrixFromRoot2AllProc(std::vector<std::vector<T>> &matrix, boost::mpi::communicator &world){
        int nRows = 0;
        int nCols = 0;

        if (world.rank() == 0) {
            nRows = static_cast<int>(matrix.size());
            nCols = nRows > 0 ? static_cast<int>(matrix[0].size()) : 0;
        }

        sendScalarFromRoot2AllProc(nRows, world);
        sendScalarFromRoot2AllProc(nCols, world);

        std::vector<T> flatMatrix;
        if (world.rank() == 0){
            flatMatrix.reserve(nRows * nCols);
            for (const auto& row : matrix){
                flatMatrix.insert(flatMatrix.end(), row.begin(), row.end());
            }
        }
        sendVectorFromRoot2AllProc(flatMatrix, world);
        if (world.rank() != 0){
            matrix.assign(nRows, std::vector<T>(nCols));
            for (int i = 0; i < nRows; ++i){
                std::copy(flatMatrix.begin() + i * nCols, flatMatrix.begin() + (i + 1) * nCols, matrix[i].begin());
            }
        }
    }

    template<typename T>
    bool RootReadsMatrixFileDistrib(std::string filename, std::vector<std::vector<T>> &M,
                                    int nCols, boost::mpi::communicator &world, int freq = 500000){
        int readSuccess = 0;
        if (world.rank() == 0){
            bool tf = readMatrix<T>(filename, M,  nCols, freq);
            readSuccess = static_cast<int>(tf);
        }
        world.barrier(); //Wait for the root processor to read_v0

        // Check if reading was successful
        sendScalarFromRoot2AllProc<int>(readSuccess, world);
        if (readSuccess == 0){
            return false;
        }
        else{
            sendFlatMatrixFromRoot2AllProc(M, world);
            return true;
        }

    }


    template<typename T>
    bool RootReadsMatrixFileDistrib(std::string filename, std::vector<std::vector<T>> &M,
                                    int nCols, bool doTranspose, boost::mpi::communicator &world, int freq = 500000){
        int readSuccess = 0;
        std::vector<std::vector<T>> Mtmp;
        std::vector<std::vector<T>> Mtmp1;
        if (world.rank() == 0){
            bool tf = readMatrix<T>(filename, Mtmp,  nCols, freq);
            readSuccess = static_cast<int>(tf);
        }
        world.barrier(); //Wait for the root processor to read_v0


        // Check if reading was successful
        sendScalarFromRoot2AllProc<int>(readSuccess, world);
        if (readSuccess == 0){
            return false;
        }
        else{
            int preTranspose = 0;
            if (world.rank() == 0){
                // We check which dimension is larger
                if (Mtmp.size() > Mtmp[0].size()){
                    transposeMatrix<T>(Mtmp, Mtmp1);
                    preTranspose = 1;
                }
                else {
                    Mtmp1 = Mtmp;
                }
            }
            // Send if the matrix has been transposed before sending
            sendScalarFromRoot2AllProc<int>(preTranspose, world);
            sentMatrixFromRoot2AllProc<T>(Mtmp1, world);
            if (preTranspose == 1){
                if (!doTranspose){
                    // We have to transpose back to its original size
                    transposeMatrix<T>(Mtmp1, M);
                }
                else {
                    M = Mtmp1;
                }
            }
            else{
                if (doTranspose){
                    transposeMatrix<T>(Mtmp1, M);
                }
                else{
                    M = Mtmp1;
                }
            }
        }
        return true;
    }

    template<typename T>
    void printMatrixForAllProc(std::vector<std::vector<T>> & M, boost::mpi::communicator& world,
                               int startRow = -9, int endRow = -9,
                               int startCol = -9, int endCol = -9){
        if (startRow < 0){startRow = 0;}
        if (endRow < 0){endRow = M.size();}
        if (startCol < 0){startCol = 0;}
        if (endCol < 0){endCol = M[0].size();}

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                for (unsigned int i = startRow; i < endRow; ++i){
                    for (unsigned int j = startCol; j < endCol; ++j){
                        std::cout << M[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            world.barrier();
        }
    }

    template<typename T>
    void printVectorForAllProc(std::vector<T> &V, boost::mpi::communicator& world,
                               int startRow = -9, int endRow = -9){
        if (startRow < 0){startRow = 0;}
        if (endRow < 0){endRow = V.size();}

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " [";
                for (unsigned int i = startRow; i < endRow; ++i){
                    std::cout << V[i] << " ";
                }
                std::cout << "]" << std::endl;
            }
            world.barrier();
        }
    }

    bool RootReadsNPSATinfo(std::string filename, std::vector<int> &V, boost::mpi::communicator &world){
        int readSuccess = 0;
        V.clear();
        world.barrier();
        if (world.rank() == 0){
            bool tf = readNPSATinfo(filename, V);
            readSuccess = static_cast<int>(tf);
        }
        world.barrier(); //Wait for the root processor to read_v0

        // Check if reading was successful
        sendScalarFromRoot2AllProc<int>(readSuccess, world);
        if (readSuccess == 0){
            return false;
        }
        else{
            sendVectorFromRoot2AllProc<int>(V,world);
            return true;
        }
    }


}

#endif //MANTISSA_MS_MPI_UTILS_H
