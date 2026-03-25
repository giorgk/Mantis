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
            if (matrix.empty()) {
                throw std::runtime_error("sendFlatMatrixFromRoot2AllProc: matrix is empty");
            }

            if (matrix.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("sendFlatMatrixFromRoot2AllProc: too many rows");
            }

            nRows = static_cast<int>(matrix.size());

            if (matrix[0].size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("sendFlatMatrixFromRoot2AllProc: too many columns");
            }

            nCols = static_cast<int>(matrix[0].size());

            if (nCols <= 0) {
                throw std::runtime_error("sendFlatMatrixFromRoot2AllProc: matrix has zero columns");
            }

            for (int i = 1; i < nRows; ++i) {
                if (static_cast<int>(matrix[i].size()) != nCols) {
                    throw std::runtime_error("sendFlatMatrixFromRoot2AllProc: matrix is not rectangular");
                }
            }
        }

        sendScalarFromRoot2AllProc(nRows, world);
        sendScalarFromRoot2AllProc(nCols, world);

        std::vector<T> flatMatrix;

        if (world.rank() == 0){
            const std::size_t totalSize = static_cast<std::size_t>(nRows) * static_cast<std::size_t>(nCols);

            flatMatrix.resize(totalSize);

            std::size_t k = 0;
            for (int i = 0; i < nRows; ++i) {
                for (int j = 0; j < nCols; ++j) {
                    flatMatrix[k++] = matrix[i][j];
                }
            }
        }

        sendVectorFromRoot2AllProc(flatMatrix, world);

        if (world.rank() != 0){
            matrix.assign(nRows, std::vector<T>(nCols));

            for (int i = 0; i < nRows; ++i) {
                std::copy(flatMatrix.begin() + static_cast<std::size_t>(i) * nCols,
                          flatMatrix.begin() + static_cast<std::size_t>(i + 1) * nCols,
                          matrix[i].begin());
            }
        }
    }

    template<typename T>
    bool RootReadsMatrixFileDistribFlat(const std::string& filename,
                                    std::vector<std::vector<T>>& M,
                                    int nCols,
                                    boost::mpi::communicator& world,
                                    int freq = 500000)
    {
        int readSuccess = 0;

        if (world.rank() == 0) {
            try {
                const bool tf = readMatrix<T>(filename, M, nCols, freq);
                readSuccess = static_cast<int>(tf);

                if (readSuccess) {
                    if (M.empty()) {
                        std::cout << "Error: matrix file " << filename
                                  << " contains no valid data rows" << std::endl;
                        readSuccess = 0;
                    }
                    else if (M[0].empty()) {
                        std::cout << "Error: matrix file " << filename
                                  << " has zero columns" << std::endl;
                        readSuccess = 0;
                    }
                    else {
                        const std::size_t nCols0 = M[0].size();
                        for (std::size_t i = 1; i < M.size(); ++i) {
                            if (M[i].size() != nCols0) {
                                std::cout << "Error: matrix file " << filename
                                          << " is not rectangular" << std::endl;
                                readSuccess = 0;
                                break;
                            }
                        }
                    }
                }
            }
            catch (const std::exception& e) {
                std::cout << "Error in RootReadsMatrixFileDistrib while reading "
                          << filename << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
            catch (...) {
                std::cout << "Unknown error in RootReadsMatrixFileDistrib while reading "
                          << filename << std::endl;
                readSuccess = 0;
            }
        }

        sendScalarFromRoot2AllProc(readSuccess, world);

        if (readSuccess == 0) {
            M.clear();
            return false;
        }

        try {
            sendFlatMatrixFromRoot2AllProc(M, world);
        }
        catch (const std::exception& e) {
            std::cout << "Error in RootReadsMatrixFileDistrib while distributing "
                      << filename << ": " << e.what() << std::endl;
            M.clear();
            return false;
        }
        catch (...) {
            std::cout << "Unknown error in RootReadsMatrixFileDistrib while distributing "
                      << filename << std::endl;
            M.clear();
            return false;
        }

        return true;
    }

    template<typename T>
    bool RootReadsMatrixFileDistrib(const std::string& filename, std::vector<std::vector<T>>& M,
                                    int nCols, bool doTranspose, boost::mpi::communicator& world, int freq = 500000){
        int readSuccess = 0;
        int preTranspose = 0;

        std::vector<std::vector<T>> Mtmp;
        std::vector<std::vector<T>> Mtmp1;

        if (world.rank() == 0){
            try
            {
                const bool tf = readMatrix<T>(filename, Mtmp, nCols, freq);
                readSuccess = static_cast<int>(tf);
                if (readSuccess) {
                    if (Mtmp.empty()) {
                        std::cout << "Error: matrix is empty after reading " << filename << std::endl;
                        readSuccess = 0;
                    }
                    else {
                        if (Mtmp.size() > Mtmp[0].size()) {
                            transposeMatrix<T>(Mtmp, Mtmp1);
                            preTranspose = 1;
                        }
                        else {
                            Mtmp1 = Mtmp;
                        }
                    }
                }
            }
            catch (const std::exception& e) {
                std::cout << "Error in RootReadsMatrixFileDistrib while reading "
                          << filename << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
            catch (...) {
                std::cout << "Unknown error in RootReadsMatrixFileDistrib while reading "
                          << filename << std::endl;
                readSuccess = 0;
            }
        }

        world.barrier(); //Wait for the root processor to read_v0


        // Check if reading was successful
        sendScalarFromRoot2AllProc<int>(readSuccess, world);
        if (readSuccess == 0){
            M.clear();
            return false;
        }

        sendScalarFromRoot2AllProc<int>(preTranspose, world);
        sentMatrixFromRoot2AllProc<T>(Mtmp1, world);


        try {
            if (preTranspose == 1) {
                if (!doTranspose) {
                    transposeMatrix<T>(Mtmp1, M);
                } else {
                    M = Mtmp1;
                }
            }
            else {
                if (doTranspose) {
                    transposeMatrix<T>(Mtmp1, M);
                } else {
                    M = Mtmp1;
                }
            }
        }
        catch (const std::exception& e) {
            std::cout << "Error in RootReadsMatrixFileDistrib after distribution: "
                      << e.what() << std::endl;
            return false;
        }
        catch (...) {
            std::cout << "Unknown error in RootReadsMatrixFileDistrib after distribution"
                      << std::endl;
            return false;
        }
        return true;
    }

    template<typename T>
    void printMatrixForAllProc(std::vector<std::vector<T>> & M, boost::mpi::communicator& world,
                               int startRow = -9, int endRow = -9,
                               int startCol = -9, int endCol = -9){
        const std::size_t nRows = M.size();
        const std::size_t nCols = M[0].size();

        std::size_t i0 = (startRow < 0) ? 0 : static_cast<std::size_t>(startRow);
        std::size_t i1 = (endRow   < 0) ? nRows : static_cast<std::size_t>(endRow);
        std::size_t j0 = (startCol < 0) ? 0 : static_cast<std::size_t>(startCol);
        std::size_t j1 = (endCol   < 0) ? nCols : static_cast<std::size_t>(endCol);

        if (i0 > nRows) i0 = nRows;
        if (i1 > nRows) i1 = nRows;
        if (j0 > nCols) j0 = nCols;
        if (j1 > nCols) j1 = nCols;


        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << "  ----------" << std::endl;
                 for (std::size_t i = i0; i < i1; ++i){
                    for (std::size_t j = j0; j < j1; ++j){
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
        const std::size_t nRows = V.size();

        std::size_t i0 = (startRow < 0) ? 0 : static_cast<std::size_t>(startRow);
        std::size_t i1 = (endRow   < 0) ? nRows : static_cast<std::size_t>(endRow);

        if (i0 > nRows) i0 = nRows;
        if (i1 > nRows) i1 = nRows;

        for (int irank = 0; irank < world.size(); ++irank){
            if (irank == world.rank()){
                std::cout << "Rank " << world.rank() << " [";
                for (std::size_t i = i0; i < i1; ++i){
                    std::cout << V[i] << " ";
                }
                std::cout << "]" << std::endl;
            }
            world.barrier();
        }
    }

    inline bool RootReadsNPSATinfo(const std::string& filename, std::vector<int>& V, boost::mpi::communicator& world){
        int readSuccess = 0;
        V.clear();

        world.barrier();
        if (world.rank() == 0){
            try {
                const bool tf = readNPSATinfo(filename, V);
                readSuccess = static_cast<int>(tf);
            }
            catch (const std::exception& e) {
                std::cout << "Error reading " << filename << ": "
                          << e.what() << std::endl;
                readSuccess = 0;
            }
            catch (...) {
                std::cout << "Unknown error reading " << filename << std::endl;
                readSuccess = 0;
            }
        }
        world.barrier(); //Wait for the root processor to read_v0

        // Check if reading was successful
        sendScalarFromRoot2AllProc<int>(readSuccess, world);

        if (readSuccess == 0){
            V.clear();
            return false;
        }

        sendVectorFromRoot2AllProc<int>(V,world);
        return true;
    }


}

#endif //MANTISSA_MS_MPI_UTILS_H
