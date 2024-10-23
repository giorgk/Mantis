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
        // Send the size to calculate the displacements
        if (world.rank() == 0){
            boost::mpi::gather(world, static_cast<int>(v.size()), vSize, 0);
            displ.push_back(0);
            for (int i = 0; i < vSize.size(); ++i){
                if (i!=0){
                    displ.push_back(displ[i-1] + vSize[i-1]);
                }
                totSize = totSize + vSize[i];
                //std::cout << "Size from rank:" << i << " " <<  vSize[i] << " displ " <<  displ[i] << std::endl;
            }
        }
        else{
            boost::mpi::gather(world, static_cast<int>(v.size()),0);
        }
        world.barrier();

        //Send the data
        std::vector<T> allValues;
        if (world.rank() == 0){
            allValues.resize(totSize,-1);
            boost::mpi::gatherv(world, v, &allValues[0], vSize, displ,0);
        }
        else{
            boost::mpi::gatherv(world, v, 0);
        }
        world.barrier();

        // Root processor gathers the data
        out.clear();
        out.resize(world.size());
        displ.push_back(totSize);
        int idx = 0;
        if (world.rank() == 0){
            for (int i = 0; i < world.size(); ++i){
                for (int j = displ[i]; j < displ[i+1]; ++j){
                    out[i].push_back(allValues[idx]);
                    idx = idx + 1;;
                }
            }
            //for (int i = 0; i < world.size(); ++i){
            //    for (int j = 0; j < out[i].size(); ++j){
            //        std::cout << out[i][j] << " ";
            //    }
            //    std::cout << std::endl;
            //}
        }
        world.barrier();
    }

    //void send

    void sendWellEids(WELLS & wells, std::vector<std::vector<int>> & wellids, boost::mpi::communicator & world){
        //std::cout << "Rank " << world.rank() << " nWells: " << wells.size() << std::endl;
        std::vector<int> proc_eids;
        WELLS::iterator itw;
        int cnt = 0;
        for (itw = wells.begin(); itw != wells.end(); ++itw){
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
    void sentFromRootToAllProc(std::vector<T> &v, boost::mpi::communicator &world){
        boost::mpi::broadcast(world, &v[0], static_cast<int>(v.size()), 0);
    }
}

#endif //MANTISSA_MS_MPI_UTILS_H
