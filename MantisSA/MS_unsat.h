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
                      double wc_in, double d_in, double r_in, boost::mpi::communicator& world);

        int traveltime(int IJ);
        double getDepth(int IJ);
        double getRch(int IJ);
        double getSurfPerc(int IJ);
    private:
        double wc;
        double minD;
        double minR;
        std::vector<double> Depth;
        std::vector<double> Rch;
        std::vector<double> SurfPerc;
        int Ncells;
    };

    bool UNSAT::readdata(std::string dpth_file, std::string dpth_name, std::string rch_file,
                         double wc_in, double d_in, double r_in,
                         boost::mpi::communicator& world) {
        wc = wc_in;
        minD = d_in;
        minR = r_in;
        std::string ext = getExtension(dpth_file);
        if (ext.compare("h5") == 0) { // Read the Depth file
#if _USEHF>0
            if (world.rank() == 0) {
                std::cout << "Reading HDF5 file " << dpth_file << std::endl;
            }
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

            int readSuccess = 0;
            if (world.rank() == 0){
                int depth_idx = 0;
                std::vector<std::vector<double>> data;
                std::vector<std::string> names;
                bool tf = readLinearData(dpth_file, names, data, 5000000);
                if (tf) {
                    readSuccess = 1;
                    for (unsigned int i = 0; i < names.size(); ++i){
                        if (names[i] == dpth_name){
                            depth_idx = static_cast<int>(i);
                            break;
                        }
                    }
                    Ncells = static_cast<int>(data[0].size());
                    Depth.resize(Ncells, 0);
                    for (unsigned int i = 0; i < data[0].size(); ++i){
                        Depth[i] = data[depth_idx][i];
                    }
                }
            }
            world.barrier();

            sendScalarFromRoot2AllProc<int>(readSuccess, world);
            if (readSuccess == 0){
                return false;
            }
            else{
                // We send the number of rows for the data initialization
                sendScalarFromRoot2AllProc<int>(Ncells, world);
                if (world.rank() > 0){
                    // All other processor initialize the data matrix
                    Depth.resize(Ncells, 0);
                }
                world.barrier();
                sendVectorFromRoot2AllProc<double>(Depth, world);
            }
        }
        if (PrintMatrices){
            printVectorForAllProc<double>(Depth,world, 0, 30);
        }

        // Read Recharge
        std::string ext1 = getExtension(rch_file);
        if (ext1.compare("h5") == 0)
        {   // Read Recharge
#if _USEHF>0
            if (world.rank() == 0) {
                std::cout << "Reading HDF5 file " << rch_file << std::endl;
            }
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
            SurfPerc = tmp[1];
            if (Ncells != Rch.size()){
                std::cout << "The size of Depth is different than the size of Recharge" << std::endl;
                return false;
            }
#endif
        }
        else{
            int readSuccess = 0;
            if (world.rank() == 0){
                std::vector<std::vector<double>> data;
                std::vector<std::string> names;
                bool tf = readLinearData(rch_file, names, data, 5000000);
                if (tf){
                    if (data[0].size() == Ncells){
                        readSuccess = 1;
                        Rch.resize(Ncells, 0);
                        SurfPerc.resize(Ncells, 0);
                        for (unsigned int i = 0; i < data[0].size(); ++i){
                            Rch[i] = data[0][i];
                            SurfPerc[i] = data[1][i];
                        }
                    }
                }
            }
            world.barrier();
            sendScalarFromRoot2AllProc<int>(readSuccess, world);
            if (readSuccess == 0){
                return false;
            }
            else{
                if (world.rank() > 0){
                    Rch.resize(Ncells, 0);
                    SurfPerc.resize(Ncells, 0);
                }
                world.barrier();
                sendVectorFromRoot2AllProc<double>(Rch, world);
                sendVectorFromRoot2AllProc<double>(SurfPerc, world);
            }

            if(PrintMatrices){
                printVectorForAllProc<double>(Rch,world, 3483753, 3483753+20);
                printVectorForAllProc<double>(SurfPerc,world, 3483753, 3483753+20);
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

    double UNSAT::getDepth(int IJ) {
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return Depth[IJ];
    }

    double UNSAT::getRch(int IJ) {
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return Rch[IJ];
    }
    double UNSAT::getSurfPerc(int IJ){
        if (IJ < 0 || IJ >= Ncells ){
            return 0.0;
        }
        return SurfPerc[IJ];
    }

}

#endif //MANTISSA_MS_UNSAT_H
