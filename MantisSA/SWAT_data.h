//
// Created by giorg on 10/16/2024.
//

#ifndef MANTISSA_SWAT_DATA_H
#define MANTISSA_SWAT_DATA_H

#include <fstream>
#include <vector>
#include <map>
#include <boost/mpi.hpp>

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

#include "MS_mpi_utils.h"

namespace MS {

    /*struct SWAT_row{
        int hru_num = 0;
        int year = 0;
        double irrtotal_mm = 0.0;
        double irrSW_mm = 0.0;
        double irrGW_mm = 0.0;
        double irrsaltSW_kgha = 0.0;
        double irrsaltGW_kgha = 0.0;
        double totpercsalt_kgha = 0.0;
        double perc_mm = 0.0;
    };*/

    class SWAT_data{
    public:
        SWAT_data(int n_hrus_in, int n_years_in)
        : n_hrus(n_hrus_in), n_yrears(n_years_in) {}
        bool read(const std::string &filename, int swat_version, boost::mpi::communicator &world);
        //bool read_v0(const std::string &filename, int NSwatYears, boost::mpi::communicator &world);
        bool read_v1(const std::string &filename, boost::mpi::communicator &world);
        bool read_HRU_idx_Map(const std::string &filename, boost::mpi::communicator &world);
        int hru_index(int HRU) const;


        //std::vector<SWAT_row> SWAT_TAB;
        //std::vector<std::vector<int>> HRUS;
        Matrix<double> irrtotal_mm;
        Matrix<double> irrSW_mm;
        Matrix<double> irrGW_mm;
        Matrix<double> irrsaltSW_kgha;
        Matrix<double> irrsaltGW_kgha;
        Matrix<double> fertsalt_kgha;
        Matrix<double> dssl_kgha;
        Matrix<double> Qsalt_kgha;
        Matrix<double> uptk_kgha;
        Matrix<double> pGW;
        Matrix<double> dSoilSalt_kgha;
        Matrix<double> totpercsalt_kgha;
        Matrix<double> perc_mm;
        Matrix<double> Salt_perc_ppm;
        Matrix<double> Trgt_AW_ppm;
        std::map<int,int> hru_idx_map;
    private:
        int n_hrus = 0;
        int n_yrears = 0;
        struct SWATField {
            const char* name;
            Matrix<double>* data;
        };
        std::vector<SWATField> v1_fields();
        //bool readASCIIset(const std::string &filename, Matrix<double> &data, int nHRU, int Nyrears);
        //bool readASCIIHRUS(std::string filename, Matrix<int> &data);
#if _USEHF > 0
        bool readHDF5PropertyDistrib(HighFive::File* file,
                                 const std::string& dataset_name,
                                 Matrix<double>& data,
                                 boost::mpi::communicator& world) const;
#endif
        bool checkRowsMatchHRU(const Matrix<double>& data,
                           const std::string& name,
                           boost::mpi::communicator& world) const;
    };

    inline std::vector<SWAT_data::SWATField> SWAT_data::v1_fields()
    {
        return {
            {"irrSW_mm",       &irrSW_mm},
            {"irrGW_mm",       &irrGW_mm},
            {"irrsaltSW_kgha", &irrsaltSW_kgha},
            {"irrsaltGW_kgha", &irrsaltGW_kgha},
            {"fertsalt_kgha",  &fertsalt_kgha},
            {"dssl_kgha",      &dssl_kgha},
            {"Qsalt_kgha",     &Qsalt_kgha},
            {"uptk_kgha",      &uptk_kgha},
            {"pGW",            &pGW},
            {"dSoilSalt_kgha", &dSoilSalt_kgha},
            {"perc_mm",        &perc_mm},
            {"Salt_perc_ppm",  &Salt_perc_ppm}
        };
    }


//     bool SWAT_data::read_v0(const std::string &filename, int NSwatYears, boost::mpi::communicator &world) {
//
//         std::string ext = getExtension(filename);
//
//         if (ext.compare("h5") == 0){
// #if _USEHF > 0
//                 const std::string HRUSNameSet("HRUS");
//                 const std::string irrtotal_mm_NameSet("irrtotal_mm");
//                 const std::string irrSW_mm_NameSet("irrSW_mm");
//                 const std::string irrGW_mm_NameSet("irrGW_mm");
//                 const std::string irrsaltSW_kgha_NameSet("irrsaltSW_kgha");
//                 const std::string irrsaltGW_Kgha_NameSet("irrsaltGW_kgha");
//                 const std::string totpercsalt_kgha_NameSet("totpercsalt_kgha");
//                 const std::string perc_mm_NameSet("perc_mm");
//                 const std::string Salt_perc_ppm_NameSet("Salt_perc_ppm");
//
//                 HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
//
//                 //HighFive::DataSet dataset_HRUS = HDFNfile.getDataSet(HRUSNameSet);
//                 HighFive::DataSet dataset_irrtotal_mm = HDFNfile.getDataSet(irrtotal_mm_NameSet);
//                 HighFive::DataSet dataset_irrSW_mm = HDFNfile.getDataSet(irrSW_mm_NameSet);
//                 HighFive::DataSet dataset_irrGW_mm = HDFNfile.getDataSet(irrGW_mm_NameSet);
//                 HighFive::DataSet dataset_irrsaltSW_kgha = HDFNfile.getDataSet(irrsaltSW_kgha_NameSet);
//                 HighFive::DataSet dataset_irrsaltGW_Kgha = HDFNfile.getDataSet(irrsaltGW_Kgha_NameSet);
//                 HighFive::DataSet dataset_totpercsalt_kgha = HDFNfile.getDataSet(totpercsalt_kgha_NameSet);
//                 HighFive::DataSet dataset_perc_mm = HDFNfile.getDataSet(perc_mm_NameSet);
//                 HighFive::DataSet dataset_Salt_perc_ppm = HDFNfile.getDataSet(Salt_perc_ppm_NameSet);
//
//                 //dataset_HRUS.read(HRUS);
//                 dataset_irrtotal_mm.read(irrtotal_mm);
//                 dataset_irrSW_mm.read(irrSW_mm);
//                 dataset_irrGW_mm.read(irrGW_mm);
//                 dataset_irrsaltSW_kgha.read(irrsaltSW_kgha);
//                 dataset_irrsaltGW_Kgha.read(irrsaltGW_kgha);
//                 dataset_totpercsalt_kgha.read(totpercsalt_kgha);
//                 dataset_perc_mm.read(perc_mm);
//                 dataset_Salt_perc_ppm.read(Salt_perc_ppm);
//
//                 return true;
// #endif
//         }
//         else{
//             //bool tf = RootReadsMatrixFileDistrib<int>(filename + "hrus.dat", HRUS,1,false,world);
//             //if (!tf){ return false;}
//             //int nHRUs = HRUS[0].size();
//
//             bool tf = RootReadsMatrixFileDistrib<double>(filename + "irrtotal_mm.dat", irrtotal_mm, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(irrtotal_mm,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "irrSW_mm.dat", irrSW_mm, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(irrSW_mm,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "irrGW_mm.dat", irrGW_mm, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(irrGW_mm,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltSW_kgha.dat", irrsaltSW_kgha, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(irrsaltSW_kgha,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "irrsaltGW_kgha.dat", irrsaltGW_kgha, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(irrsaltGW_kgha, world, 0, 10, 34, 40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "totpercsalt_kgha.dat", totpercsalt_kgha, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(totpercsalt_kgha,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "perc_mm.dat", perc_mm, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(perc_mm,world,0,10,34,40);}
//             tf = RootReadsMatrixFileDistrib<double>(filename + "Salt_perc_ppm.dat", Salt_perc_ppm, NSwatYears, true, world);
//             if (!tf){ return false;}
//             if (PrintMatrices){printMatrixForAllProc(Salt_perc_ppm,world,0,10,34,40);}
//
//             return true;
//         }
//         return false;
//     }

    bool SWAT_data::read_v1(const std::string& filename, boost::mpi::communicator &world) {
        const std::string ext = getExtension(filename);

        if (n_hrus <= 0) {
            if (world.rank() == 0) {
                std::cout << "SWAT_data::read_v1 error: hru_idx_map is empty." << std::endl;
            }
            return false;
        }

        if (n_yrears <= 0) {
            if (world.rank() == 0) {
                std::cout << "SWAT_data::read_v1 error: NSwatYears must be > 0." << std::endl;
            }
            return false;
        }

        const std::vector<SWATField> fields = v1_fields();

        if (ext == "h5" || ext == "H5"){
#if _USEHF > 0
            if (world.rank() == 0) {
                std::cout << "Reading HDF5 SWAT file " << filename << std::endl;
            }
            int fileOpenSuccess = 1;
            std::unique_ptr<HighFive::File> h5file;
            if (world.rank() == 0) {
                try {
                    h5file = std::make_unique<HighFive::File>(filename, HighFive::File::ReadOnly);
                }
                catch (const std::exception& e) {
                    std::cout << "Failed to open HDF5 file " << filename
                              << ": " << e.what() << std::endl;
                    fileOpenSuccess = 0;
                }
            }

            sendScalarFromRoot2AllProc(fileOpenSuccess, world);
            if (!fileOpenSuccess) {
                return false;
            }

            HighFive::File* file_ptr = (world.rank() == 0) ? h5file.get() : nullptr;

            for (const auto& field : fields) {
                const bool tf = readHDF5PropertyDistrib(file_ptr,
                                                        field.name,
                                                        *(field.data),
                                                        world);
                if (!tf) {
                    return false;
                }
                if (PrintMatrices) {
                    if (world.rank() == 0) {
                        std::cout << "Data for: " << field.name << std::endl;
                    }

                    printMatrixForAllProc(*(field.data), world, 46, 58, 1, 10);
                }
            }
            return true;
#else
            if (world.rank() == 0) {
                std::cout << "Cannot read HDF5 SWAT file because _USEHF == 0." << std::endl;
            }
            return false;
#endif
        }
        else{
            for (const auto& field : fields) {
                const std::string ascii_file = filename + field.name + ".dat";
                if (world.rank() == 0) {
                    std::cout << "Reading " << ascii_file << std::endl;
                }
                const bool tf = RootReadsMatrixFileDistrib<double>(ascii_file,
                                                                   *(field.data),
                                                                   n_yrears,
                                                                   world, 1000000, n_hrus);
                if (!tf) {
                    return false;
                }
                if (!checkRowsMatchHRU(*(field.data), ascii_file, world)) {
                    return false;
                }
                if (PrintMatrices) {
                    printMatrixForAllProc(*(field.data), world, 46, 58, 0, 10);
                }
            }
            return true;
        }
    }

    bool SWAT_data::read(const std::string &filename, int swat_version,
                         boost::mpi::communicator &world) {
        bool tf = false;
        // if (swat_version == 0){
        //     tf = read_v0(filename, NSwatYears,world);
        // }
        // else
        if(swat_version == 1){
            tf = read_v1(filename, world);
        }
        return tf;
    }

    bool SWAT_data::read_HRU_idx_Map(const std::string &filename, boost::mpi::communicator &world){
        std::vector<int> tmp;
        bool tf = RootReadsVectorFileDistrib<int>(filename, tmp, world,5000000, n_hrus);
        if (!tf){ return false;}
        if (PrintMatrices){printVectorForAllProc(tmp,world,10,30);}
        if (tmp.size() != n_hrus) {
            if (world.rank() == 0) {
                std::cout << "Error: the number of HRUs in the file " << filename << " does not match the number of HRUs in the model." << std::endl;
            }
            return false;
        }

        for (int i = 0; i < tmp.size(); ++i){
            hru_idx_map.insert(std::pair<int,int>(tmp[i], i));
        }
        return true;
    }

    int SWAT_data::hru_index(int HRU) const{
        if (HRU < 0){
            return -9;
        }
        std::map<int,int>::const_iterator  it = hru_idx_map.find(HRU);
        if (it != hru_idx_map.end()){
            return it->second;
        }
        else{
            return -9;
        }
    }

    // bool SWAT_data::readASCIIset(const std::string &filename, Matrix<double> &data, int nHRU, int Nyrears){
    //     data.allocate(Nyrears, nHRU);
    //
    //     std::ifstream ifile(filename);
    //     if (!ifile.is_open()){
    //         std::cout << "Cant open file: " << filename << std::endl;
    //         return false;
    //     }
    //     else{
    //         std::cout << "Reading " << filename << std::endl;
    //         std::string line;
    //         double val;
    //
    //         for (int i = 0; i < nHRU; ++i){
    //             getline(ifile, line);
    //             std::istringstream inp(line.c_str());
    //             for (int j = 0; j < Nyrears; ++j){
    //                 inp >> val;
    //                 data(j,i) = val;
    //             }
    //         }
    //         ifile.close();
    //         return true;
    //     }
    // }

    // bool SWAT_data::readASCIIHRUS(const std::string &filename, Matrix<int> &data){
    //     std::ifstream datafile(filename.c_str());
    //     if (!datafile.good()) {
    //         std::cout << "Can't open the file " << filename << std::endl;
    //         return false;
    //     }
    //     else{
    //         std::string line;
    //         std::vector<int> tmp;
    //         int v;
    //
    //         while (getline(datafile, line)){
    //             std::istringstream inp(line.c_str());
    //             inp >> v;
    //             tmp.push_back(v);
    //         }
    //         data.push_back(tmp);
    //         datafile.close();
    //     }
    //     return true;
    // }
#if _USEHF > 0
    inline bool SWAT_data::readHDF5PropertyDistrib(HighFive::File* file,
                                                   const std::string& dataset_name,
                                                   Matrix<double>& data,
                                                   boost::mpi::communicator& world) const {

        int readSuccess = 1;
        if (world.rank() == 0) {
            data.clear();
            try {
                if (file == nullptr) {
                    std::cout << "Internal error: null HDF5 file pointer on rank 0." << std::endl;
                    readSuccess = 0;
                }
                else if (!file->exist(dataset_name)) {
                    std::cout << "Missing HDF5 dataset: " << dataset_name << std::endl;
                    readSuccess = 0;
                }
                else {
                    HighFive::DataSet ds = file->getDataSet(dataset_name);
                    const std::vector<std::size_t> dims = ds.getSpace().getDimensions();

                    if (dims.size() != 2) {
                        std::cout << "Dataset " << dataset_name << " is not 2D" << std::endl;
                        readSuccess = 0;
                    }
                    else {
                        if (!((dims[0] == n_hrus && dims[1] == n_yrears) || (dims[1] == n_hrus && dims[0] == n_yrears))) {
                            std::cout << "The dataset dimensions are not [" << n_hrus << " x " << n_yrears << "] or ["
                                      << n_yrears << " x " << n_hrus << "] but [" << dims[0] << " x " << dims[1] << "]" << std::endl;
                            readSuccess = 0;
                        }
                        else {
                            std::vector<std::vector<double>> tmp;
                            ds.read(tmp);
                            data.allocate(n_hrus, n_yrears);
                            if (dims[0] == n_hrus && dims[1] == n_yrears) {
                                for (int i = 0; i < n_hrus; ++i) {
                                    for (int j = 0; j < n_yrears; ++j) {
                                        data(i, j) = tmp[i][j];
                                    }
                                }
                            }
                            else if (dims[0] == n_yrears && dims[1] == n_hrus) {
                                for (int i = 0; i < n_hrus; ++i) {
                                    for (int j = 0; j < n_yrears; ++j) {
                                        data(i, j) = tmp[j][i];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (const std::exception& e) {
                std::cout << "Exception while reading HDF5 dataset "
                      << dataset_name << ": " << e.what() << std::endl;
                readSuccess = 0;
            }
        }

        sendScalarFromRoot2AllProc(readSuccess, world);
        if (!readSuccess) {
            data.clear();
            return false;
        }

        sendMatrixFromRoot2AllProc(data, world);
        return true;
    }
#endif

    inline bool SWAT_data::checkRowsMatchHRU(const Matrix<double>& data,
                                             const std::string& name,
                                             boost::mpi::communicator& world) const
    {
        if (static_cast<int>(data.num_rows()) != n_hrus) {
            if (world.rank() == 0) {
                std::cout << name << " row count does not match hru_idx_map. "
                          << "Expected " << n_hrus << ", got " << data.num_rows() << std::endl;
            }
            return false;
        }
        return true;
    }
}

#endif //MANTISSA_SWAT_DATA_H
