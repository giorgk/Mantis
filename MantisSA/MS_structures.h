//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_STRUCTURES_H
#define MANTISSA_MS_STRUCTURES_H
#include <string>
#include <fstream>
namespace MS{

    struct SelectedWellsGroup{
        std::string groupName;
        std::vector<int> idVI;
        std::vector<int> idVD;
    };


    struct RasterOptions {
        int Ncells;
        int Nrows;
        int Ncols;
        double Xorig;
        double Yorig;
        double CellSize;
        std::string File;
    };

    struct SwatOptions{
        int Nyears;
        int StartYear;
        int version;
        std::string HRU_raster_file;
        std::string HRU_index_file;
        std::string Data_file;
    };

    struct HistoricOptions{
        int BlendStart;
        int BlendEnd;
        std::string filename;
        std::string ext;

    };

    struct NpsatOptions{
        bool bUseInitConcVI;
        bool bUseInitConcVD;
        //int version;
        std::string VIdataFile;
        std::string VDdataFile;
        std::string InitSaltVIFile;
        std::string InitSaltVDFile;
        std::string DistribuPumpFile;
    };

    struct RiverOptions{
        double StartDist;
        double EndDist;
        double ConcValue;
    };

    struct SimOptions{
        bool EnableFeedback;
        bool OutofAreaUseInitConc;
        int StartYear;
        int Nyears;
        int Porosity;
        int nBuffer;
        double MaxConc;
        double MaxAge;
        double SurfConcValue;
        double OutofAreaConc;
        int nYears_historic;
        int nYears_blendEnd;

    };

    struct UnsatOptions{
        double wc;
        double minDepth;
        double minRch;
        std::string Depth_file;
        std::string Depth_name;
        std::string Rch_file;
    };

    struct OutputOptions{
        bool printLoad;
        bool printURFs;
        bool printBTCs;
        bool printSelectedWells;
        bool compress;
        std::string OutFile;
        std::string SelectedWells;
        std::string SelectedWellsGroups;
    };

    struct MiscOptions {
        double dp_mult;
    };

    struct STRML{
        int Sid;
        int urfI;
        int urfJ;
        int IJ;
        bool inRiv;
        double rivInfl;
        //int hru_idx;
        double W;
        double Len;
        double m;
        double s;
        double a;
        int wellSourceId;
        std::vector<double> urf;
        std::vector<double> lf_conc;
        //std::vector<double> lf_mass;
        std::vector<double> gw_conc;
        //std::vector<double> gw_conc;
        std::vector<double> btc;

    };

    struct WELL{
        double initConc = 0.0;
        std::vector<STRML> strml;
        std::vector<double> wellBtc;
    };

    struct NPSATTMP{
        double sumW = 0;
        std::vector<int> Sid;
        std::vector<int> urfI;
        std::vector<int> urfJ;
        std::vector<bool> inRiv;
        std::vector<double> rivDist;
        std::vector<int> hru_idx;
        std::vector<double> W;
        std::vector<double> Len;
        std::vector<double> m;
        std::vector<double> s;
        std::vector<double> a;
    };

    struct WELLS{
        std::map<int, WELL> wells;
        int version;
    };
    typedef std::map<int,NPSATTMP> tmpWELLS;

    //typedef std::map<int, WELL> WELLS;

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

    inline bool readSelectedWellsGroupInfo(std::string filename, std::map<int,SelectedWellsGroup>& swg_map,
                                           boost::mpi::communicator &world){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }

        swg_map.clear();

        if (world.rank() == 0){
            std::cout << "Reading " << filename << std::endl;
        }

        std::string line;
        while (std::getline(datafile, line)) {
            const std::size_t first = line.find_first_not_of(" \t\r\n");
            // skip blank / whitespace-only lines
            if (first == std::string::npos) {
                continue;
            }

            // skip comment lines (including lines with leading spaces before #)
            if (line[first] == '#') {
                continue;
            }

            std::istringstream inp(line);
            int groupId;
            std::string groupName;

            if (!(inp >> groupId >> groupName)) {
                if (world.rank() == 0) {
                    std::cout << "Invalid format in file " << filename
                              << ". Expected: ID name" << std::endl;
                }
                return false;
            }

            if (swg_map.find(groupId) != swg_map.end()) {
                if (world.rank() == 0) {
                    std::cout << "Selected Wells Group with id "
                              << groupId << " was found more than once" << std::endl;
                }
                return false;
            }

            SelectedWellsGroup swg;
            swg.groupName = groupName;
            swg_map[groupId] = swg;
        }
        return true;
    }

    template <typename T>
    bool parse_value_space(const char*& p, const char* end, T& value);

    template <>
    inline bool parse_value_space<int>(const char*& p, const char* end, int& value)
    {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
        if (p >= end) return false;

#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
        auto res = std::from_chars(p, end, value);
        if (res.ec != std::errc()) return false;
        p = res.ptr;
        return true;
#else
        char* q = nullptr;
        long v = std::strtol(p, &q, 10);
        if (q == p) return false;
        if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) return false;
        value = static_cast<int>(v);
        p = q;
        return true;
#endif
    }

    template <>
    inline bool parse_value_space<double>(const char*& p, const char* end, double& value)
    {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
        if (p >= end) return false;

        char* q = nullptr;
        value = std::strtod(p, &q);
        if (q == p) return false;
        p = q;
        return true;
    }

    template<typename T>
    bool readMatrix(const std::string& filename, std::vector<std::vector<T>>& data, int nCols, int freq = 500000)
    {
        if (nCols <= 0) {
            std::cout << "readMatrix: nCols must be > 0" << std::endl;
            return false;
        }
        if (freq <= 0) {
            freq = 500000;
        }

        std::ifstream datafile(filename);
        if (!datafile.is_open()) {
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }

        std::cout << "Reading " << filename << std::endl;

        data.clear();

        std::string line;
        int lineCount = 0;
        int nextPrint = (freq > 0) ? freq : std::numeric_limits<int>::max();

        while (std::getline(datafile, line))
        {
            const char* p = line.data();
            const char* end = p + line.size();

            std::vector<T> row(static_cast<std::size_t>(nCols));

            for (int j = 0; j < nCols; ++j) {
                if (!parse_value_space<T>(p, end, row[static_cast<std::size_t>(j)])) {
                    std::cout << "Error reading file " << filename
                              << " at data row " << lineCount
                              << ", column " << j << std::endl;
                    return false;
                }
            }

            // Allow trailing spaces/tabs only
            while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
            if (p != end) {
                std::cout << "Error reading file " << filename
                          << " at data row " << lineCount
                          << ": expected " << nCols
                          << " columns, got more" << std::endl;
                return false;
            }

            data.push_back(std::move(row));
            ++lineCount;

            if (lineCount >= nextPrint) {
                std::cout << lineCount << " lines read..." << std::endl;
                nextPrint += freq;
            }
        }
        if (data.empty()) {
            std::cout << "Error: matrix file " << filename
                      << " contains no valid data rows" << std::endl;
            return false;
        }
        return true;
    }

    template<typename T>
    void transposeMatrix(const std::vector<std::vector<T>> & data_in, std::vector<std::vector<T>> & data_out){
        data_out.clear();
        if (data_in.empty()) {
            return;
        }

        const std::size_t nRows = data_in.size();
        const std::size_t nCols = data_in[0].size();

        // Check rectangular matrix
        for (std::size_t i = 1; i < nRows; ++i) {
            if (data_in[i].size() != nCols) {
                throw std::runtime_error("transposeMatrix: input is not rectangular");
            }
        }

        data_out.assign(nCols, std::vector<T>(nRows));

        for (std::size_t i = 0; i < nRows; ++i) {
            const auto& row_in = data_in[i];
            for (std::size_t j = 0; j < nCols; ++j) {
                data_out[j][i] = row_in[j];
            }
        }
    }


    bool readLinearData(std::string filename,
                        std::vector<std::string> &names,
                        std::vector<std::vector<double>> &data,
                        int freq = 500000){

        std::ifstream ifile(filename);

        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }

        std::cout << "Reading " << filename << " ..." << std::endl;
        std::string line;
        int nData = 0;

        // --- Read number of data series
        if (!std::getline(ifile, line)) {
            std::cerr << "File format error: missing nData" << std::endl;
            return false;
        }
        nData = std::stoi(line);
        data.clear();
        data.resize(nData);
        names.clear();
        names.reserve(nData);

        // --- Read names
        for (int i = 0; i < nData; ++i){
            if (!std::getline(ifile, line)){
                std::cerr << "File format error: missing name at line " << i + 2 << std::endl;
                return false;
            }
            names.push_back(line);  // assume each name is a line
        }

        // --- Read values
        int ir = 0;
        int nextLog = freq;
        std::vector<double> rowValues(nData);
        while (std::getline(ifile, line)){
            if (line.empty()) continue;

            std::istringstream iss(line);
            for (int j = 0; j < nData; ++j){
                if (!(iss >> rowValues[j])) {
                    std::cerr << "Parsing error at row " << ir << ", column " << j << std::endl;
                    return false;
                }
            }
            for (int j = 0; j < nData; ++j) {
                data[j].push_back(rowValues[j]);
            }

            ++ir;
            if (ir >= nextLog) {
                std::cout << ir << " rows read..." << std::endl;
                nextLog += freq;
            }
        }
        return true;
    }

    inline double calcRiverInfluence(double dist, const RiverOptions &ropt){
        const double d0 = ropt.StartDist;
        const double d1 = ropt.EndDist;

        // Handle degenerate or invalid cases
        if (!(d1 > d0)) {
            // If equal or invalid ordering, treat as step function
            return (dist <= d0) ? 1.0 : 0.0;
        }

        if (dist <= d0) return 1.0;
        if (dist >= d1) return 0.0;

        const double u = 1.0 - (dist - d0) / (d1 - d0);

        // Optional clamp (very cheap, avoids numerical drift)
        return std::max(0.0, std::min(1.0, u));
    }

    inline bool readNPSATinfo(const std::string& filename, std::vector<int>& M){
        std::ifstream datafile(filename);
        if (!datafile.good()){
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }

        std::cout << "Reading " << filename << std::endl;

        M.clear();

        std::string line;
        std::string propname;

        int version = -1;
        bool foundVersion = false;
        bool foundPor = false;
        bool foundINT = false;
        bool foundDBL = false;
        bool foundMSA = false;

        std::vector<int> por_values;
        std::vector<int> int_size;
        std::vector<int> dbl_size;
        std::vector<int> msa_size;

        int lineCount = 0;

        while (getline(datafile, line)) {
            ++lineCount;

            const std::size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) {
                continue; // blank line
            }
            if (line[first] == '#') {
                continue; // comment line
            }

            std::istringstream inp(line);
            if (!(inp >> propname)) {
                continue;
            }

            if (propname == "Version") {
                if (!(inp >> version)) {
                    std::cout << "Invalid Version line in " << filename
                              << " at line " << lineCount << std::endl;
                    return false;
                }
                foundVersion = true;
                continue;
            }

            if (propname == "Por") {
                int value;
                por_values.clear();
                while (inp >> value) {
                    por_values.push_back(value);
                }
                if (por_values.empty()) {
                    std::cout << "Invalid Por line in " << filename
                              << " at line " << lineCount << std::endl;
                    return false;
                }
                foundPor = true;
                continue;
            }

            if (propname == "INT") {
                int value;
                int_size.clear();
                while (inp >> value) {
                    int_size.push_back(value);
                }
                if (int_size.size() != 2) {
                    std::cout << "INT line in " << filename
                              << " must contain exactly 2 integers" << std::endl;
                    return false;
                }
                foundINT = true;
                continue;
            }

            if (propname == "DBL") {
                int value;
                dbl_size.clear();
                while (inp >> value) {
                    dbl_size.push_back(value);
                }
                if (dbl_size.size() != 2) {
                    std::cout << "DBL line in " << filename
                              << " must contain exactly 2 integers" << std::endl;
                    return false;
                }
                foundDBL = true;
                continue;
            }

            if (propname == "MSA") {
                int value;
                msa_size.clear();
                while (inp >> value) {
                    msa_size.push_back(value);
                }
                if (msa_size.size() != 2) {
                    std::cout << "MSA line in " << filename
                              << " must contain exactly 2 integers" << std::endl;
                    return false;
                }
                foundMSA = true;
                continue;
            }

            std::cout << "Unknown property name \"" << propname
                  << "\" in " << filename
                  << " at line " << lineCount << std::endl;
            return false;
        }

        if (!foundVersion) {
            std::cout << "Missing Version line in " << filename << std::endl;
            return false;
        }
        if (!foundPor) {
            std::cout << "Missing Por line in " << filename << std::endl;
            return false;
        }
        if (!foundINT) {
            std::cout << "Missing INT line in " << filename << std::endl;
            return false;
        }
        if (!foundDBL) {
            std::cout << "Missing DBL line in " << filename << std::endl;
            return false;
        }
        if (!foundMSA) {
            std::cout << "Missing MSA line in " << filename << std::endl;
            return false;
        }

        M.push_back(version);
        M.push_back(static_cast<int>(por_values.size()));
        M.insert(M.end(), por_values.begin(), por_values.end());
        M.insert(M.end(), int_size.begin(), int_size.end());
        M.insert(M.end(), dbl_size.begin(), dbl_size.end());
        M.insert(M.end(), msa_size.begin(), msa_size.end());

        return true;
    }

    int calcYearIndex(const int &yr, const int &startYr, const int &nYrs){
        int ii;
        if (yr >= startYr){
            ii = 0;
            for (int i = 0; i < 1000; ++i){
                if (startYr + i == yr){
                    break;
                }
                ii = ii + 1;
                if (ii > nYrs){
                    ii = 0;
                }
            }
        }
        else{
            ii = nYrs;
            for (int i = -1; i > -1000; --i){
                if (startYr + i == yr){
                    break;
                }
                ii = ii - 1;
                if (ii < 0){
                    ii = nYrs;
                }
            }
        }
        return ii;
    }

    inline void yearMap(std::map<int, int> &yr_id, const int &startSimYr, const int &endSimYr, const int &startYr, const int &nYrs){
        yr_id.clear();

        if (nYrs <= 0) {
            throw std::invalid_argument("yearMap: nYrs must be > 0");
        }
        if (endSimYr < startSimYr) {
            throw std::invalid_argument("yearMap: endSimYr < startSimYr");
        }

        for (int y = startSimYr; y <= endSimYr; ++y) {
            int idx = (y - startYr) % nYrs;
            if (idx < 0) {
                idx += nYrs;
            }
            yr_id[y] = idx;
        }
    }

    bool copy_file(std::string src, std::string dst){
        std::ifstream src_file(src.c_str(), std::ios::binary);
        std::ofstream dst_file(dst.c_str(), std::ios::binary);
        if (!src_file.is_open()) {
            std::cout << "Cant open file: " << src << std::endl;
            return false;
        }
        if (!dst_file.is_open()) {
            std::cout << "Cant open file: " << dst << std::endl;
            return false;
        }

        try {
            dst_file << src_file.rdbuf();
        }
        catch (std::exception &e) {
            std::cout << "Error: " << e.what() << std::endl;
            return false;
        }
        return true;
    }

}

#endif //MANTISSA_MS_STRUCTURES_H
