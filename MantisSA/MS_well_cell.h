//
// Created by giorgk on 3/26/2026.
//

#ifndef MANTISSA_MS_WELL_CELL_H
#define MANTISSA_MS_WELL_CELL_H

#include <unordered_map>

namespace MS {
    class WELL_CELLS {
    public:
        WELL_CELLS(){}
        bool read_build(const std::string& filename,
                        boost::mpi::communicator& world);
        void clear()
        {
            well_to_pos.clear();
            offsets.clear();
            cells.clear();
        }

        bool empty() const
        {
            return well_to_pos.empty();
        }

        std::size_t n_wells() const
        {
            return well_to_pos.size();
        }

        std::size_t n_cells() const
        {
            return cells.size();
        }

        bool has_well(int wellid) const
        {
            return well_to_pos.find(wellid) != well_to_pos.end();
        }

        // Efficient access: returns [begin, end) pointers into internal storage.
        // Returns false if wellid does not exist.
        bool get_cells_ptr(int wellid, const int*& begin, const int*& end) const
        {
            begin = nullptr;
            end   = nullptr;

            const std::unordered_map<int, int>::const_iterator it = well_to_pos.find(wellid);
            if (it == well_to_pos.end()) {
                return false;
            }

            const int pos = it->second;
            const int ibeg = offsets[static_cast<std::size_t>(pos)];
            const int iend = offsets[static_cast<std::size_t>(pos + 1)];

            begin = cells.data() + ibeg;
            end   = cells.data() + iend;
            return true;
        }

        // Convenience method: returns a copy of the cell list.
        std::vector<int> get_cells_copy(int wellid) const
        {
            const int* begin = nullptr;
            const int* end   = nullptr;

            if (!get_cells_ptr(wellid, begin, end)) {
                return std::vector<int>();
            }

            return std::vector<int>(begin, end);
        }

    private:
        std::unordered_map<int, int> well_to_pos; // well id -> compressed position
        std::vector<int> offsets;                 // size Nwells + 1
        std::vector<int> cells;                   // flattened cell indices

        bool read_ascii(const std::string& filename, std::vector<int>& cellWell, int freq = 5000000);
        bool read_hdf5(const std::string& filename, std::vector<int>& cellWell);
        bool build_from_cellWell(const std::vector<int>& cellWell);
        void print_progress_10(std::size_t i,
                                  std::size_t n,
                                  std::size_t& next_mark,
                                  int& next_pct);
    };

    inline bool WELL_CELLS::read_build(const std::string& filename, boost::mpi::communicator& world) {
        clear();
        if (world.rank() != 0) {
            return true;
        }

        std::vector<int> cellWell;
        const std::string ext = getExtension(filename);

        bool tf = false;
        if (ext == "h5" || ext == "H5") {
            tf = read_hdf5(filename, cellWell);
        }
        else {
            tf = read_ascii(filename, cellWell);
        }

        if (!tf) {
            clear();
            return false;
        }

        tf = build_from_cellWell(cellWell);
        if (!tf) {
            clear();
            return false;
        }

        return true;
    }

    inline bool WELL_CELLS::read_ascii(const std::string& filename, std::vector<int>& cellWell, int freq) {
        std::ifstream fin(filename.c_str());
        if (!fin.is_open()) {
            std::cout << "Cannot open file " << filename << std::endl;
            return false;
        }

        std::cout << "Reading " << filename << std::endl;

        cellWell.clear();

        std::string line;
        int value = 0;
        std::size_t lineCount = 0;
        std::size_t nextPrint = freq;

        while (std::getline(fin, line)) {
            std::istringstream iss(line);

            if (!(iss >> value)) {
                std::cout << "Error reading " << filename
                          << " at row " << lineCount << std::endl;
                return false;
            }

            cellWell.push_back(value);
            ++lineCount;

            if (lineCount >= nextPrint) {
                std::cout << lineCount << " lines read..." << std::endl;
                nextPrint += freq;
            }
        }

        if (cellWell.empty()) {
            std::cout << "File " << filename << " contains no data." << std::endl;
            return false;
        }

        return true;
    }

    inline bool WELL_CELLS::read_hdf5(const std::string& filename, std::vector<int>& cellWell) {
#if _USEHF > 0
        try {
            std::cout << "Reading " << filename << std::endl;
            HighFive::File file(filename, HighFive::File::ReadOnly);

            if (!file.exist("WellID")) {
                std::cout << "Dataset WellID does not exist in " << filename << std::endl;
                return false;
            }

            HighFive::DataSet ds = file.getDataSet("WellID");
            ds.read(cellWell);

            if (cellWell.empty()) {
                std::cout << "Dataset WellID in " << filename << " is empty." << std::endl;
                return false;
            }

            return true;
        }
        catch (const std::exception& e) {
            std::cout << "Error reading HDF5 file " << filename
                      << ": " << e.what() << std::endl;
            return false;
        }
#else
        (void)filename;
        (void)cellWell;
        std::cout << "Cannot read HDF5 file because _USEHF == 0." << std::endl;
        return false;
#endif
    }

    inline bool WELL_CELLS::build_from_cellWell(const std::vector<int>& cellWell)
    {
        if (cellWell.empty()) {
            std::cout << "build_from_cellwell: input vector is empty." << std::endl;
            return false;
        }

        const std::size_t n = cellWell.size();
        // -------------------------------
        // Pass 1: count how many cells belong to each well
        // -------------------------------
        std::cout << "Building well-cell structure: counting wells..." << std::endl;

        std::unordered_map<int, int> counts;
        counts.reserve(n);

        std::size_t next_mark = n / 10;
        int next_pct = 1;

        for (std::size_t i = 0; i < n; ++i) {
            print_progress_10(i, n, next_mark, next_pct);
            ++counts[cellWell[i]];
        }
        if (n >= 10) {
            std::cout << "100%";
        }
        std::cout << std::endl;

        // -------------------------------
        // Assign compressed positions
        // -------------------------------
        std::cout << "Building well-cell structure: assigning well positions..." << std::endl;

        well_to_pos.clear();
        well_to_pos.reserve(counts.size());

        int pos = 0;
        for (std::unordered_map<int, int>::const_iterator it = counts.begin();
             it != counts.end(); ++it)
        {
            well_to_pos[it->first] = pos;
            ++pos;
        }

        const int nWells = static_cast<int>(counts.size());

        // -------------------------------
        // Build offsets
        // -------------------------------
        std::cout << "Building well-cell structure: building offsets..." << std::endl;
        offsets.assign(static_cast<std::size_t>(nWells + 1), 0);

        for (std::unordered_map<int, int>::const_iterator it = counts.begin();
         it != counts.end(); ++it)
        {
            const int p = well_to_pos[it->first];
            offsets[static_cast<std::size_t>(p + 1)] = it->second;
        }

        for (int i = 0; i < nWells; ++i) {
            offsets[static_cast<std::size_t>(i + 1)] += offsets[static_cast<std::size_t>(i)];
        }

        // -------------------------------
        // Allocate flat array
        // -------------------------------
        cells.assign(n, -1);

        // Working write positions
        std::vector<int> next = offsets;

        // -------------------------------
        // Pass 2: fill cell indices
        // -------------------------------
        std::cout << "Building well-cell structure: filling cell lists..." << std::endl;
        next_mark = n / 10;
        next_pct = 1;

        for (std::size_t i = 0; i < n; ++i) {
            print_progress_10(i, n, next_mark, next_pct);

            const int wellid = cellWell[i];
            const std::unordered_map<int, int>::const_iterator it = well_to_pos.find(wellid);

            if (it == well_to_pos.end()) {
                std::cout << "Internal error: well id " << wellid
                          << " not found during fill." << std::endl;
                return false;
            }

            const int p = it->second;
            cells[static_cast<std::size_t>(next[static_cast<std::size_t>(p)]++)] =
                static_cast<int>(i);
        }
        if (n >= 10) {
            std::cout << "100%";
        }
        std::cout << std::endl;

        std::cout << "Finished building well-cell structure for "
              << nWells << " wells and " << n << " cells." << std::endl;

        return true;
    }

    inline void WELL_CELLS::print_progress_10(std::size_t i,
                                              std::size_t n,
                                              std::size_t& next_mark,
                                              int& next_pct)
    {
        if (n < 10) {
            return;
        }

        while (i >= next_mark && next_pct <= 100) {
            std::cout << next_pct << "% " << std::flush;
            ++next_pct;
            next_mark = (n * static_cast<std::size_t>(next_pct)) / 10;
        }
    }

}

#endif //MANTISSA_MS_WELL_CELL_H