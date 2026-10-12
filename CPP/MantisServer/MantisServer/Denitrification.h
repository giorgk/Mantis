#ifndef MANTISSERVER_DENITRIFICATION_H
#define MANTISSERVER_DENITRIFICATION_H

#include <string>

#include "MShelper.h"
#include "BRaster.h"

namespace mantisServer {

    /**
     * Read-only collection of named denitrification maps.  This class only
     * supplies spatial values; simulation policy is intentionally handled by
     * the scenario and Region classes.
     */
    class DenitrificationMaps {
    public:
        bool readData(const std::string& filename, int nCells);
        bool hasMap(const std::string& name) const;
        int mapIndex(const std::string& name) const;
        bool getValue(int mapIndex, int activeCellIndex, double& value) const;
        bool getValue(int mapIndex, int row, int col,
                      const BackgroundRaster& raster, double& value) const;

    private:
        LinearData maps;
    };

    bool DenitrificationMaps::readData(const std::string& filename, int nCells) {
        maps.clear();
        maps.setNoDataValue(-1.0);
        if (!maps.readData(filename, nCells)) {
            return false;
        }
        int clampedValues = 0;
        if (!maps.clampRange(0.0, 1.0, clampedValues)) {
            std::cout << "Denitrification data in " << filename
                      << " have invalid dimensions or non-finite values"
                      << std::endl;
            return false;
        }
        if (clampedValues > 0) {
            std::cout << "WARNING: Clamped " << clampedValues
                      << " denitrification value(s) outside [0,1] in "
                      << filename << std::endl;
        }
        return true;
    }

    bool DenitrificationMaps::hasMap(const std::string& name) const {
        return maps.hasScenario(name);
    }

    int DenitrificationMaps::mapIndex(const std::string& name) const {
        return maps.ScenarioIndex(name);
    }

    bool DenitrificationMaps::getValue(int selectedMap, int activeCellIndex,
                                       double& value) const {
        value = maps.getValue(selectedMap, activeCellIndex);
        return std::isfinite(value) && value >= 0.0 && value <= 1.0;
    }

    bool DenitrificationMaps::getValue(int selectedMap, int row, int col,
                                       const BackgroundRaster& raster,
                                       double& value) const {
        return getValue(selectedMap, raster.IJ(row, col), value);
    }

}

#endif // MANTISSERVER_DENITRIFICATION_H
