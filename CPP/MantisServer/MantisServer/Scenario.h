//
// Created by giorgk on 5/8/2023.
//

#ifndef MANTISSERVER_SCENARIO_H
#define MANTISSERVER_SCENARIO_H

#include "MShelper.h"

namespace mantisServer {

/**
	 * @brief Scenario is a struct variable that contains all the information needed for each simulation scenario.
	 *
	 * This struct is a convenient way to transfer around all inputs at once.
	 */
    class Scenario {
    public:
        Scenario() {};

        void clear();

        bool parse_incoming_msg(std::string &msg, std::string &outmsg);


        std::string region;

        //! mapID is the id of selected background map.
        std::string mapID;

        //! flowScen is the code name of the selected flow scenario.
        std::string flowScen;

        //! This is a combination of Flow scenario and Welltype
        std::string flowWellScen;

        std::string wellType;

        bool bUseMonitoringWells;
        //bool bUseInitRFsets;
        //std::vector<std::string> RFSets;
        //std::string RFnameSet;

        //! The loading scenario from the list of initialization data
        std::string loadScen;
        int loadScenID;

        // Options for Raster loading type
        //! The name of Scenario Subset
        std::string loadSubScen;
        int loadSubScenID;

        //! The filename with the uploaded raster loading
        std::string modifierName;
        //! The operation type with respect to Scenario subscenario load option
        RasterOperation modifierType;
        //int modReplace;
        LinearData userRasterLoad;
        LoadUnits modifierUnit;

        ///! regionIDs is a list of regions to compute the breakthrough curves
        std::vector<std::string> regionIDs;

        ///! LoadReductionMap is a map that sets the nitrate loading reduction for selected land use categories.
        ///! The key value of the map is the land use category and the value is the percentage of loading reduction.
        std::map<int, double> LoadReductionMap;

        ///! If the input message contains crop id with -9 then this is applied to all crops except to the ones defined
        ///! in the LoadReductionMap
        double globalReduction;

        ///! ReductionYear is the year to start the reduction.
        int startReductionYear;

        ///! The end reduction is the year when the loading has reduced to the desired amount
        int endReductionYear;

        ///! This is a constant reduction coefficient that is applied over the simulation area on all loading types
        float constReduction;
        bool userSuppliedConstRed;

        //! This is the number of years to simulate starting from 1945
        int startSimulationYear;
        int endSimulationYear;

        int LoadTransitionStart;
        int LoadTransitionEnd;
        std::string LoadTransitionName;
        bool buseLoadTransition;

        //! This is the name of the unsaturated flow scenario
        std::string unsatScenario;
        int unsatScenarioID;

        //! This is the unsaturated zone mobile water content (m3/m3)
        //! Typical values are 0.05,0.1 ,0.15 and 0.20
        double unsatZoneMobileWaterContent;

        /**
         * This is an id for the current scenario. e.g. test001.
         * If this variable is present then during the simulation the following files will be printed
         * - debugID_urfs.dat contains the expanded unit respondse functions
         * - debugID_lf.dat containts the loading functions
         * - debugID_btc.dat contains the breakthrough curves. The convolution between urfs and lf
         * - debugID_well_btc.dat containts the averages BTC for the wells
         */
        std::string debugID;
        //! If the #debugID is present then this is set to true;
        bool printAdditionalInfo = false;

        /**
         * In GNLM the loading is defined as kg/ha. To convert to cocnentration we divide by the groundwater
         * recharge. However sometimes this has very low values due to errors in particle tracing.
         * In those cases letting the NO3 mas divided by a very small number results in unrealistically
         * large concentration. To avoid this we set loading to zero the the recharge is less than minRecharge.
         *
         * By defaults it gets 0.000027 which corresponds to 10 mm/year
         */
        double minRecharge;
        double maxConc;

        bool bUseFlowRch;
        std::string rchName;
        int rchScenID;

        int PixelRadius;

        radialSelection RadSelect;
        rectSelection RectSelect;
        bool useRadSelect;
        bool useRectSelect;

        numericRange DepthRange;
        numericRange ScreenLengthRange;
        bool useDepthRange;
        bool useScreenLenghtRange;
        bool bNarrowSelection;
        bool printWellIds;

        sourceArea SourceArea;
    };

/**
     * @brief clear is making sure that the scenario has no data from a previous run.
     *
     */
    void Scenario::clear() {
        region = "";
        mapID = "";
        flowScen = "";
        loadScen = "";
        loadSubScen = "";
        wellType = "";
        flowWellScen = "";
        modifierName = "";
        modifierType = RasterOperation::UNKNOWN;
        userRasterLoad.clear();
        regionIDs.clear();
        startReductionYear = 2020;
        endReductionYear = 2030;
        startSimulationYear = 1945;
        endSimulationYear = 2100;
        LoadReductionMap.clear();
        printAdditionalInfo = false;
        globalReduction = 1.0;
        constReduction = 1.0;
        userSuppliedConstRed = false;
        unsatZoneMobileWaterContent = 0.0;
        minRecharge = 0.000027; // 10 mm/year
        maxConc = 250;
        PixelRadius = 0;
        RadSelect.clear();
        RectSelect.clear();
        DepthRange.clear();
        ScreenLengthRange.clear();
        useRectSelect = false;
        useRadSelect = false;
        useDepthRange = false;
        useScreenLenghtRange = false;
        bNarrowSelection = false;
        buseLoadTransition = true;
        LoadTransitionName = "GNLM";
        LoadTransitionStart = 2005;
        LoadTransitionEnd = 2015;
        bUseMonitoringWells = false;
        bUseFlowRch = true;
        SourceArea.clear();
    }

    bool Scenario::parse_incoming_msg(std::string &msg, std::string &outmsg){
        bool out = false;
        clear();
        // convert the string to a stringstream
        std::istringstream ss(msg);

        int counter = 0;
        while (true){
            std::string test, tmp;
            counter++;
            ss >> test;
        }

        if (test == "Region") {
            ss >> region;
            continue;
        }

        if (test == "bMap") {
            ss >> mapID;
            continue;
        }

    }

}

#endif //MANTISSERVER_SCENARIO_H
