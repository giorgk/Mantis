//
// Created by giorgk on 5/8/2023.
//

#ifndef MANTISSERVER_SCENARIO_H
#define MANTISSERVER_SCENARIO_H

#include <string>
#include <string_view>

#include "MShelper.h"

namespace mantisServer {

    struct InitConcParam{
        double mean = 1;
        double std = 1;
        double minv = 1;
        double maxv = 1;
        bool bUse = false;
        void reset(){
            mean = 1;
            std = 1;
            minv = 1;
            maxv = 1;
            bUse = false;
        }
    };

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
        bool parse_incoming_msg2(std::string &msg, std::string &outmsg);


        std::string modelArea;

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

        int porosity = 10;
        int porosityIndex = 0;

        //! This is the unsaturated zone mobile water content (m3/m3)
        //! Typical values are 0.05,0.1 ,0.15 and 0.20
        double unsatZoneMobileWaterContent;

        double unsatMinDepth = 1.0;
        double unsatDepthOffset = 0.0;
        double unsatDepthMult = 1.0;
        double unsatMinRch = 10.0;
        double unsatRchOffset = 0.0;
        double unsatRchMult = 1;
        double maxAge = 400;

        /**
         * This is an id for the current scenario. e.g. test001.
         * If this variable is present then during the simulation the following files will be printed
         * - debugID_urfs.dat contains the expanded unit respondse functions
         * - debugID_lf.dat containts the loading functions
         * - debugID_btc.dat contains the breakthrough curves. The convolution between urfs and lf
         * - debugID_well_btc.dat containts the averages BTC for the wells
         */
        std::string debugID;
        std::string debugPath;
        bool printLF = false;
        bool printBTC = false;
        bool printURF = false;
        bool printWellBTC = false;
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
        double surfConc = 0.0;
        double surfRivDist = 50.0;
        double surfRivInfl = 200.0;


        bool bUseFlowRch;
        std::string rchName;
        int rchScenID;

        int PixelRadius;

        radialSelection RadSelect;
        rectSelection RectSelect;
        bool useRadSelect;
        bool useRectSelect;

        numericRange DepthRange; // Depth is the total depth from GSE to the bottom of the screen
        numericRange wt2tRange;
        numericRange ScreenLengthRange;
        numericRange unsatRange;

        bool bNarrowSelection;
        bool printWellIds;
        bool printWellinfo;

        sourceArea SourceArea;
        int maxSourceCells;
        URFTYPE urfType = URFTYPE::LGNRM;
        double adeLambda = 0.0;
        double adeR = 1.0;
        InitConcParam initcondparam;
    };

/**
     * @brief clear is making sure that the scenario has no data from a previous run.
     *
     */
    void Scenario::clear() {
        modelArea = "";
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
        globalReduction = 1.0;
        constReduction = 1.0;
        userSuppliedConstRed = false;
        unsatZoneMobileWaterContent = 0.0;
        minRecharge = 10;// mm/year
        maxConc = 250;
        surfConc = 0.0;
        surfRivDist = 50.0;
        surfRivInfl = 200.0;
        maxAge = 400;
        PixelRadius = 0;
        RadSelect.clear();
        RectSelect.clear();
        DepthRange.clear();
        ScreenLengthRange.clear();
        useRectSelect = false;
        useRadSelect = false;
        bNarrowSelection = false;
        buseLoadTransition = true;
        LoadTransitionName = "NONE";
        LoadTransitionStart = 2005;
        LoadTransitionEnd = 2015;
        bUseMonitoringWells = false;
        bUseFlowRch = true;
        porosity = 10;
        porosityIndex = 0;
        SourceArea.clear();
        printWellIds = false;
        printWellinfo = false;
        debugID = "";
        printLF = false;
        printURF = false;
        printBTC = false;
        printWellBTC = false;
        printAdditionalInfo = false;
        maxSourceCells = 1000;
        urfType = URFTYPE::LGNRM;
        adeLambda = 0.0;
        adeR = 1.0;
        initcondparam.reset();
    }

    bool Scenario::parse_incoming_msg(std::string &msg, std::string &outmsg){
        bool out = false;
        clear();
        // convert the string to a stringstream
        std::istringstream ss(msg);

        int counter = 0;
        while (true){
            std::string test, tmp;
            counter = counter + 1;
            ss >> test;

            if (test == "modelArea") {
                ss >> modelArea;
                continue;
            }
            if (test == "bMap") {
                ss >> mapID;
                continue;
            }
            if (test == "Nregions") {
                int Nregions;
                ss >> Nregions;
                for (int i = 0; i < Nregions; ++i) {
                    ss >> test;
                    regionIDs.push_back(test);
                }
                continue;
            }

            if (test == "flowScen") {
                ss >> flowScen;
                continue;
            }

            if (test == "wellType"){
                ss >> wellType;
                continue;
            }
            if (test == "unsatScen") {
                ss >> unsatScenario;
                continue;
            }
            if (test == "unsatWC") {
                ss >> unsatZoneMobileWaterContent;
                continue;
            }
            if (test == "rchMap"){
                bUseFlowRch = false;
                ss >> rchName;
                continue;
            }
            if (test == "minRch") {
                ss >> minRecharge;
                continue;
            }
            if (test == "startSimYear") {
                ss >> startSimulationYear;
                continue;
            }
            if (test == "endSimYear") {
                ss >> endSimulationYear;
                continue;
            }
            if (test == "startRed" ) {
                ss >> startReductionYear;
                continue;
            }
            if (test == "endRed") {
                ss >> endReductionYear;
                continue;
            }
            if (test == "loadScen") {
                ss >> loadScen;
                continue;
            }
            if (test == "loadSubScen") {
                ss >> loadSubScen;
                continue;
            }
            if (test == "modifierName") {
                ss >> modifierName;
                continue;
            }
            if (test == "modifierType") {
                ss >> tmp;
                modifierType = string2RasterOperation(tmp);
                continue;
            }
            if (test == "modifierUnit") {
                ss >> tmp;
                modifierUnit = string2LoadUnits(tmp);
                continue;
            }
            if (test == "Ncrops") {
                int Ncrops, cropid;
                double perc;
                ss >> Ncrops;
                for (int i = 0; i < Ncrops; ++i) {
                    ss >> cropid;
                    ss >> perc;
                    if (cropid == -9){
                        globalReduction = perc;
                    }
                    else{
                        LoadReductionMap.insert(std::pair<int, double>(cropid, perc));
                    }
                }
                continue;
            }

            if (test == "maxConc") {
                ss >> maxConc;
                continue;
            }

            if (test == "surfConc") {
                ss >> surfConc;
                continue;
            }

            if (test == "surfRivDist") {
                ss >> surfRivDist;
                continue;
            }

            if (test == "surfRivInfl") {
                ss >> surfRivInfl;
                continue;
            }

            if (test == "maxAge"){
                ss >> maxAge;
                continue;
            }

            if (test == "constRed"){
                userSuppliedConstRed = true;
                ss >> constReduction;
                continue;
            }

            if (test == "loadTrans"){
                ss >> LoadTransitionName;
                ss >> LoadTransitionStart;
                ss >> LoadTransitionEnd;
                if (LoadTransitionName.compare("NONE") == 0){
                    buseLoadTransition = false;
                }
                else{
                    buseLoadTransition = true;
                }

                continue;
            }

            if (test == "por"){
                ss >> porosity;
                continue;
            }

            if (test == "urfType"){
                ss >> tmp;
                urfType = string2URFTYPE(tmp);
                continue;
            }

            if (test == "ADELR"){
                ss >> adeLambda;
                ss >> adeR;
                continue;
            }

            if (test == "maxSourceCells") {
                ss >> maxSourceCells;
                continue;
            }

            if (test == "RadSelect"){
                bNarrowSelection = true;
                useRadSelect = true;
                double cx, cy, r;
                ss >> cx;
                ss >> cy;
                ss >> r;
                RadSelect.setData(cx,cy,r);
                continue;
            }

            if (test == "RectSelect"){
                bNarrowSelection = true;
                useRectSelect = true;
                double xmin, ymin, xmax, ymax;
                ss >> xmin;
                ss >> ymin;
                ss >> xmax;
                ss >> ymax;
                RectSelect.setData(xmin, ymin, xmax, ymax);
                continue;
            }

            if (test == "DepthRange"){
                bNarrowSelection = true;
                double dmin, dmax;
                ss >> dmin;
                ss >> dmax;
                DepthRange.setData(dmin, dmax);
                continue;
            }

            if (test == "UnsatRange"){
                bNarrowSelection = true;
                double umin, umax;
                ss >> umin;
                ss >> umax;
                unsatRange.setData(umin, umax);
                continue;
            }

            if (test == "Wt2tRange"){
                bNarrowSelection = true;
                double wtmin, wtmax;
                ss >> wtmin;
                ss >> wtmax;
                wt2tRange.setData(wtmin, wtmax);
                continue;
            }

            if (test == "ScreenLenRange"){
                bNarrowSelection = true;
                double slmin, slmax;
                ss >> slmin;
                ss >> slmax;
                ScreenLengthRange.setData(slmin, slmax);
                continue;
            }

            if (test == "SourceArea"){//TODO
                int nPix, minPix, maxPix;
                double percPix;
                ss >> nPix;
                ss >> minPix;
                ss >> maxPix;
                ss >> percPix;
                SourceArea.setParameters(nPix, minPix, maxPix, percPix);
                continue;
            }
            if (test ==  "PixelRadius") {//TODO
                ss >> PixelRadius;
                continue;
            }
            if (test == "getids"){
                int tmp;
                ss >> tmp;
                printWellIds = tmp != 0;
                continue;
            }

            if (test == "getotherinfo"){
                int tmp;
                ss >> tmp;
                printWellinfo = tmp != 0;
                continue;
            }

            if (test == "initConcParam"){
                ss >> initcondparam.mean;
                ss >> initcondparam.std;
                ss >> initcondparam.minv;
                ss >> initcondparam.maxv;
                if (initcondparam.mean < -90 || initcondparam.std < 0 || initcondparam.minv < 0 || initcondparam.maxv < 0){
                    initcondparam.reset();
                }
                else{
                    initcondparam.bUse = true;
                }

                continue;
            }

            if (test ==  "DebugID") {
                ss >> debugID;
                //if (!debugID.empty())
                //    printAdditionalInfo = true;
                continue;
            }
            if (test == "printLF"){
                int tmp;
                ss >> tmp;
                printLF = tmp != 0;
                continue;
            }
            if (test == "printURF"){
                int tmp;
                ss >> tmp;
                printURF = tmp != 0;
                continue;
            }
            if (test == "printBTC"){
                int tmp;
                ss >> tmp;
                printBTC = tmp != 0;
                continue;
            }
            if (test == "printWellBTC"){
                int tmp;
                ss >> tmp;
                printWellBTC = tmp != 0;
                continue;
            }
            /*
            if (test == "RFset"){
                scenario.bUseRuntimeRFsets = true;
                scenario.bUseInitRFsets = false;
                int Nsets;
                ss >> Nsets;
                for (int i = 0; i < Nsets; ++i){
                    ss >> test;
                    scenario.RFSets.push_back(test);
                }
            }
            */
            if (test == "ENDofMSG") {
                out = true;
                break;
            }
            if (test.empty()) {
                outmsg += "0 ERROR: Empty message was found\n";
                out = false;
                break;
            }
            if (counter > 500) {
                outmsg += "0 ERROR: After 200 iterations I cant find ENDofMSG flag\n";
                out = false;
                break;
            }
            outmsg += "0 ERROR: UNKNOWN option [" + test + "]\n";
            out = false;
            break;
        }
        return out;
    }

    inline bool Scenario::parse_incoming_msg2(std::string &msg, std::string &outmsg)
    {
        clear();
        outmsg.clear();

        const char *p   = msg.c_str();
        const char *end = p + msg.size();

        struct Tok { const char *b; std::size_t n; };

        auto skip_ws = [&]() {
            while (p < end && std::isspace(static_cast<unsigned char>(*p))) ++p;
        };

        auto next_tok = [&]() -> Tok {
            skip_ws();
            const char *b = p;
            while (p < end && !std::isspace(static_cast<unsigned char>(*p))) ++p;
            return Tok{b, static_cast<std::size_t>(p - b)};
        };

        auto tok_empty = [&](const Tok &t) -> bool { return t.n == 0; };

        auto tok_eq = [&](const Tok &t, const char *lit) -> bool {
            const std::size_t ln = std::strlen(lit);
            return (t.n == ln) && (std::memcmp(t.b, lit, ln) == 0);
        };

        auto tok_to_string = [&](const Tok &t) -> std::string {
            return std::string(t.b, t.n);
        };

        // Null-terminate token for strtod/strtol without mutating msg.
        auto with_cstr = [&](const Tok &t, const char *&s, char buf[128], std::string &tmp) {
            if (t.n < 128) {
                std::memcpy(buf, t.b, t.n);
                buf[t.n] = '\0';
                s = buf;
            } else {
                tmp.assign(t.b, t.n);
                s = tmp.c_str();
            }
        };

        auto parse_int = [&](const Tok &t, int &v) -> bool {
            if (t.n == 0) return false;
            char buf[128]; std::string tmp; const char *s = nullptr;
            with_cstr(t, s, buf, tmp);
            char *ep = nullptr;
            long r = std::strtol(s, &ep, 10);
            if (ep == s || *ep != '\0') return false;
            v = static_cast<int>(r);
            return true;
        };

        auto parse_double = [&](const Tok &t, double &v) -> bool {
            if (t.n == 0) return false;
            char buf[128]; std::string tmp; const char *s = nullptr;
            with_cstr(t, s, buf, tmp);
            char *ep = nullptr;
            double r = std::strtod(s, &ep);
            if (ep == s || *ep != '\0') return false;
            v = r;
            return true;
        };

        auto parse_bool01 = [&](const Tok &t, bool &v) -> bool {
            int tmp = 0;
            if (!parse_int(t, tmp)) return false;
            v = (tmp != 0);
            return true;
        };

        auto require_tok = [&](Tok &t, const char *ctx_key) -> bool {
            t = next_tok();
            if (tok_empty(t)) {
                outmsg += "0 ERROR: Malformed value(s) after key [";
                outmsg += ctx_key;
                outmsg += "]\n";
                return false;
            }
            return true;
        };

        bool out = false;
        int counter = 0;

        while (true)
        {
            if (++counter > 500) {
                outmsg += "0 ERROR: After 200 iterations I cant find ENDofMSG flag\n";
                return false;
            }

            Tok key = next_tok();

            if (tok_empty(key)) {
                outmsg += "0 ERROR: Empty message was found\n";
                return false;
            }

            if (tok_eq(key, "ENDofMSG")) {
                out = true;
                break;
            }

            Tok v1, v2, v3, v4;

            // ---- string fields ----
            if (tok_eq(key, "modelArea")) {
                if (!require_tok(v1, "modelArea")) return false;
                modelArea = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "bMap")) {
                if (!require_tok(v1, "bMap")) return false;
                mapID = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "flowScen")) {
                if (!require_tok(v1, "flowScen")) return false;
                flowScen = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "wellType")) {
                if (!require_tok(v1, "wellType")) return false;
                wellType = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "unsatScen")) {
                if (!require_tok(v1, "unsatScen")) return false;
                unsatScenario = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "loadScen")) {
                if (!require_tok(v1, "loadScen")) return false;
                loadScen = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "loadSubScen")) {
                if (!require_tok(v1, "loadSubScen")) return false;
                loadSubScen = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "modifierName")) {
                if (!require_tok(v1, "modifierName")) return false;
                modifierName = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "DebugID")) {
                if (!require_tok(v1, "DebugID")) return false;
                debugID = tok_to_string(v1);
                // If your old behavior was to enable extra prints when debugID exists:
                // printAdditionalInfo = !debugID.empty();
                continue;
            }

            // ---- vector<string> ----
            if (tok_eq(key, "Nregions")) {
                if (!require_tok(v1, "Nregions")) return false;
                int n = 0;
                if (!parse_int(v1, n) || n < 0) return false;

                regionIDs.clear();
                regionIDs.reserve(static_cast<std::size_t>(n));

                for (int i = 0; i < n; ++i) {
                    Tok rid = next_tok();
                    if (tok_empty(rid)) return false;
                    regionIDs.push_back(tok_to_string(rid));
                }
                continue;
            }

            // ---- numeric fields ----
            if (tok_eq(key, "unsatWC")) {
                if (!require_tok(v1, "unsatWC") || !parse_double(v1, unsatZoneMobileWaterContent)) return false;
                continue;
            }

            if (tok_eq(key, "rchMap")) {
                bUseFlowRch = false;
                if (!require_tok(v1, "rchMap")) return false;
                rchName = tok_to_string(v1);
                continue;
            }

            if (tok_eq(key, "minRch")) {
                if (!require_tok(v1, "minRch") || !parse_double(v1, minRecharge)) return false;
                continue;
            }

            if (tok_eq(key, "startSimYear")) {
                if (!require_tok(v1, "startSimYear") || !parse_int(v1, startSimulationYear)) return false;
                continue;
            }

            if (tok_eq(key, "endSimYear")) {
                if (!require_tok(v1, "endSimYear") || !parse_int(v1, endSimulationYear)) return false;
                continue;
            }

            if (tok_eq(key, "startRed")) {
                if (!require_tok(v1, "startRed") || !parse_int(v1, startReductionYear)) return false;
                continue;
            }

            if (tok_eq(key, "endRed")) {
                if (!require_tok(v1, "endRed") || !parse_int(v1, endReductionYear)) return false;
                continue;
            }

            if (tok_eq(key, "maxConc")) {
                if (!require_tok(v1, "maxConc") || !parse_double(v1, maxConc)) return false;
                continue;
            }

            if (tok_eq(key, "surfConc")) {
                if (!require_tok(v1, "surfConc") || !parse_double(v1, surfConc)) return false;
                continue;
            }

            if (tok_eq(key, "surfRivDist")) {
                if (!require_tok(v1, "surfRivDist") || !parse_double(v1, surfRivDist)) return false;
                continue;
            }

            if (tok_eq(key, "surfRivInfl")) {
                if (!require_tok(v1, "surfRivInfl") || !parse_double(v1, surfRivInfl)) return false;
                continue;
            }

            if (tok_eq(key, "maxAge")) {
                if (!require_tok(v1, "maxAge") || !parse_double(v1, maxAge)) return false;
                continue;
            }

            if (tok_eq(key, "constRed")) {
                userSuppliedConstRed = true;
                if (!require_tok(v1, "constRed")) return false;
                double tmp = 0.0;
                if (!parse_double(v1, tmp)) return false;
                constReduction = static_cast<float>(tmp);
                continue;
            }

            if (tok_eq(key, "loadTrans")) {
                if (!require_tok(v1, "loadTrans")) return false;
                LoadTransitionName = tok_to_string(v1);

                if (!require_tok(v2, "loadTrans") || !parse_int(v2, LoadTransitionStart)) return false;
                if (!require_tok(v3, "loadTrans") || !parse_int(v3, LoadTransitionEnd)) return false;

                buseLoadTransition = (LoadTransitionName.compare("NONE") != 0);
                continue;
            }

            if (tok_eq(key, "por")) {
                if (!require_tok(v1, "por") || !parse_int(v1, porosity)) return false;
                continue;
            }

            if (tok_eq(key, "ADELR")) {
                if (!require_tok(v1, "ADELR") || !parse_double(v1, adeLambda)) return false;
                if (!require_tok(v2, "ADELR") || !parse_double(v2, adeR)) return false;
                continue;
            }

            if (tok_eq(key, "maxSourceCells")) {
                if (!require_tok(v1, "maxSourceCells") || !parse_int(v1, maxSourceCells)) return false;
                continue;
            }

            if (tok_eq(key, "PixelRadius")) {
                if (!require_tok(v1, "PixelRadius") || !parse_int(v1, PixelRadius)) return false;
                continue;
            }

            // ---- enums / converters ----
            if (tok_eq(key, "modifierType")) {
                if (!require_tok(v1, "modifierType")) return false;
                modifierType = string2RasterOperation(tok_to_string(v1));
                continue;
            }

            if (tok_eq(key, "modifierUnit")) {
                if (!require_tok(v1, "modifierUnit")) return false;
                modifierUnit = string2LoadUnits(tok_to_string(v1));
                continue;
            }

            if (tok_eq(key, "urfType")) {
                if (!require_tok(v1, "urfType")) return false;
                urfType = string2URFTYPE(tok_to_string(v1));
                continue;
            }

            // ---- Ncrops map<int,double> + globalReduction ----
            if (tok_eq(key, "Ncrops")) {
                if (!require_tok(v1, "Ncrops")) return false;
                int Ncrops = 0;
                if (!parse_int(v1, Ncrops) || Ncrops < 0) return false;

                for (int i = 0; i < Ncrops; ++i) {
                    Tok t_crop = next_tok();
                    Tok t_perc = next_tok();
                    if (tok_empty(t_crop) || tok_empty(t_perc)) return false;

                    int cropid = 0;
                    double perc = 0.0;
                    if (!parse_int(t_crop, cropid) || !parse_double(t_perc, perc)) return false;

                    if (cropid == -9) globalReduction = perc;
                    else              LoadReductionMap.insert(std::make_pair(cropid, perc));
                }
                continue;
            }

            // ---- selections / ranges ----
            if (tok_eq(key, "RadSelect")) {
            bNarrowSelection = true;
            useRadSelect = true;

            double cx=0.0, cy=0.0, r=0.0;
            if (!require_tok(v1, "RadSelect") || !parse_double(v1, cx)) return false;
            if (!require_tok(v2, "RadSelect") || !parse_double(v2, cy)) return false;
            if (!require_tok(v3, "RadSelect") || !parse_double(v3, r )) return false;

            RadSelect.setData(cx, cy, r);
            continue;
        }

        if (tok_eq(key, "RectSelect")) {
            bNarrowSelection = true;
            useRectSelect = true;

            double xmin=0.0, ymin=0.0, xmax=0.0, ymax=0.0;
            if (!require_tok(v1, "RectSelect") || !parse_double(v1, xmin)) return false;
            if (!require_tok(v2, "RectSelect") || !parse_double(v2, ymin)) return false;
            if (!require_tok(v3, "RectSelect") || !parse_double(v3, xmax)) return false;
            if (!require_tok(v4, "RectSelect") || !parse_double(v4, ymax)) return false;

            RectSelect.setData(xmin, ymin, xmax, ymax);
            continue;
        }

        if (tok_eq(key, "DepthRange")) {
            bNarrowSelection = true;
            double dmin=0.0, dmax=0.0;
            if (!require_tok(v1, "DepthRange") || !parse_double(v1, dmin)) return false;
            if (!require_tok(v2, "DepthRange") || !parse_double(v2, dmax)) return false;
            DepthRange.setData(dmin, dmax);
            continue;
        }

        if (tok_eq(key, "UnsatRange")) {
            bNarrowSelection = true;
            double umin=0.0, umax=0.0;
            if (!require_tok(v1, "UnsatRange") || !parse_double(v1, umin)) return false;
            if (!require_tok(v2, "UnsatRange") || !parse_double(v2, umax)) return false;
            unsatRange.setData(umin, umax);
            continue;
        }

        if (tok_eq(key, "Wt2tRange")) {
            bNarrowSelection = true;
            double wtmin=0.0, wtmax=0.0;
            if (!require_tok(v1, "Wt2tRange") || !parse_double(v1, wtmin)) return false;
            if (!require_tok(v2, "Wt2tRange") || !parse_double(v2, wtmax)) return false;
            wt2tRange.setData(wtmin, wtmax);
            continue;
        }

        if (tok_eq(key, "ScreenLenRange")) {
            bNarrowSelection = true;
            double slmin=0.0, slmax=0.0;
            if (!require_tok(v1, "ScreenLenRange") || !parse_double(v1, slmin)) return false;
            if (!require_tok(v2, "ScreenLenRange") || !parse_double(v2, slmax)) return false;
            ScreenLengthRange.setData(slmin, slmax);
            continue;
        }

            if (tok_eq(key, "SourceArea")) {
                int nPix=0, minPix=0, maxPix=0;
                double percPix=0.0;

                if (!require_tok(v1, "SourceArea") || !parse_int(v1, nPix)) return false;
                if (!require_tok(v2, "SourceArea") || !parse_int(v2, minPix)) return false;
                if (!require_tok(v3, "SourceArea") || !parse_int(v3, maxPix)) return false;
                if (!require_tok(v4, "SourceArea") || !parse_double(v4, percPix)) return false;

                SourceArea.setParameters(nPix, minPix, maxPix, percPix);
                continue;
            }

            // ---- flags ----
            if (tok_eq(key, "getids")) {
                if (!require_tok(v1, "getids")) return false;
                if (!parse_bool01(v1, printWellIds)) return false;
                continue;
            }

            if (tok_eq(key, "getotherinfo")) {
                if (!require_tok(v1, "getotherinfo")) return false;
                if (!parse_bool01(v1, printWellinfo)) return false;
                continue;
            }

            if (tok_eq(key, "printLF")) {
                if (!require_tok(v1, "printLF")) return false;
                if (!parse_bool01(v1, printLF)) return false;
                continue;
            }

            if (tok_eq(key, "printURF")) {
                if (!require_tok(v1, "printURF")) return false;
                if (!parse_bool01(v1, printURF)) return false;
                continue;
            }

            if (tok_eq(key, "printBTC")) {
                if (!require_tok(v1, "printBTC")) return false;
                if (!parse_bool01(v1, printBTC)) return false;
                continue;
            }

            if (tok_eq(key, "printWellBTC")) {
                if (!require_tok(v1, "printWellBTC")) return false;
                if (!parse_bool01(v1, printWellBTC)) return false;
                continue;
            }

            // ---- initConcParam ----
            if (tok_eq(key, "initConcParam")) {
                Tok t1 = next_tok(), t2 = next_tok(), t3 = next_tok(), t4 = next_tok();
                if (tok_empty(t1) || tok_empty(t2) || tok_empty(t3) || tok_empty(t4)) return false;

                if (!parse_double(t1, initcondparam.mean)) return false;
                if (!parse_double(t2, initcondparam.std )) return false;
                if (!parse_double(t3, initcondparam.minv)) return false;
                if (!parse_double(t4, initcondparam.maxv)) return false;

                if (initcondparam.mean < -90 || initcondparam.std < 0 ||
                    initcondparam.minv < 0  || initcondparam.maxv < 0) {
                    initcondparam.reset();
                    } else {
                        initcondparam.bUse = true;
                    }
                continue;
            }

            // ---- unknown key ----
            outmsg += "0 ERROR: UNKNOWN option [";
            outmsg.append(key.b, key.n);
            outmsg += "]\n";
            return false;
        }
        return out;


    }

}

#endif //MANTISSERVER_SCENARIO_H
