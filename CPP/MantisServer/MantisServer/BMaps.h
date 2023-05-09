//
// Created by giorg on 2/24/2023.
//

#ifndef MANTISSERVER_BMAPS_H
#define MANTISSERVER_BMAPS_H

namespace mantisServer{

    class GeoUnit{
    public:
        GeoUnit(){}
        void addWell(std::string flowScenarioName, int wellid);
        void getWells(std::string flowScenarioName, std::vector<int> &wellids);
        bool hasFlowWellMap(std::string &fswName);
    private:
        std::map<std::string, std::vector<int> > FlowScenWellsMap;
    };
    bool GeoUnit::hasFlowWellMap(std::string &fswName) {
        std::map<std::string, std::vector<int> >::iterator it;
        it = FlowScenWellsMap.find(fswName);
        if (it != FlowScenWellsMap.end()){
            return true;
        }
        else{
            return false;
        }
    }

    void GeoUnit::addWell(std::string flowScenarioName, int wellid) {
        std::map<std::string, std::vector<int> >::iterator it;
        it = FlowScenWellsMap.find(flowScenarioName);
        if (it == FlowScenWellsMap.end()){
            std::vector<int> tmp;
            tmp.push_back(wellid);
            FlowScenWellsMap.insert(std::pair<std::string, std::vector<int>>(flowScenarioName, tmp));
        }
        else{
            it->second.push_back(wellid);
        }
    }

    void GeoUnit::getWells(std::string flowScenarioName, std::vector<int> &wellids) {
        std::map<std::string, std::vector<int> >::iterator it;
        it = FlowScenWellsMap.find(flowScenarioName);
        if (it != FlowScenWellsMap.end()){
            wellids = it->second;
        }
    }

    class BMapLayer{
    public:
        BMapLayer(){}
        bool addwell(std::string geoUnitName, std::string flowScenarioName, int wellid);
        void getWells(std::vector<std::string> geoUnitNames,
                      std::string flowScenarioName, std::vector<int> &wellids);
        void addGeoUnit(std::string geoUnitName);
        bool hasGeoUnit(std::string geoUnitName);
        bool hasGeoUnitFlowScenWell(std::string &gunitName, std::string &fswName, std::string &outmsg);

    private:
        std::map<std::string, GeoUnit> geoUnits;
    };

    bool BMapLayer::hasGeoUnitFlowScenWell(std::string &gunitName, std::string &fswName, std::string &outmsg) {
        std::map<std::string, GeoUnit>::iterator it;
        it = geoUnits.find(gunitName);
        if (it != geoUnits.end()){
            if (it->second.hasFlowWellMap(fswName)){
                return true;
            }
            else{
                outmsg = "subregion [" + gunitName + "]";
                outmsg += "does not have the flow scenario wells [" + fswName + "]";
                return false;
            }
        }
        else{
            outmsg = "does not have the subregion [" + gunitName + "]";
            return false;
        }
    }

    bool BMapLayer::hasGeoUnit(std::string geoUnitName) {
        std::map<std::string, GeoUnit>::iterator it;
        it = geoUnits.find(geoUnitName);
        if (it == geoUnits.end())
            return false;
        else
            return true;
    }

    void BMapLayer::addGeoUnit(std::string geoUnitName) {
        std::map<std::string, GeoUnit>::iterator it;
        it = geoUnits.find(geoUnitName);
        if (it == geoUnits.end()){
            GeoUnit GU;
            geoUnits.insert(std::pair<std::string, GeoUnit>(geoUnitName, GU));
        }
    }

    bool BMapLayer::addwell(std::string geoUnitName, std::string flowScenarioName, int wellid){
        std::map<std::string, GeoUnit>::iterator it;
        it = geoUnits.find(geoUnitName);
        if (it != geoUnits.end()){
            it->second.addWell(flowScenarioName, wellid);
            return true;
        }
        else{
            std::cout << "GeoUnit with name [ " << geoUnitName << " ] not found in the Background Maps" << std::endl;
            return false;
        }

    }

    void BMapLayer::getWells(std::vector<std::string> geoUnitNames,
                             std::string flowScenarioName, std::vector<int> &wellids){
        std::map<std::string, GeoUnit>::iterator it;
        for (int i = 0; i < geoUnitNames.size(); ++i){
            it = geoUnits.find(geoUnitNames[0]);
            if (it != geoUnits.end()){
                std::vector<int> tmp;
                it->second.getWells(flowScenarioName,tmp);
                wellids.insert(wellids.end(), tmp.begin(), tmp.end());
            }
        }
    }


    class BMapCollection{
    public:
        BMapCollection(){}
        bool addwell(std::string GeoUnitName, std::string flowScenName, int wellid);
        void getWells(std::string LayerName, std::vector<std::string> GeoUnitNames,
                      std::string flowScenName, std::vector<int> &wellids);

        bool readData(std::string filename);
        int Nbmaps(){return BMaps.size();}
        bool hasMap(std::string &mapName);
        bool hasMapRegion(std::string &mapName, std::string &regName);
        bool validateMapsWells(std::string &mapName, std::vector<std::string> &regNames,
                               std::string &flowWellScen, std::string &outmsg);

    private:
        std::map<std::string, BMapLayer> BMaps;
        void addGeoUnit(std::string LayerName, std::string GeoUnitName);
    };

    bool BMapCollection::hasMap(std::string &mapName) {
        std::map<std::string, BMapLayer>::iterator it;
        it = BMaps.find(mapName);
        if (it != BMaps.end()){
            return true;
        }
        else{
            return false;
        }
    }

    bool BMapCollection::validateMapsWells(std::string &mapName, std::vector<std::string> &regNames,
                                          std::string &flowWellScen, std::string &outmsg) {
        std::map<std::string, BMapLayer>::iterator mapit;
        mapit = BMaps.find(mapName);
        if (mapit != BMaps.end()){
            std::vector<std::string>::iterator it;
            for (unsigned int i = 0; i < regNames.size(); ++i){
                bool tf = mapit->second.hasGeoUnitFlowScenWell(regNames[i], flowWellScen, outmsg);
                if (!tf){
                    outmsg = "Map [" + mapName + "] " + outmsg;
                    return false;
                }
            }
        }
        else{
            outmsg = "does not have the map [" + mapName + "]";
            return false;
        }
        return true;
    }

    bool BMapCollection::hasMapRegion(std::string &mapName, std::string &regName) {
        std::map<std::string, BMapLayer>::iterator it;
        it = BMaps.find(mapName);
        if (it != BMaps.end()){
            return it->second.hasGeoUnit(regName);
        }
        else{
            return false;
        }
    }

    void BMapCollection::addGeoUnit(std::string LayerName, std::string GeoUnitName) {
        std::map<std::string, BMapLayer>::iterator it;
        it = BMaps.find(LayerName);
        if (it == BMaps.end()){
            BMapLayer BML;
            BML.addGeoUnit(GeoUnitName);
            BMaps.insert(std::pair<std::string, BMapLayer>(LayerName, BML));
        }
        else{
            it->second.addGeoUnit(GeoUnitName);
        }
    }

    bool BMapCollection::addwell(std::string GeoUnitName,
                                 std::string flowScenName, int wellid){
        std::map<std::string, BMapLayer>::iterator it;
        bool tf = false;
        for (it = BMaps.begin(); it != BMaps.end(); ++it){
            if (it->second.hasGeoUnit(GeoUnitName)){
                tf = it->second.addwell(GeoUnitName,flowScenName,wellid);
                if (tf){
                    break;
                }
            }
        }
        return tf;
    }

    void BMapCollection::getWells(std::string LayerName, std::vector<std::string> GeoUnitNames,
                                  std::string flowScenName, std::vector<int> &wellids){
        std::map<std::string, BMapLayer>::iterator it;
        it = BMaps.find(LayerName);
        if (it != BMaps.end()){
            it->second.getWells(GeoUnitNames,flowScenName,wellids);
        }
    }

    bool BMapCollection::readData(std::string filename) {
        std::ifstream MAPSdatafile;
        MAPSdatafile.open(filename);
        if (!MAPSdatafile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            int Nmaps;
            std::string line;
            { // Get the number of background maps e.g. All area, Basins, counties, farms
                getline(MAPSdatafile, line);
                std::istringstream inp(line.c_str());
                inp >> Nmaps;
            }
            for (int imap = 0; imap < Nmaps; ++imap){
                std::string MapLayerName;
                int NGeoUnits;
                { // Get the layer name and the number of Subregions of this layer
                    std::getline(MAPSdatafile, line);
                    std::istringstream inp(line.c_str());
                    inp >> MapLayerName;
                    inp >> NGeoUnits;
                    for (int igu = 0; igu < NGeoUnits; ++igu){
                        std::string geoUnitName;
                        {// Get the name of the geoUnit
                            std::getline(MAPSdatafile, line);
                            std::istringstream inp(line.c_str());
                            inp >> geoUnitName;
                            addGeoUnit(MapLayerName, geoUnitName);
                        }
                    }
                }
            }
        }
        MAPSdatafile.close();
        return true;
    }
}

#endif //MANTISSERVER_BMAPS_H
