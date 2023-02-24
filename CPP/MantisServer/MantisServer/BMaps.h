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
    private:
        std::map<std::string, std::vector<int> > FlowScenWellsMap;
    };

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
        wellids.clear();
        std::map<std::string, std::vector<int> >::iterator it;
        it = FlowScenWellsMap.find(flowScenarioName);
        if (it != FlowScenWellsMap.end()){
            wellids = it->second;
        }
    }

    class BMapLayer{
    public:
        BMapLayer(){}
        void addwell(std::string geoUnitName, std::string flowScenarioName, int wellid);
        void getWells(std::vector<std::string> geoUnitNames, std::string flowScenarioName, std::vector<int> &wellids);

    private:
        std::map<std::string, GeoUnit> geoUnits;
    };


}

#endif //MANTISSERVER_BMAPS_H
