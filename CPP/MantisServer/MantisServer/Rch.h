//
// Created by giorg on 2/24/2023.
//

#ifndef MANTISSERVER_RCH_H
#define MANTISSERVER_RCH_H

#include "MShelper.h"


namespace mantisServer{
    /**
 * RechargeScenario holds the data for one recharge scenario
 */
    class RechargeScenario{
    public:
        //! A usefull do nothing constructor
        RechargeScenario(){}

        /**
         * Reads the data for one recharge scenario.
         * The recharge scenario consists of two values.
         * The total recharge amount and the percentage of clean recharge
         * @param filename
         * @return
         */
        bool readData(std::string filename, int Nr);

        bool getValues(int lin_idx, double &rch, double &clean_prc);

    private:
        LinearData Recharge;
    };

    bool RechargeScenario::readData(std::string filename, int Nr) {
        Recharge.setNoDataValue(0.0);
        return Recharge.readData(filename, Nr);

    }

    bool RechargeScenario::getValues(int lin_idx, double &rch, double &clean_prc) {
        rch = Recharge.getValue(0,lin_idx);
        clean_prc = Recharge.getValue(1,lin_idx);
        if (rch == Recharge.getNoDataValue()){
            return false;
        }
        return true;
    }

    class RechargeScenarioList{
    public:
        RechargeScenarioList(){}

        bool readData(std::string path, std::string filename, int Nr);

        bool getValue(std::string scenarioName, int lin_idx, double &rch, double &clean_prc);
        bool hasRechargeMaps(std::string &rchmap);

        std::map<std::string, RechargeScenario> RechargeList;
    };

    bool RechargeScenarioList::hasRechargeMaps(std::string &rchmap) {
        std::map<std::string, RechargeScenario>::iterator it;
        it = RechargeList.find(rchmap);
        return it != RechargeList.end();
    }

    bool RechargeScenarioList::readData(std::string path, std::string filename, int Nr) {
        auto start = std::chrono::high_resolution_clock::now();

        std::ifstream rchMainFile;

        if (!path.empty()){
            filename = path + filename;
        }
        std::cout << "Reading " << filename << std::endl;
        rchMainFile.open(filename);
        if (!rchMainFile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            while (getline(rchMainFile, line)){
                std::string scenarioName, scenarioFile;
                std::istringstream inp(line.c_str());

                inp >> scenarioName;
                if (scenarioName.empty())
                    continue;
                if (scenarioName.front() == '#')
                    continue;

                inp >> scenarioFile;

                if (!path.empty()){
                    scenarioFile = path + scenarioFile;
                }
                RechargeScenario rch;
                bool tf = rch.readData(scenarioFile,Nr);
                if (tf){
                    RechargeList.insert(std::pair<std::string, RechargeScenario>(scenarioName, rch));
                }
                else{
                    return false;
                }
            }
        }
        rchMainFile.close();
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Recharge data in " << elapsed.count() << std::endl;
        return true;
    }

    bool RechargeScenarioList::getValue(std::string scenarioName, int lin_idx, double &rch, double &clean_prc) {
        std::map<std::string, RechargeScenario>::iterator it;
        it = RechargeList.find(scenarioName);
        if (it == RechargeList.end()){
            return false;
        }
        else{
            return it->second.getValues(lin_idx,rch,clean_prc);
        }
    }

}

#endif //MANTISSERVER_RCH_H
