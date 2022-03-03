//
// Created by giorg on 2/28/2022.
//

#ifndef MANTISSERVER_RUNTIMEWELLS_H
#define MANTISSERVER_RUNTIMEWELLS_H
#include "string"

namespace mantisServer{
    struct URFParam{
        int urfI;
        int urfJ;
        int inRiver;
        double m;
        double s;
        double w;
        double depth;
    };
    struct Entity{
        std::vector<URFParam> S;
        double X;
        double Y;
    };

    class runtimeURFSet{
    public:
        runtimeURFSet(){}
        bool readWellSet(std::string &nameSet, std::vector<std::string> &fileList);
        void reset();
        int getNurfs(){return Nurfs;};
        int getNwells(){return static_cast<int>(urfSet.size());};
        bool checkNameSet(std::string nm){return nameSet.compare(nm) == 0;}
        int getNsid(unsigned int id);
        bool getParam(unsigned int eid, unsigned int sid,
                      double &m, double &s, double &w, double &d, int &I, int &J, int &riv);
        bool getWellCoords(int iw, double &x, double &y);
        int getlife(){return Nlife;}
        void reduceLife(){Nlife--;}
        void increaseLife(){Nlife++;}
        void setlife(int l){Nlife = l;}

    private:
        std::string nameSet;
        std::vector<std::string> fileList;
        std::vector<Entity> urfSet;
        bool readSet(std::string filename);
        int Nurfs = 0;
        int Nlife = 0;
    };

    void runtimeURFSet::reset() {
        fileList.clear();
        nameSet.clear();
        urfSet.clear();
        Nurfs = 0;
    }

    int runtimeURFSet::getNsid(unsigned int id) {
        if (id < urfSet.size()){
            return urfSet[id].S.size();
        }
        else{
            return 0;
        }
    }

    bool runtimeURFSet::getWellCoords(int iw, double &x, double &y) {
        if (iw < urfSet.size() && iw >= 0){
            x = urfSet[iw].X;
            y = urfSet[iw].Y;
            return true;
        }
        return false;
    }

    bool runtimeURFSet::getParam(unsigned int eid, unsigned int sid,
                                 double &m, double &s, double &w, double &d,
                                 int &I, int &J, int &riv) {
        if (eid < urfSet.size()){
            if (sid < urfSet[eid].S.size()){
                m = urfSet[eid].S[sid].m;
                s = urfSet[eid].S[sid].s;
                w = urfSet[eid].S[sid].w;
                d = urfSet[eid].S[sid].depth;
                I = urfSet[eid].S[sid].urfI;
                J = urfSet[eid].S[sid].urfJ;
                riv = urfSet[eid].S[sid].inRiver;
                return true;
            }
            else{
                return false;
            }
        }
        else{
            return false;
        }
    }

    bool runtimeURFSet::readWellSet(std::string &nameSet_in, std::vector<std::string> &fileList_in){
        if (nameSet.compare(nameSet_in) == 0){
            return true;
        }
        else{
            try{
                for (unsigned int i = 0; i < fileList_in.size(); ++i) {
                    bool tf = readSet(fileList_in[i]);
                    if (!tf) {
                        return false;
                    }
                    fileList.push_back(fileList_in[i]);
                }
                nameSet = nameSet_in;
                return true;
            }
            catch(...) {
                return false;
            }
        }
        return false;
    }

    bool runtimeURFSet::readSet(std::string filename) {
#if _USEHF>0
        std::string file_ext = filename + ".h5";
        const std::string EIDNameSet("EID");
        const std::string EXYNameSet("EXY");
        const std::string SEIJRNameSet("SEIJR");
        const std::string SMSWDNameSet("SMSWD");
        HighFive::File HDFfile(file_ext, HighFive::File::ReadOnly);
        HighFive::DataSet datasetEID = HDFfile.getDataSet(EIDNameSet);
        HighFive::DataSet datasetEXY = HDFfile.getDataSet(EXYNameSet);
        HighFive::DataSet datasetSEIJR = HDFfile.getDataSet(SEIJRNameSet);
        HighFive::DataSet datasetSMSWD = HDFfile.getDataSet(SMSWDNameSet);
        std::vector<int> eid;
        std::vector<std::vector<double>> exy;
        std::vector<std::vector<int>> seijr;
        std::vector<std::vector<double>> smswd;
        datasetEID.read(eid);
        datasetEXY.read(exy);
        datasetSEIJR.read(seijr);
        datasetSMSWD.read(smswd);
        Entity e;

        URFParam prm;
        int prevEid = -9999;
        int idx_well = 0;
        for (unsigned int i = 0; i < seijr[0].size(); ++i){
            prm.m = smswd[0][i];
            prm.s = smswd[1][i];
            prm.w = smswd[2][i];
            prm.depth = smswd[3][i];
            prm.urfI = seijr[1][i];
            prm.urfJ = seijr[2][i];
            prm.inRiver = seijr[3][i];
            if (prevEid != seijr[0][i]){
                if (i > 0){
                    urfSet.push_back(e);
                    e.S.clear();
                }
                e.X = exy[0][idx_well];
                e.Y = exy[1][idx_well];
                e.S.push_back(prm);
                Nurfs++;
                prevEid = seijr[0][i];
                idx_well++;
            }
            else if (i == seijr[0].size()-1){
                e.S.push_back(prm);
                Nurfs++;
                urfSet.push_back(e);
            }
            else{
                e.S.push_back(prm);
                Nurfs++;
            }
        }
        return true;
#endif
        return false;
    }
}


#endif //MANTISSERVER_RUNTIMEWELLS_H
