#include <iostream>
#include <boost/mpi.hpp>


#include "MS_input.h"
#include "MS_raster.h"
#include "NPSAT_data.h"
#include "MS_mpi_utils.h"
#include "MS_unsat.h"

int main(int argc, char* argv[]) {
    boost::mpi::environment env( argc, argv );
    boost::mpi::communicator world;

    /*int N=30000000;
    std::vector<double> tmp(N,0);
    if (world.rank() == 0) {
        for (int i = 0; i < tmp.size(); ++i)
            tmp[i] = static_cast<double>(i)*2.36547;
    }
    //boost::mpi::broadcast(world, &tmp[0],N, 0);
    MS::sentFromRootToAllProc<double>(tmp, world);

    for (int k = 0; k < world.size(); ++k){
        for (int i = N-20; i < N; ++i){
            std::cout << tmp[i] << " ";
        }
        std::cout << std::endl;
    }
    return 0;*/


    MS::UserInput UI(world);

    UI.read(argc,argv);

    MS::UNSAT UZ;
    bool tf1 = UZ.readdata(UI.depth_input_file,UI.depth_name,UI.rch_input_file,UI.wc, UI.minDepth,UI.minRch);

    MS::BackgroundRaster backRaster;
    backRaster.readData(UI.rasteroptions.File,UI.rasteroptions.Nrows,UI.rasteroptions.Ncols,UI.rasteroptions.Ncells);

    MS::SWAT_data swat;
    bool tf = swat.read(UI.swat_input_file);



    //This is a map between Eid and cell that receiving water from this well
    // Only root processor knows this
    MS::WELL_CELLS well_cells;
    if (world.rank() == 0){
        tf = MS::readDistribPumping(UI.cell_well_file, well_cells);
    }
    world.barrier();

    MS::WELLS VI;
    MS::WELLS VD;
    tf = MS::readNPSATdata(UI.npsat_VI_file, VI, UI.porosity, UI.NsimYears, backRaster, world);
    tf = MS::readNPSATdata(UI.npsat_VD_file, VD, UI.porosity, UI.NsimYears, backRaster, world);

    std::vector<std::vector<int>> Eid_proc;
    sendWellEids(VI, Eid_proc, world);


    tf = MS::readInitSaltConc(UI.init_salt_VI_file, VI);
    tf = MS::readInitSaltConc(UI.init_salt_VD_file, VD);

    std::cout << "Proc: " << world.rank() << " has " << VI.size() << " VI wells and " << VD.size() << " VD wells" << std::endl;


    std::vector<double> ConcFromPump(backRaster.Ncell(), 0);

    // Main simulation loop
    MS::WELLS ::iterator itw;
    int hruidx;
    std::vector<double> totMfeed(UI.NsimYears, 0);
    for (int iyr = 0; iyr < UI.NsimYears; ++iyr){
        std::vector<double> wellConc;
        world.barrier();
        for (itw = VI.begin(); itw != VI.end(); ++itw){
            std::cout << itw->first << std::endl;
            std::vector<double> wellbtc(iyr + 1,0.0);
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                //std::cout << i << " " << std::flush;
                hruidx = itw->second.strml[i].hru_idx - 1;
                if (hruidx < 0){
                    continue;
                }

                // Calculate the groundwater concentration
                double gw_ratio = 0.0;
                if (std::abs(swat.irrGW_mm[iyr][hruidx] + swat.irrSW_mm[iyr][hruidx]) > 0.00000001){
                    gw_ratio = swat.irrGW_mm[iyr][hruidx] / (swat.irrGW_mm[iyr][hruidx] + swat.irrSW_mm[iyr][hruidx]);
                }

                // irrigated Salt Mass
                double m_gw = swat.irrsaltGW_Kgha[iyr][hruidx];

                //Calculate the percolated groundwater volume
                double v_gw = swat.perc_mm[iyr][hruidx] * gw_ratio;

                // Calculate concentration from NPSAT spreading (Feedback)
                double m_npsat = 0.0;
                if (itw->second.strml[i].IJ > 0){
                    m_npsat = ConcFromPump[itw->second.strml[i].IJ] * v_gw / 100.0;
                }

                double Mfeed = m_npsat - m_gw;
                if (Mfeed < 0) {
                    Mfeed = 0.0;
                }
                totMfeed[iyr] = totMfeed[iyr] + Mfeed;

                //Calculate the C_SWAT
                double c_swat = 0.0;
                if (swat.perc_mm[iyr][hruidx] > 0.00000001){
                    c_swat = (swat.totpercsalt_kgha[iyr][hruidx] + Mfeed) * 100 / swat.perc_mm[iyr][hruidx];
                }

                itw->second.strml[i].lf.push_back(c_swat);

                std::vector<double> btc(itw->second.strml[i].lf.size(), 0);
                std::vector<double> prebtc(itw->second.strml[i].lf.size(), 0);
                MS::convolute(itw->second.strml[i].urf, itw->second.strml[i].lf, btc, prebtc, itw->second.initConc);
                int shift = UZ.traveltime(itw->second.strml[i].IJ);

                for (unsigned int j = 0; j < btc.size(); ++j){
                    if (j <= shift){
                        wellbtc[j] = wellbtc[j] + itw->second.initConc;
                    }
                    else{
                        wellbtc[j] = wellbtc[j] + itw->second.strml[i].W * (btc[j-shift] + prebtc[j-shift]);
                    }
                }
            }// Loop streamlines
            std::cout << std::endl;
            wellConc.push_back(wellbtc.back());
        }// Loop wells

        // Go through the domestic wells but built only this year concentration
        for (itw = VD.begin(); itw != VD.end(); ++itw){
            std::cout << itw->first << std::endl;
            for (unsigned int i = 0; i < itw->second.strml.size(); ++i){
                //std::cout << i << " " << std::flush;
                hruidx = itw->second.strml[i].hru_idx - 1;
                if (hruidx < 0){
                    continue;
                }

                // Calculate the groundwater concentration
                double gw_ratio = 0.0;
                if (std::abs(swat.irrGW_mm[iyr][hruidx] + swat.irrSW_mm[iyr][hruidx]) > 0.00000001){
                    gw_ratio = swat.irrGW_mm[iyr][hruidx] / (swat.irrGW_mm[iyr][hruidx] + swat.irrSW_mm[iyr][hruidx]);
                }

                // irrigated Salt Mass
                double m_gw = swat.irrsaltGW_Kgha[iyr][hruidx];

                //Calculate the percolated groundwater volume
                double v_gw = swat.perc_mm[iyr][hruidx] * gw_ratio;

                // Calculate concentration from NPSAT spreading (Feedback)
                double m_npsat = 0.0;
                if (itw->second.strml[i].IJ > 0){
                    m_npsat = ConcFromPump[itw->second.strml[i].IJ] * v_gw / 100.0;
                }

                double Mfeed = m_npsat - m_gw;
                if (Mfeed < 0) {
                    Mfeed = 0.0;
                }

                //Calculate the C_SWAT
                double c_swat = 0.0;
                if (swat.perc_mm[iyr][hruidx] > 0.00000001){
                    c_swat = (swat.totpercsalt_kgha[iyr][hruidx] + Mfeed) * 100 / swat.perc_mm[iyr][hruidx];
                }

                itw->second.strml[i].lf.push_back(c_swat);
            }
        }

        // Send this year pumped concentration to root processor
        {
            std::vector<std::vector<double>> AllwellsConc(world.size());
            MS::sendVec2Root<double>(wellConc, AllwellsConc, world);
            // Put the well concentrations in the right cells
            if (world.rank() == 0){
                MS::WELL_CELLS::iterator it; // well_cells
                for (int i = 0; i < world.size(); ++i){
                    for (int j = 0; j < static_cast<int>(AllwellsConc[i].size()); ++j){
                        it = well_cells.find(Eid_proc[i][j]);
                        if (it != well_cells.end()){
                            for (int k = 0; k < it->second.size(); ++k){
                                ConcFromPump[it->second[k]] = AllwellsConc[i][j];
                            }
                        }
                    }
                }
            }

        }
        world.barrier();
        MS::sentFromRootToAllProc(ConcFromPump, world);

        /*{
            for (int k = 0; k < world.size(); ++k){
                for (int i = ConcFromPump.size()-20; i < ConcFromPump.size(); ++i){
                    std::cout << ConcFromPump[i] << " ";
                }
                std::cout << std::endl;
            }
        }
        return 0;*/
    }


    return 0;
}
