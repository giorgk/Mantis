//
// Created by giorg on 10/17/2024.
//

#ifndef MANTISSA_MS_URF_CALC_H
#define MANTISSA_MS_URF_CALC_H

#include<algorithm>

namespace MS{
    struct ADEoptions{
        double alpha = 0.32;
        double beta = 0.83;
        double R = 1.0;
        double lambda = 0.0;
        double dm = 1e-7;
    };
    const ADEoptions opt;
    const double sqrt2pi = std::sqrt(2*std::acos(-1.0));
    const double pi = std::atan(1)*4;
    std::vector<double> urfOne {0.48443252, 0.41642340, 0.08307405, 0.01338364, 0.00219913, 0.00039049, 0.00007584, 0.00001608, 0.00000370, 0.00000092, 0.00000023};
    std::vector<double> urfTwo {0.00095171, 0.19513243, 0.41957050, 0.25244126, 0.09333424, 0.02803643, 0.00771540, 0.00206063, 0.00055016, 0.00014916, 0.00004141, 0.00001182, 0.00000347, 0.00000106, 0.00000032};

    void calcURFs(int N, double m, double s, double a, double l, std::vector<double>& urf){
        double sumurf = 0.0;
        urf.clear();
        urf.resize(N,0.0);
        if (std::abs(m + 1) < 0.0001){
            for (unsigned int i = 0; i < std::min(static_cast<int>(urfOne.size()),N); ++i){
                urf[i] = urfOne[i];
                sumurf = sumurf + urf[i];
            }
        }
        else if (std::abs(m + 2) < 0.0001){
            for (unsigned int i = 0; i < std::min(static_cast<int>(urfTwo.size()),N); ++i){
                urf[i] = urfTwo[i];
                sumurf = sumurf + urf[i];
            }
        }
        else if (m < -65 || std::abs(m) < 0.0001){
            double x = l;
            double v = l/a;
            double aL = opt.alpha * std::pow(x, opt.beta);
            double D = aL*v + opt.dm;

            double t = 1;
            double c_prev = 0.0;
            double c_next = 0.0;
            for (int i = 0; i < N; ++i){
                double x_vt = x - v*t;
                double sqrDtx2 = 2*std::sqrt(D*t);
                c_next = 0.5*std::erfc(x_vt/sqrDtx2) +
                      std::sqrt(t*v*v/(pi*D)) *
                      exp(-1.0*(x_vt*x_vt)/(4*D*t)) -
                      0.5 * (1+ v*x/D + t*v*v/D) *
                      exp(v*x/D)*
                      std::erfc((x + v*t )/sqrDtx2);

                urf[i] = std::max(0.0, c_next - c_prev);
                sumurf = sumurf + urf[i];
                //std::cout << urf[i] << " ";
                t = t + 1.0;
                c_prev = c_next;
            }
            //std::cout << std::endl;

        }
        else{
            double t = 1;
            for (int i = 0; i < N; ++i){
                double p1 = (1 / (t*s*sqrt2pi));
                double p2 = (log(t) - m);
                p2 = -p2 * p2;
                p2 = p2 / (2 * s*s);
                urf[i] = p1*exp(p2);
                sumurf = sumurf + urf[i];

                t = t + 1.0;
                //std::cout << urf[i] << " ";
            }
            //std::cout << std::endl;
        }

        if (sumurf > 1.0){
            for (int i = 0; i < static_cast<int>(urf.size()); ++i){
                urf[i] = urf[i]/sumurf;
            }
        }

    }

    void convolute(std::vector<double>& urf, std::vector<double>& lf, std::vector<double>& btc, std::vector<double>& prebtc, double initconc){
        /*bool tmp = false;
        int iii = 0;
        iii++;
        if (tmp){
            for (int i = 0; i < urf.size(); ++i ){
                std::cout << urf[i] << " ";
            }
            std::cout << std::endl;
            for (int i = 0; i < lf_conc.size(); ++i ){
                std::cout << lf_conc[i] << " ";
            }
            std::cout << std::endl;
        }*/

        const double EPS = 1e-12;
        const std::size_t N = lf.size();

        double cumurf = 0.0;

        for (std::size_t i = 0; i < N; ++i){
            cumurf += urf[i];
            prebtc[i] = initconc * (1.0 - cumurf);
            if (std::abs(lf[i]) > EPS){
                const double lf_i = lf[i];
                for (std::size_t k = i; k < N; ++k){
                    btc[k] += urf[k - i] * lf_i;
                }
            }
        }


        /*unsigned int shift = 0;
        for (unsigned int i = 0; i < lf.size(); ++i){
            cumurf = cumurf + urf[i];
            prebtc[i] = initconc*(1.0 - cumurf);
            if (std::abs(lf[i]) > 0.00000000001){
                for (unsigned int k = shift; k < lf.size(); ++k){
                    btc[k] = btc[k] + urf[k - shift] * lf[i];
                }
            }
            shift = shift + 1;
        }*/
    }
}

#endif //MANTISSA_MS_URF_CALC_H
