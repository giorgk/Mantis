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
    inline const ADEoptions opt{};
    inline const double sqrt2pi = std::sqrt(2.0 * std::acos(-1.0));
    inline const double pi = 4.0 * std::atan(1.0);
    inline const std::vector<double> urfOne {
        0.48443252, 0.41642340, 0.08307405, 0.01338364, 0.00219913,
        0.00039049, 0.00007584, 0.00001608, 0.00000370, 0.00000092, 0.00000023
    };

    inline const std::vector<double> urfTwo {
        0.00095171, 0.19513243, 0.41957050, 0.25244126, 0.09333424,
        0.02803643, 0.00771540, 0.00206063, 0.00055016, 0.00014916,
        0.00004141, 0.00001182, 0.00000347, 0.00000106, 0.00000032
    };

    inline void calcURFs(int N, double m, double s, double a, double l, std::vector<double>& urf){
        urf.clear();

        if (N <= 0) {
            return;
        }

        const std::size_t NN = static_cast<std::size_t>(N);
        urf.assign(NN, 0.0);

        if (!std::isfinite(m) || !std::isfinite(s) || !std::isfinite(a) || !std::isfinite(l)) {
            return;
        }

        double sumurf = 0.0;

        if (std::abs(m + 1.0) < 1e-4) {
            const std::size_t ncopy = std::min<std::size_t>(urfOne.size(), NN);
            for (std::size_t i = 0; i < ncopy; ++i) {
                urf[i] = urfOne[i];
                sumurf += urf[i];
            }
        }
        else if (std::abs(m + 2.0) < 1e-4) {
            const std::size_t ncopy = std::min<std::size_t>(urfTwo.size(), NN);
            for (std::size_t i = 0; i < ncopy; ++i) {
                urf[i] = urfTwo[i];
                sumurf += urf[i];
            }
        }
        else if (m < -65.0 || std::abs(m) < 1e-4) {
            // ADE / Ogata-Banks-like cumulative response, converted to discrete increments
            if (a <= 0.0 || l <= 0.0) {
                return;
            }
            const double x = l;
            const double v = l/a;
            const double aL = opt.alpha * std::pow(x, opt.beta);
            const double D = aL*v + opt.dm;

            double t = 1;
            double c_prev = 0.0;

            for (std::size_t i = 0; i < NN; ++i) {
                const double sqrDtx2 = 2.0 * std::sqrt(D * t);
                if (!(sqrDtx2 > 0.0) || !std::isfinite(sqrDtx2)) {
                    return;
                }

                const double x_vt = x - v * t;

                const double term1 = 0.5 * std::erfc(x_vt / sqrDtx2);

                const double term2 = std::sqrt(t * v * v / (pi * D)) * std::exp(-(x_vt * x_vt) / (4.0 * D * t));

                const double exp_arg = v * x / D;
                const double term3_factor = 0.5 * (1.0 + v * x / D + t * v * v / D);

                const double term3 = term3_factor * std::exp(exp_arg) * std::erfc((x + v * t) / sqrDtx2);

                const double c_next = term1 + term2 - term3;

                if (!std::isfinite(c_next)) {
                    return;
                }

                urf[i] = std::max(0.0, c_next - c_prev);
                sumurf += urf[i];

                c_prev = c_next;
                t += 1.0;
            }
        }
        else{
            // Lognormal PDF sampled at integer times t = 1,2,...,N
            if (s <= 0.0) {
                return;
            }
            double t = 1.0;
            const double s2 = s * s;
            const double denom_const = s * sqrt2pi;

            for (std::size_t i = 0; i < NN; ++i) {
                const double lt = std::log(t);
                const double z = lt - m;
                const double expo = -(z * z) / (2.0 * s2);
                const double val = (1.0 / (t * denom_const)) * std::exp(expo);

                if (!std::isfinite(val)) {
                    return;
                }

                urf[i] = val;
                sumurf += val;
                t += 1.0;
            }
            //std::cout << std::endl;
        }

        if (sumurf > 1.0  && std::isfinite(sumurf)){
            for (std::size_t i = 0; i < NN; ++i) {
                urf[i] /= sumurf;
            }
        }
    }

    inline void convolute(std::vector<double>& urf, std::vector<double>& lf, std::vector<double>& btc, std::vector<double>& prebtc, double initconc){
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
