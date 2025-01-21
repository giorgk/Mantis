//
// Created by giorgk on 5/9/2023.
//
#include <iomanip>

#include "MShelper.h"

#ifndef MANTISSERVER_URFS_H
#define MANTISSERVER_URFS_H

namespace mantisServer{

    /*!
	 * This is a container for ADE options. In future this will be deleted as we use only LGNRM type of URF
	 */
    struct ADEoptions{
        double alpha = 0.32;
        double beta = 0.83;
        double R = 1.0;
        double lambda = 0.0;
        double dm = 1e-7;
    };

    /*! \class URF
		\brief The URF class contains functionality related to the Unit Response Functions.

		At the moment we use only the LGNRM type of URF which is a function of the following form: \n

		\htmlonly
		<a href="https://www.codecogs.com/eqnedit.php?latex=\huge&space;URF(t)&space;=&space;\frac{1}{t&space;\cdot&space;s&space;\sqrt{2\pi}}\cdot&space;e^{-\frac{(\ln{(t)}-m)^2}{2\cdot&space;s^2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\huge&space;URF(t)&space;=&space;\frac{1}{t&space;\cdot&space;s&space;\sqrt{2\pi}}\cdot&space;e^{-\frac{(\ln{(t)}-m)^2}{2\cdot&space;s^2}}" title="\huge URF(t) = \frac{1}{t \cdot s \sqrt{2\pi}}\cdot e^{-\frac{(\ln{(t)}-m)^2}{2\cdot s^2}}" /></a>
		\endhtmlonly

		\n
		where \a m is the mean and \a s is the standard deviation. These parameters calculated by fitting the numerically generated URF from the
		<a href="https://gwt.ucdavis.edu/research-tools-and-applications/npsat-engine">NPSAT engine</a>. Here we store only the two values.
		This has two advantages. First for each urf we store just 2 double values instead ~200.
		Secondly we can adjust the simulation time at runtime as we can generate urfs for as longer periods.
		However extending the simulation to let's say 1000 years may not be appropriate.

		The units of the time parameter \a are in years. URF(20) is the value after 20 years of simulation.
	*/
    class URF {
    public:
        URF(){}
        /**
         * @brief Construct a new URF object
         *
         * The constructor expands the Unit respons Function using the input parameters for Nyrs.
         *
         * @param Nyrs is the number of years to expand the URF
         * @param paramA This is the first parameter.
         * @param paramB  This is the second parameter.
         * @param type If the type is URFTYPE::LGNRM paramA is the mean and paramB is the standard deviation
         * If the type is URFTYPE::ADE, paramA is the screen length and paramB is the velocity.
         */
        void init(int Nyrs, double paramA_in, double paramB_in, URFTYPE type_in, ADEoptions ade_opt_in = ADEoptions());

        //! calc_conc returns the concentration for a given year. This is actuall used internaly by \ref calc_urf method.
        double calc_conc(double t);

        /**
        \brief convolute executes the convolution of the urf with the loading function.
        \param LF is the loading function. The size of the loading function must be equal to \ref urf.
        \param BTC is the output of the convolution. The size of BTC must be equal to \LF.
        Since this function will be called million times the function doesn't perform and check or allocation for the BTC.

        To improve performance the function does not even loop for loading values less than \ref zeroLoading.
        */
        void convolute(std::vector<double> &LF, std::vector<double> &BTC);

        void convolute(std::vector<double> &LF, std::vector<double> &BTC, std::vector<double> &preBTC, double initconc);

        /**
         *
         * @param urf_file prints the full expanded urf into file
         */
        void print_urf(std::ofstream& urf_file);

    private:
        //! The mean fitted value.
        double paramA;
        //! The standard deviation fitted value
        double paramB;
        //! This is a vector where the expanded values of the urfs will be stored.
        std::vector<double> urf;
        std::vector<double> invert_cumurf;
        //! calc_urf expands the urf. This takes place during contruction.
        void calc_urf();
        /**
         A hard coded threshold. Loading values lower than this threshold are treated as zero
         therefore the convolution can skip an entire loop associated with the particular loading value.
         */
        double zeroLoading = 0.00000000001;

        //! THis is the URF type
        URFTYPE type;
        ADEoptions ade_opt;

        std::vector<double> TwoYears{0.00095171, 0.19513243, 0.41957050, 0.25244126, 0.09333424, 0.02803643, 0.00771540, 0.00206063, 0.00055016, 0.00014916, 0.00004141, 0.00001182, 0.00000347, 0.00000106, 0.00000032};
        std::vector<double> OneYear{0.48443252, 0.41642340, 0.08307405, 0.01338364, 0.00219913, 0.00039049, 0.00007584, 0.00001608, 0.00000370, 0.00000092, 0.00000023};
    };

    void URF::init(int Nyrs, double paramA_in, double paramB_in, URFTYPE type_in, ADEoptions ade_opt_in)
    {
        paramA = paramA_in; // mean or Length
        paramB = paramB_in; // std or age
        type = type_in;
        ade_opt = ade_opt_in;

        urf.resize(Nyrs, 0);
        invert_cumurf.resize(Nyrs,1);

        if (type == URFTYPE::LGNRM){
            if (std::abs(paramA) < 0.000000001){
                return;
            }
            if (std::abs(paramA + 1) < 0.0000001){
                double cumurf = 0;
                for (unsigned int i = 0; i < OneYear.size(); ++i){
                    if (i < urf.size()){
                        urf[i] = OneYear[i];
                        cumurf = cumurf + urf[i];
                        invert_cumurf[i] = 1.0 - cumurf;
                    }
                }
                return;
            }
            else if(std::abs(paramA + 2) < 0.0000001){
                double cumurf = 0;
                for (unsigned int i = 0; i < TwoYears.size(); ++i){
                    if (i < urf.size()){
                        urf[i] = TwoYears[i];
                        cumurf = cumurf + urf[i];
                        invert_cumurf[i] = 1.0 - cumurf;
                    }
                }
                return;
            }
        }
        else if (type == URFTYPE::ADE & std::abs(paramB) < 0.000000001){
            return;
        }

        calc_urf();
    }

    double URF::calc_conc(double t) {
        //paramA mean or Length
        //paramB std or age
        double out = 0.0;
        switch (type)
        {
            case URFTYPE::LGNRM:
            {
                //(1 / (x*b*sqrt(2 * pi)))*exp((-(log(x) - a) ^ 2) / (2 * b ^ 2))
                //std::cout << "t: (m,s)" << t << " (" << m << "," << s << ")" << std::endl;
                double p1 = (1 / (t*paramB*sqrt2pi));
                //std::cout << "p1: " << p1 << std::endl;
                //std::cout << "log(t): " << log(t) << std::endl;
                double p2 = (log(t) - paramA);
                //std::cout << "p2: " << p2 << std::endl;
                p2 = -p2 * p2;
                //std::cout << "p2: " << p2 << std::endl;
                p2 = p2 / (2 * paramB*paramB);
                //std::cout << "p2: " << p2 << std::endl;
                //std::cout << "return: " << p1 * exp(p2) << std::endl;
                out = p1*exp(p2);
            }
                break;
            case URFTYPE::ADE:
            {
                double x = paramA;
                double v = paramA/paramB;
                double aL = ade_opt.alpha * std::pow(x,ade_opt.beta);
                double D = aL*v + ade_opt.dm;

                if (std::abs(ade_opt.lambda) < 0.000000001 & std::abs(ade_opt.R - 1.0) < 0.0001){
                    // No decay, no retardation
                    double x_vt = x - v*t;
                    double sqrDtx2 = 2*std::sqrt(D*t);
                    out = 0.5*std::erfc(x_vt/sqrDtx2) +
                          std::sqrt(t*v*v/(pi*D)) *
                          exp(-1.0*(x_vt*x_vt)/(4*D*t)) -
                          0.5 * (1+ v*x/D + t*v*v/D) *
                          exp(v*x/D)*
                          std::erfc((x + v*t )/sqrDtx2);
                }
                else if (std::abs(ade_opt.lambda) < 0.000000001){
                    //No decay
                    double Rx_vt = ade_opt.R*x - v*t;
                    double sqrDRtx2 = 2*std::sqrt(D*ade_opt.R*t);

                    out = 0.5*std::erfc(Rx_vt/sqrDRtx2) +
                          std::sqrt(t*v*v/(pi*D*ade_opt.R)) *
                          exp(-1.0*(Rx_vt*Rx_vt)/(4*D*ade_opt.R*t)) -
                          0.5 * (1+ v*x/D + (t*v*v)/(D*ade_opt.R)) *
                          exp(v*x/D)*
                          std::erfc((ade_opt.R*x + v*t )/sqrDRtx2);
                }
                else{
                    double U = sqrt(v*v+4*D*ade_opt.R*ade_opt.lambda);
                    double sqrDRtx2 = 2*std::sqrt(D*ade_opt.R*t);
                    double vmU = v - U;
                    double vpU = v + U;
                    double Rx = ade_opt.R * x;
                    double Ut = U*t;

                    out = v/vpU * exp(x*vmU/(2*D)) * std::erfc((Rx - Ut)/sqrDRtx2) +
                            v/vmU * exp(x*vpU/(2*D)) * std::erfc((Rx + Ut)/sqrDRtx2) +
                            (v*v)/(2*D*ade_opt.R*ade_opt.lambda) *
                            exp(v*x/D - ade_opt.lambda*t) * std::erfc((Rx + v*t)/sqrDRtx2);
                }
            }
                break;
            default:
                out = std::numeric_limits<double>::quiet_NaN();
                break;
        }
        return out;
    }

    void URF::calc_urf() {
        double sumurf = 0.0;
        double v_prev = 0.0;
        double v_curr = 0.0;
        for (int i = 0; i < static_cast<int>(urf.size()); ++i){
            if (type == URFTYPE::LGNRM){
                urf[i] = calc_conc(static_cast<double>(i + 1));
            }
            else if (type == URFTYPE::ADE){
                if (i == 0){
                    v_prev = calc_conc(static_cast<double>(i + 1));
                    urf[i] = v_prev;
                }
                else{
                    v_curr = calc_conc(static_cast<double>(i + 1));
                    urf[i] = v_curr - v_prev;
                    v_prev = v_curr;
                }
            }
            sumurf += urf[i];
            invert_cumurf[i] = 1.0 - sumurf;
            if (sumurf > 1.0){
                invert_cumurf[i] = 0.0;
            }
        }
        if (sumurf > 1.0){
            for (int i = 0; i < static_cast<int>(urf.size()); ++i){
                urf[i] = urf[i]/sumurf;
            }
        }
    }

    void URF::print_urf(std::ofstream& urf_file) {

        urf_file << std::scientific<< std::setprecision(10) << paramA << " " << paramB << " ";

        for (int i = 0; i < static_cast<int>(urf.size()); ++i)
            urf_file << std::scientific << std::setprecision(10) << urf[i] << " ";
        urf_file << std::endl;
    }

    void URF::convolute(std::vector<double> &LF, std::vector<double> &BTC) {
        //std::cout << "LF size: " << LF.size() << std::endl;
        //BTC.resize(LF.size(), 0.0);
        //std::cout << "BTC size: " << BTC.size() << std::endl;
        int shift = 0;
        for (int i = 0; i < static_cast<int>(LF.size()); ++i) {
            if (std::abs(LF[i] - 0) > zeroLoading) {
                for (int k = shift; k < static_cast<int>(LF.size()); ++k) {
                    //std::cout << k - shift << " : " << urf[k - shift] << " " << LF[i] << std::endl;
                    BTC[k] = BTC[k] + urf[k - shift] * LF[i];
                }
            }
            shift = shift + 1;
        }
    }

    void URF::convolute(std::vector<double> &LF, std::vector<double> &BTC, std::vector<double> &preBTC,
                        double initconc) {
        double cumurf = 0.0;
        unsigned int shift = 0;
        for (int i = 0; i < static_cast<int>(LF.size()); ++i){
            cumurf = cumurf + urf[i];
            preBTC[i] = initconc * (1.0 - cumurf);
            if (std::abs(LF[i] - 0) > zeroLoading){
                for (int k = shift; k < static_cast<int>(LF.size()); ++k){
                    BTC[k] = BTC[k] + urf[k - shift] * LF[i];
                }
            }
            shift = shift + 1;
        }
    }
}

#endif //MANTISSERVER_URFS_H
