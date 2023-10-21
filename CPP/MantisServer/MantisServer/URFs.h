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
        double R = 1;
        double lambda = 0.0;
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
        URF(int Nyrs, double paramA, double paramB, URFTYPE type, ADEoptions ade_opt = ADEoptions());

        //! calc_conc returns the concentration for a given year. This is actuall used internaly by \ref calc_urf method.
        double calc_conc(double t, ADEoptions ade_opt = ADEoptions());

        /**
        \brief convolute executes the convolution of the urf with the loading function.
        \param LF is the loading function. The size of the loading function must be equal to \ref urf.
        \param BTC is the output of the convolution. The size of BTC must be equal to \LF.
        Since this function will be called million times the function doesn't perform and check or allocation for the BTC.

        To improve performance the function does not even loop for loading values less than \ref zeroLoading.
        */
        void convolute(std::vector<double> &LF, std::vector<double> &BTC);

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
        //! calc_urf expands the urf. This takes place during contruction.
        void calc_urf(ADEoptions ade_opt = ADEoptions());
        /**
         A hard coded threshold. Loading values lower than this threshold are treated as zero
         therefore the convolution can skip an entire loop associated with the particular loading value.
         */
        double zeroLoading = 0.00000000001;

        //! THis is the URF type
        URFTYPE type;

        std::vector<double> TwoYears{0.00095171, 0.19513243, 0.41957050, 0.25244126, 0.09333424, 0.02803643, 0.00771540, 0.00206063, 0.00055016, 0.00014916, 0.00004141, 0.00001182, 0.00000347, 0.00000106, 0.00000032};
        std::vector<double> OneYear{0.48443252, 0.41642340, 0.08307405, 0.01338364, 0.00219913, 0.00039049, 0.00007584, 0.00001608, 0.00000370, 0.00000092, 0.00000023};
    };

    URF::URF(int Nyrs, double paramA_in, double paramB_in, URFTYPE type_in, ADEoptions ade_opt)
            :
            paramA(paramA_in),
            paramB(paramB_in),
            type(type_in)
    {
        urf.resize(Nyrs, 0);
        if (std::abs(paramA + 1) < 0.0000001){
            for (unsigned int i = 0; i < OneYear.size(); ++i){
                if (i < urf.size())
                    urf[i] = OneYear[i];
            }
        }
        else if(std::abs(paramA + 2) < 0.0000001){
            for (unsigned int i = 0; i < TwoYears.size(); ++i){
                if (i < urf.size())
                    urf[i] = TwoYears[i];
            }
        }
        else{
            calc_urf(ade_opt);
            //printVector<double>(urf, "URF");
        }


    }

    double URF::calc_conc(double t, ADEoptions ade_opt) {
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
                double D = std::pow(ade_opt.alpha*paramA, ade_opt.beta);
                double DRt = 2*std::sqrt(D*ade_opt.R*t);
                if (ade_opt.lambda < 0.000001){
                    out = 0.5*std::erfc(( ade_opt.R * paramA - paramB * t )/DRt) +
                          std::sqrt(t*paramB*paramB/pi*D*ade_opt.R) *
                          exp(-1*(ade_opt.R*paramA - paramB*t)*(ade_opt.R*paramA - paramB*t)/(4*D*ade_opt.R*t)) -
                          0.5 * (1+ paramB*paramA/D + t*paramB*paramB/D*ade_opt.R ) *
                          exp(paramB*paramA/D)*
                          std::erfc((ade_opt.R*paramA + paramB*t )/DRt);
                }
                else{
                    std::cout << "Not implemented for lambda > 0" << std::endl;
                }
            }
                break;
            default:
                out = std::numeric_limits<double>::quiet_NaN();
                break;
        }
        return out;
    }

    void URF::calc_urf(ADEoptions ade_opt) {
        double sumurf = 0.0;
        for (int i = 0; i < static_cast<int>(urf.size()); ++i){
            urf[i] = calc_conc(static_cast<double>(i + 1), ade_opt);
            sumurf += urf[i];
        }
        if (sumurf > 1.0){
            for (int i = 0; i < static_cast<int>(urf.size()); ++i){
                urf[i] = urf[i]/sumurf;
            }
        }
    }

    void URF::print_urf(std::ofstream& urf_file) {
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
}

#endif //MANTISSERVER_URFS_H
