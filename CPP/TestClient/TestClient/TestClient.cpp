// TestClient.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#ifdef _WIN64
#include "pch.h"
//#endif

#include <iostream>
#include <string>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/algorithm/string.hpp>
//#include <boost/thread.hpp>
#include <math.h>
#include <cstdlib>

std::vector<std::string> split(const std::string& str, const std::string& delim)
{
	std::vector<std::string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == std::string::npos) pos = str.length();
		std::string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}

int main(int argc, char* argv[])
{
	std::cout << "Running version 2.0.04" << std::endl;
	// Test random numbers
	//for (int i = 0; i < 10; ++i)
	//	std::cout << static_cast<double>(std::rand() % 100)/10 << std::endl;
	//return 0;
	//std::vector<double> reductionCoeff(23, 1);
	//for (int i = 1945; i < 2100; i++) {
	//	//std::cout << (i-1940) % 25 << std::endl;
	//	std::cout << i << " -> " << std::floor((i-1945)/15) << std::endl;
	//}

	bool quit = false;
	bool ping = false;
	std::string infile, outfile;
	if (argc > 1) {
		infile = argv[1];
		if (infile.compare("quit") == 0) {
            quit = true;
        }
		else if (infile.compare("ping") == 0){
		    ping = true;
		}
	}

	if (argc > 2) {
		outfile = argv[2];
	}
	else {
		outfile = "testClientResults.dat";
	}

	std::string msg;
	if (quit) {
		msg = "quit\n";
	}
	else if (ping){
        msg = "ping\n";
	}
	else {
		if (argc > 1) {
			std::ifstream indata(infile.c_str());
			if (!indata.good()) {
				std::cout << "Can't open the file " << infile << std::endl;
			}
			else {
				std::string line, line1;
				bool firstLine = true;
				while (getline(indata, line)) {
				    if (line.empty())
				        continue;
				    if (line.front() == '#')
				        continue;

				    std::vector<std::string> words;
				    boost::split(words, line, boost::is_any_of(" "));

				    for (unsigned int i = 0; i < words.size(); ++i ){
				        boost::trim(words[i]);
				        if (words[i].empty())
                            continue;
                        if (firstLine){
                            msg += words[i];
                            firstLine = false;
                        }
                        else{
                            msg += " " + words[i];
                        }
				    }
				}
 				msg += " ENDofMSG\n";
			}
		}
		else {
			// Default message
			msg = "endSimYear 2100";
            msg += " startRed 2025";
            msg += " endRed 2035";
            msg += " flowScen C2VSIM_MAR24_Padj";
            msg += " wellType VI";
            msg += " loadScen TEST01";
            msg += " loadSubScen test100";
            msg += " unsatScen C2VSIM_MAR24";
            msg += " unsatWC 0.1";
            msg += " bMap Subregions";
            msg += " Nregions 1 Subregion11";
            msg += " por 20";
            msg += " urfType LGNRM";
            //msg += " ADELR 0 1";
            //msg += " isLoadConc 1";
            //msg += " RadSelect -91791 -33506 10000"; // X Y radius
            //msg += " RectSelect -81323.0 -52870.0 -71424.0 -43179.0"; // Xmin Ymin Xmax Ymax
            //msg += " DepthRange 0 100"; // min max
            //msg += " ScreenLenRange 10 70"; // min max
            //msg += " SourceArea 1 -9 -9 -9 "; Npixels minPixels maxPixels percPixels
			msg += " Ncrops 1";
			//msg += " Ncrops 13";
			//msg += " 301 0.5";
			//msg += " 302 0.5";
			//msg += " 303 0.4";
			//msg += " 400 0.6";
			//msg += " 401 0.45";
			//msg += " 402 0.6";
			//msg += " 605 0.3";
			//msg += " 1451 0.75";
			//msg += " 1452 0.253";
			//msg += " 1460 0.5501";
			//msg += " 212910 0.102";
			//msg += " 10003 0.458";
            msg += " -9 0.7";
            //msg += " PixelRadius 0";
			//msg += " DebugID test21";
			msg += " ENDofMSG\n";
		}
	}

    boost::asio::io_service io_service;

    boost::asio::ip::tcp::socket socket(io_service);
    socket.connect(boost::asio::ip::tcp::endpoint(boost::asio::ip::address::from_string("127.0.0.1"), 1234));
		
	boost::system::error_code error;
	boost::asio::write(socket, boost::asio::buffer(msg), error);

	if (!error) {
		std::cout << "Client sent a message" << std::endl;
	}
	else {
		std::cout << "sent failed" << std::endl;
	}

	boost::asio::streambuf receive_buffer;
	boost::asio::read(socket, receive_buffer, boost::asio::transfer_all(), error);

	if (error && error != boost::asio::error::eof) {
		std::cout << "receive failed: " << error.message() << std::endl;
	}
	else {
		const char* data = boost::asio::buffer_cast<const char*>(receive_buffer.data());
		//std::cout << data << std::endl;
		std::string str(data);
		//std::cout << str << std::endl;
		//std::stringstream ss(str.c_str());

		std::vector<std::string> str1 = split(str, " ");
		int ii = 0;
		//ss << data;
		int tf = std::atoi(str1[ii].c_str()); ii++;

        std::ofstream outstream;
        outstream.open(outfile.c_str());
		//ss >> tf;
		if (tf == 1) {
			int Nbtc = std::atoi(str1[ii].c_str()); ii++;
			int NsimYears = std::atoi(str1[ii].c_str()); ii++;
			std::cout << "Data size: " << Nbtc << " x " << NsimYears << std::endl;
			//ss >> Nbtc;

			std::cout << "Printing Results..." << std::endl;

			double dd;
			outstream << Nbtc << " " << NsimYears << std::endl;
			for (int i = 0; i < Nbtc; ++i) {
				for (int j = 0; j < NsimYears; ++j) {
					//ss >> dd;
					dd = std::atof(str1[ii].c_str()); ii++;
					outstream << dd << " ";
				}
				outstream << std::endl;
			}
		}
		else {
            outstream << 0 << " " << 0 << std::endl;
            outstream << data << std::endl;
			//std::cout << data << std::endl;
		}

        outstream.close();
	}
	return 0;
}
