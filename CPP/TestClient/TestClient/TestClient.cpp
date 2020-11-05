// TestClient.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#ifdef _WIN64
#include "pch.h"
//#endif

#include <iostream>
#include <string>
#include <fstream>
#include <boost/asio.hpp>
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
	std::string infile;
	if (argc > 1) {
		infile = argv[1];
		if (infile.compare("quit") == 0)
			quit = true;
	}

	boost::asio::io_service io_service;

	boost::asio::ip::tcp::socket socket(io_service);
	socket.connect(boost::asio::ip::tcp::endpoint(boost::asio::ip::address::from_string("127.0.0.1"), 1234));
	
	std::string msg;
	int NsimYears = 300;
	if (quit) {
		msg = "quit\n";
	}
	else {
		if (argc > 1) {
			std::ifstream indata(infile.c_str());
			if (!indata.good()) {
				std::cout << "Can't open the file " << infile << std::endl;
			}
			else {
				std::string line;
				bool firstLine = true;
				while (getline(indata, line)) {
					if (firstLine) {
						msg += line;
						firstLine = false;
					}
					else {
						msg += " " + line;
					}
				}
 				msg += " ENDofMSG\n";
			}
		}
		else {
			// Number of years to simulate, The year to start the reductions, the unsaturated zone mobile water content 
			msg = std::to_string(NsimYears);
			msg += " 2025 0.002";
			// Examples of the Second line line
			//CVHM_95_99 1 1 1 -> CVHM_95_99 scenario 1st Base, 1 region, with id 1 (The 1st base map has only one polygon
			//CVHM_95_99 2 1 3 -> CVHM_95_99 scenario 2st Base, 1 region, with id 3 (The second base map has 3 polygons (Subbasins) TLB has id 3)
			//CVHM_95_99 5 2 1 19 -> CVHM_95_99 scenario 5st Base, 2 regions, with id 1 amd 19 (The 5th base map has 21 polygons (farms) )
			msg += " CVHM_92_03_bud0 5 1 3"; // Scenario Name, MapID, Nregions, Region ids,
			msg += " 12"; // Number of categories for reduction
			msg += " 301 0.5";
			msg += " 302 0.5";
			msg += " 303 0.4";
			msg += " 400 0.7";
			msg += " 401 0.7";
			msg += " 402 0.6";
			msg += " 605 0.3";
			msg += " 1451 0.75";
			msg += " 1452 0.253";
			msg += " 1460 0.5501";
			msg += " 212910 0.102";
			msg += " 10003 0.458";
			msg += " ENDofMSG\n";
		}
	}
		
	boost::system::error_code error;
	boost::asio::write(socket, boost::asio::buffer(msg), error);

	if (!error) {
		std::cout << "Client sent hello message" << std::endl;
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
		std::cout << data << std::endl;
		std::string str(data);
		std::cout << str << std::endl;
		//std::stringstream ss(str.c_str());

		std::vector<std::string> str1 = split(str, " ");
		int ii = 0;
		//ss << data;
		int tf = std::atoi(str1[ii].c_str()); ii++;
		
		//ss >> tf;
		if (tf == 1) {
			int Nbtc = std::atoi(str1[ii].c_str()); ii++;
			NsimYears = std::atoi(str1[ii].c_str()); ii++;
			//ss >> Nbtc;
			std::string filename = "testClientResults.dat";
			std::ofstream outstream;
			outstream.open(filename.c_str());
			float dd;
			for (int i = 0; i < Nbtc; ++i) {
				for (int j = 0; j < NsimYears; ++j) {
					//ss >> dd;
					dd = std::atof(str1[ii].c_str()); ii++;
					outstream << dd << " ";
				}
				outstream << std::endl;
			}
			outstream.close();
		}
		else {
			std::cout << data << std::endl;
		}
	}

	return 0;

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
