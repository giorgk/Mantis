// MantisServer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include <boost/asio.hpp>
#include <boost/bind.hpp>

#include "MSoptions.h"
#include "MantisMain.h"

namespace ba = boost::asio;

size_t read_complete(char * buff, const boost::system::error_code & err, size_t bytes) {
	if (err) return 0;
	bool found = std::find(buff, buff + bytes, '\n') < buff + bytes;
	// we read one-by-one until we get to enter, no buffering
	return found ? 0 : 1;
}


int main(int argc, char *argv[])
{
	mantisServer::options msOptions;
	bool tf = mantisServer::readInputParameters(argc, argv, msOptions);

	
	mantisServer::Mantis M(msOptions);

	std::cout << "Reading input data ..." << std::endl;
	tf = M.readInputs();

	std::cout << "Mantis Server Ready..." << std::endl;
	if (tf) {
		ba::io_service io_service;
		ba::ip::tcp::acceptor acceptor(io_service, ba::ip::tcp::endpoint(ba::ip::tcp::v4(), msOptions.port));
		char buff[1024];
		while (true) {
			ba::ip::tcp::socket socket(io_service);
			acceptor.accept(socket);
			int bytes = ba::read(socket, ba::buffer(buff), boost::bind(read_complete, buff, _1, _2));
			std::string msg(buff, bytes);
			M.parse_input(msg);
			socket.close();
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
