// MantisServer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <thread>
#include <chrono>
#include <ctime>

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
	if (!tf)
		return 0;

	
	static mantisServer::Mantis M(msOptions);

	std::cout << "Reading input data ..." << std::endl;
	tf = M.readInputs();
	if (!tf)
		return 0;


	std::cout << "Mantis Server Ready..." << std::endl;
	if (tf) {
		ba::io_service io_service;
		ba::ip::tcp::acceptor acceptor(io_service, ba::ip::tcp::endpoint(ba::ip::tcp::v4(), msOptions.port));
		char buff[1024];
		while (true) {
			ba::ip::tcp::socket socket(io_service);
			acceptor.accept(socket);
			int bytes = static_cast<int>(ba::read(socket, ba::buffer(buff), boost::bind(read_complete, buff, _1, _2)));
			std::string msg(buff, bytes);
			std::string outmsg;
			if (msg.compare("quit\n") == 0) {
				std::cout << "Received request to quit" << std::endl;
				outmsg = "Mantis Server is quiting...";
				socket.write_some(ba::buffer(outmsg));
				socket.close();
				break;
			}
				

			bool bvalidMsg = M.parse_incoming_msg(msg, outmsg);
			if (bvalidMsg) {
				M.resetReply();
				auto start = std::chrono::high_resolution_clock::now();
				std::vector<std::thread> T;
				for (int i = 0; i < msOptions.nThreads; ++i)
					T.push_back(std::thread(&mantisServer::Mantis::simulate_with_threads, std::ref(M), i));
				for (int i = 0; i < msOptions.nThreads; ++i)
					T[i].join();

				auto finish = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = finish - start;
				std::cout << "Simulation executed in " << elapsed.count() << std::endl;

				M.makeReply(outmsg);
				std::cout << "out message length: " << outmsg.size() << std::endl;
				std::cout << "Stop here " << std::endl;

			}
			//M.simulate(msg, outmsg);
			socket.write_some(ba::buffer(outmsg));
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
