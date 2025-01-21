// MantisServer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
//#ifdef _WIN64
//#include "pch.h"
//#endif

//#define _USEHF 1

#include <iostream>

#include <boost/asio.hpp>
#include <boost/bind/bind.hpp>
//#include <functional>
#include <thread>
#include <chrono>
#include <ctime>

#include "MShelper.h"
#include "MSoptions.h"
#include "Nload.h"
#include "Rch.h"
#include "Region.h"
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
	if (!tf){
        return 0;
    }
	static mantisServer::Mantis M(msOptions);

    std::cout << "Running MantisServer Version : " << msOptions.version << std::endl;
    std::time_t result = std::time(nullptr);
    std::cout << "Starting Date : " << std::asctime(std::localtime(&result)) << std::endl;

    if (_USEHF>0){
        std::cout << "USE HDF5" << std::endl;
    }
    else{
        std::cout << "DONT USE HDF" << std::endl;
    }

    std::cout << "Reading input data ..." << std::endl;
	tf = M.readInputs();
	if (!tf)
		return 0;

    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "Mantis Server is Ready..." << std::endl;
    std::cout << "   Version: " << msOptions.version << std::endl;

	if (tf) {
		ba::io_service io_service;
		//ba::ip::tcp::acceptor acceptor(io_service, ba::ip::tcp::endpoint(ba::ip::tcp::v4(), msOptions.port));
        ba::ip::tcp::acceptor acceptor(io_service, ba::ip::tcp::endpoint(ba::ip::address::from_string(msOptions.ip), msOptions.port));
        std::cout << "   IP: " << acceptor.local_endpoint().address().to_string() << std::endl;
        std::cout << "   PORT: " << acceptor.local_endpoint().port() << std::endl;
		char buff[8192];
		while (true) {
			ba::ip::tcp::socket socket(io_service);
			acceptor.accept(socket);
            //std::cout << " After SOCKET" << std::endl;
			int bytes = static_cast<int>(ba::read(socket, ba::buffer(buff), boost::bind(read_complete, buff, boost::placeholders::_1, boost::placeholders::_2)));
            //std::cout << " After read_complete" << std::endl;
            std::string msg(buff, bytes);
			std::string outmsg;
			if (msg.compare("quit\n") == 0) {
				std::cout << "Received request to quit" << std::endl;
				outmsg = "0 Mantis Server is quiting...";
				socket.write_some(ba::buffer(outmsg));
				socket.close();
				break;
			}

			if (msg.compare("ping\n") == 0){
                std::cout << "Received a ping" << std::endl;
                outmsg = "0 pong";
                socket.write_some(ba::buffer(outmsg));
                socket.close();
                continue;
			}

            //std::cout << " MSG Recieved" << std::endl;
			bool bvalidMsg = M.parse_incoming_msg(msg, outmsg);
            //std::cout << " MSG Parsed" << std::endl;
			if (bvalidMsg) {
                //std::cout << "Validate MSG" << std::endl;
				bvalidMsg = M.validate_msg(outmsg);
                //std::cout << " MSG Validated" << std::endl;
			}
			if (bvalidMsg) {
                //std::cout << " Prepare Simulation" << std::endl;
				M.resetReply();
				auto start = std::chrono::high_resolution_clock::now();
				std::vector<std::thread> T;
				//M.simulate_with_threads(0);
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
			M.postReplyActions();
		}
	}

	return 0;
}