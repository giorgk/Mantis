// TestClient.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <string>
#include <boost/asio.hpp>
#include <math.h>

int main()
{

	boost::asio::io_service io_service;

	boost::asio::ip::tcp::socket socket(io_service);
	socket.connect(boost::asio::ip::tcp::endpoint(boost::asio::ip::address::from_string("127.0.0.1"), 1234));

	std::string msg = "CVHM_95_99 5 1 10"; // Scenario Name, MapID, Nregions, Region ids,
	msg += " 12 2020"; // Number of categories for reduction year to start reduction
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
	msg += "\n";
		
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
