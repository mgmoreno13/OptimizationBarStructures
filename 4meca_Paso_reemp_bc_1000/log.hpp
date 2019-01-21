#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <cstring>

using namespace std;

class Log {
private:
	time_t t;
	string base_name;
	ofstream file;
	bool is_open;
public:
	Log(string base_name);
	~Log();
	/*
		1: Aviso, 2: Alerta, 3: Error
	*/
	void write(string msg, int type);
	void close();
};

#endif