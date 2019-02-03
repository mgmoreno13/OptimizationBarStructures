#include "log.hpp"

Log::Log(string base_name) {
	this->base_name = base_name;
	//Obtener el tiempo actual
	t = time(0);
	struct tm *now = localtime(&t);

	//Crear el archivo
	stringstream ss;
	ss << base_name << (now->tm_year + 1900) << "-";
	ss << (now->tm_mon + 1) << "-";
	ss << now->tm_mday << "-";
	ss << now->tm_hour << ".txt";
	string name = ss.str();

	file.open(name.c_str(), ofstream::out | ofstream::app);
	is_open = file.is_open();
	if(!file.is_open()) {
		cout << "Error al crear el archivo de log.\n";
	}
}

Log::~Log() {
	close();
}

void Log::write(string msg, int type) {
	//Obtener el tiempo actual
	t = time(0);
	struct tm *now = localtime(&t);
	stringstream ss;
	ss << (now->tm_year + 1900) << "-";
	ss << (now->tm_mon + 1) << "-";
	ss << now->tm_mday << " ";
	ss << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << " ";

	switch(type) {
		case 1:
			ss << "[Mensaje]: ";
			break;
		case 2:
			ss << "[Alerta]: ";
			break;
		case 3:
			ss << "[Error]: ";
			break;
		default:
			//Tipo no reconocido
			return;
	}

	ss << msg;
	file << ss.str() << "\n";
}

void Log::close() {
	if(file.is_open()) {
		file.close();
	}
}