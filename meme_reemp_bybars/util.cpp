#include "util.hpp"

string int_to_str(int n) {
	stringstream ss;
	ss << n;

	return ss.str();
}

vector<NODE> read_machinefile(string file) {
	vector<NODE> res;
	ifstream ifs(file.c_str());
	string line;
	if(ifs.is_open()) {
		while(getline(ifs, line)) {
			NODE n;
			n.name = line;
			n.in_use = false;
			n.pid = 0;
			res.push_back(n);
		}
		ifs.close();
	}
	else {
		cout << "Error al abrir el archivo: " << file << ".\n";
	}

	return res;
}

queue<pair<string, string> > read_tasksfile(string file) {
	queue<pair<string, string> > res;
	ifstream ifs(file.c_str());
	string line, first, second;
	if(ifs.is_open()) {
		while(getline(ifs, line)) {
			//Cortar la línea en dos
			string::size_type pos;
			pos = line.find(' ', 0);
			second = line.substr(pos + 1);
			first = line.substr(0, pos);

			pair<string, string> pp = make_pair(first, second);
			res.push(pp);
		}
		ifs.close();
	}
	else {
		cout << "Error al abrir el archivo: " << file << ".\n";
	}

	return res;
}

int get_next_free_node(vector<NODE> nodes) {
	for(int i = 0; i < nodes.size(); i++) {
		if(!nodes[i].in_use) {
			return i;
		}
	}

	return -1;
}