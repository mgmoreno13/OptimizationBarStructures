#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <queue>

using namespace std;

typedef struct node {
	string name;
	bool in_use;
	int pid;
} NODE;

string int_to_str(int n);

vector<NODE> read_machinefile(string file);

queue<pair<string, string> > read_tasksfile(string file);

int get_next_free_node(vector<NODE> nodes);

#endif