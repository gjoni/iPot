
#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>

struct OPTS {
	std::string rec;
	std::string lig;
	std::string table;
	double dmax;
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const std::string &exec, const double DMAX);

#endif /* OPTIONS_H_ */
