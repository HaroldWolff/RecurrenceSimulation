/// *** download the boost library from https://www.boost.org/ and make sure that the path is correct ***

#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#define BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <dirent.h>
#include <Windows.h>
#include <cassert>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "rndm.hpp"

using namespace std;

void GenerateOutputfile(int, int, int);
int GenerateNewPatient(int, int, string);
void WritePatient();
void CloseOutputfile();

#endif // SIMULATION_HPP_
