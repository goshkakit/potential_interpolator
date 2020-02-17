//
//  Extrapolation.hpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#ifndef Extrapolation_hpp
#define Extrapolation_hpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <random>
#include <ctime>
#include "PotentialCounter.h"
#include "Triangulation.hpp"

using namespace std;

class Vertice_data
{
public:
    double U, r, theta, fi;
    void ReadFromStr(const char* str);
};

struct Extrapolation{
    
    
    vector<Vertice_data> LoadFromFile(string filename);
    
    double distance(double theta, double fi, Vertice_data V);
    double triangle_counter(vector<Vertice_data> vect, double theta, double fi);
    double extrapolator(int deg, double r, double theta, double fi);
};

#endif /* Extrapolation_hpp */
