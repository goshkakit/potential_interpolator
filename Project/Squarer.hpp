//
//  Squarer.hpp
//  Project
//
//  Created by Георгий on 04.02.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#ifndef Squarer_hpp
#define Squarer_hpp

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
#include "PotentialCounter.h"

using namespace std;

#define pi 3.141592653589793238462643383279502884L

struct Squarer
{
    vector<double> data;
    double r_step = 100, th_step = pi/180, fi_step = pi/180;
    int i_max = (6578100 - r0) / r_step;
    int j_max = 180;
    int k_max = 360;
    void mesher();
    double distance(double th1, double th2, double fi1, double fi2);
    double approxer(double r, double theta, double fi);
};

#endif /* Squarer_hpp */
