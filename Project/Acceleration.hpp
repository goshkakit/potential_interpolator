//
//  Acceleration.hpp
//  Project
//
//  Created by Георгий on 08.03.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#ifndef Acceleration_hpp
#define Acceleration_hpp

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
#include "Extrapolation.hpp"

using namespace std;

const double J2 = 1.082628e-3;

class Acceleration{
public:
    double dR, dTh, dFi;
    void dSet(double dR, double dTh, double dFi);
    void countLoader();
    void extrLoader();
    PotentialCounter pc;
    Extrapolation ex;
    double U0(double r);
    double U2(double r, double theta);
    double accR_C(double r, double theta, double fi);
    double accR_E(double r, double theta, double fi);
    double accTh_C(double r, double theta, double fi);
    double accTh_E(double r, double theta, double fi);
    double accFi_C(double r, double theta, double fi);
    double accFi_E(double r, double theta, double fi);
    double g_C(double r, double theta, double fi);
    double g_E(double r, double theta, double fi);
};

#endif /* Acceleration_hpp */
