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

class Vertice
{
public:
    int index;
    double U, r, theta, fi, x, y, z;
    void ReadFromStr(const char* str);
};

class Tr
{
public:
    int index, fatherInd;
    vector<int> V = {-1, -1, -1};
    vector<int> childInd = {-1, -1, -1, -1};
    void ReadFromStr(const char* str);
};

class Extrapolation{
public:
    int stepRad;
    double maxRad;
    int trIdx;
    int vertAm;
    int trAm;
    void stepSet(int step);
    void radSet(double r);
    void verticeLoader(string filename);
    void triangleLoader(string filename);
    vector<vector<Vertice>> vert_arr;
    vector<Tr> tr_arr;
    double sphDist(double theta, double fi, int idx, int rIdx);
    void loader(int deg, int maxR);
    int zeroSearcher(int rIdx, double x, double y, double z);
    int searcher(int fthrIdx, int rIdx, double x, double y, double z);
    bool isInTr(int rIdx, int idx, double x, double y, double z);
    double layerCounter(double theta, double fi, int trIdx, int rIdx);
    double counter(double r, double theta, double fi, int tr, int rIdx);
    double extrapolator(double r, double theta, double fi);
};

#endif /* Extrapolation_hpp */
