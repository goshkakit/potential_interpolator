//
//  Triangulation.hpp
//  Project
//
//  Created by Георгий on 09.12.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#ifndef Triangulation_hpp
#define Triangulation_hpp
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
#include <algorithm>
#include "PotentialCounter.h"

using namespace std;

#define pi 3.141592653589793238462643383279502884L

class VerticeWithCompared
{
public:
    double x,y,z,r,theta,fi,U;
    int index;
    VerticeWithCompared *compared;
    VerticeWithCompared(){};
    void toSph();
    bool operator < (const VerticeWithCompared &rht) const;
};

class Triangle
{
public: vector<int> childInd ;
        Triangle();
        Triangle(int index, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3, int fatherIndex = -1);
        int index;
        int fatherInd;
        bool isDone = false;
        //vector<int> V;
        vector<VerticeWithCompared*> V;
        bool operator == (const Triangle &tr) const;
//private:
  //      Triangle(const Triangle& T);
};

class Triangulation
{
public:
    double r;
    void setRad(double r);
    static double distCounter(VerticeWithCompared V1, VerticeWithCompared V2);
    void zero_triangles();
    void zeroMeshCreator(double r);
    int vIndex;
    int globalIndex;
    int localVIndex;
    int localTrIndex;
    void mesher(double r, int degree);
    void mesherDot(double r, int degree);
    vector<VerticeWithCompared> vert_arr;
    vector<Triangle> tr_arr;
    vector<VerticeWithCompared> closestVerticesFinder(VerticeWithCompared V, vector<VerticeWithCompared> tmpVertices, int mode);
    VerticeWithCompared verticeCreator(VerticeWithCompared V1, VerticeWithCompared V2);
    bool vertDetector(VerticeWithCompared V);
    void triangledetector(Triangle* tr, vector<Triangle> trArr);
    vector<Triangle> trCreator(Triangle* fthrTr, VerticeWithCompared** V);
    void map(int degree, double r);
    void globalMapDot(int degree, int polynomDegree, int maxRad, int step);
};

#endif /* Triangulation_hpp */
