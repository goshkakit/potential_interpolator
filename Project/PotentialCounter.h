//
//  PotentialCounter.h
//  Project
//
//  Created by Георгий on 30.11.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#ifndef PotentialCounter_h
#define PotentialCounter_h
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <time.h>

using namespace std;

const long double f_c = 6.67408e-11;
const long double M_c = 5.97236473e+24L;
const long double r0 = 6378100;
#define pi 3.141592653589793238462643383279502884L

struct data_row
{
    int n, m;
    double a, b;
    void ReadFromStr(const char* str);
};

struct PotentialCounter
{
    int length;
    data_row* data;
    double* legArr;
    void LoadFromFile(string filename);
    void cleanData();
    double Cnm(int n,int m);
    double Snm(int n,int m);
    long double coefCounter(int n, int m);
    void LegendreCounter(int n, int m, long double x);
    void cleanLeg();
    long double partSumCounter(double r, double fi, double lmbd, int n, int m);
    long double potential(double r, double fi, double lmbd);
    void Umap();
};



#endif /* PotentialCounter_h */
