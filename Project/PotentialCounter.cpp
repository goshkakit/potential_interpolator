//
//  PotentialCounter.cpp
//  Project
//
//  Created by Георгий on 30.11.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#include "PotentialCounter.h"

void data_row::ReadFromStr(const char* str)
{
    sscanf(str, "%d %d %le %le", &n, &m, &a, &b);
}

void PotentialCounter::LoadFromFile(string filename){
    string line;
    ifstream in(filename);
    if (!in)
        cout << "Error!" << endl;
    else
    {
        const long len = 90, strings = 65338;
        const char ch = '\n';
        char mass[len] = {};
        this->data = new data_row[strings];
        for(int r = 0; r<strings; r++)
        {
            memset(mass, 0, len);
            in.getline(mass, len-1,ch);
            data[r].ReadFromStr(mass);
        }
    }
}

void PotentialCounter::cleanData(){
    delete[] data;
    data = nullptr;
}

double PotentialCounter::Cnm(int n,int m){
    int index = (n+3)*(n-2)/2 + m ;
    return data[index].a;
}

double PotentialCounter::Snm(int n,int m){
    int index = (n+3)*(n-2)/2 + m ;
    return data[index].b;
}

long double PotentialCounter::coefCounter(int n, int m)
{
    long double coef = 2.0*(2.0*n+1.0);
    if(m != 0){
        for(int i = n-m+1; i<=n+m; i++){
            coef/=i;
        }
    }
    return sqrt(abs(coef));
}


void PotentialCounter::LegendreCounter(int n, int m, long double x){
    double fact;
    int j;
    int k;
    legArr = new double[n+1];
    for (j = 0; j < n + 1; j++)
    {
        legArr[j] = 0.0;
    }
    if (m <= n)
    {
        legArr[m] = 1.0;
        fact = 1.0;
        for (k = 0; k < m; k++)
        {
            legArr[m] = - legArr[m] * fact * sqrt (1.0 - x * x);
            fact = fact + 2.0;
        }
    }
    if (m + 1 <= n)
    {
        legArr[m+1] = x * (double) (2*m+1) * legArr[m];
    }
    for (j = m+2; j <= n; j++)
     {
         legArr[j] = ((double)(2 * j     - 1) * x * legArr[j-1]
                     +(double)(  - j - m + 1) *     legArr[j-2] )
                     /(double)(    j - m    );
     }
}

long double PotentialCounter::partSumCounter(double r, double fi, double lmbd, int n, int m){
    this->LegendreCounter(n,m,sin(fi));
    long double partSum = 0;
    partSum = pow((r0/r),(n+1)) * legArr[n] *
              (Cnm(n, m) * cos(m * lmbd) + Snm(n, m) * sin(m * lmbd)) * this->coefCounter(n, m);
    this->cleanLeg();
    return partSum;
}

void PotentialCounter::cleanLeg(){
    delete[] legArr;
    legArr = nullptr;
}

long double PotentialCounter::potential(double r, double fi, double lmbd){
    long double sum = 0;
    for (int n=3; n<=length; n++) {
        for(int m=0; m<=n; m++){
            sum += this->partSumCounter(r, fi, lmbd, n, m);
        }
    }
    return f_c*M_c*sum/r0;
}

void PotentialCounter::Umap(){
    time_t start, end;
    time(&start);
    const int width = 200;
    const int height = 400;
    long double arr[width][height] = {0};
    ofstream infile("U_map.txt");
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            arr[i][j] = this->potential(7000000, -pi/2 + i*pi/(width-1), -pi + j*pi/((height-1)/2));
            printf("%d %d %Lf\n", i, j, arr[i][j]);
            infile << arr[i][j]<< " ";
        }
        infile << endl;
    }
    infile.close();
    time(&end);
    printf("Completed in %f min\n", difftime(end, start)/60.0);
}
