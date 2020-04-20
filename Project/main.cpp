//
//  main.cpp
//  Project
//
//  Created by Георгий on 07.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include <stdio.h>
#include "Triangulation.hpp"
#include "PotentialCounter.h"
#include "Extrapolation.hpp"
#include "ParallelTriangulation.hpp"

int main(int argc, const char * argv[]) {
    /*
    std::uniform_real_distribution<> r_gen(r0, r0+100000000);
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(0, 2*pi);
    std::mt19937_64 gen((int)time(0));
    PotentialCounter pc = {};
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    ofstream infile("dots_N.txt");
    vector<double> arr;
    for(int len = 5; len<=100; len+=1){
        cout<<len<<endl;
        pc.length = len;
        clock_t t;
        t = clock();
        for(int i=0; i<100000; i++){
            double r = r_gen(gen);
            double th = th_gen(gen);
            double fi = fi_gen(gen);
            pc.potential(r, th, fi);
        }
        t = clock() - t;
        double dt = (double)t / CLOCKS_PER_SEC;
        dt/=100000;
        arr.push_back(dt);
    }
    for(int i=0; i<arr.size(); i++){
        infile<<arr.at(i)<<" ";
    }
     */
    /*
    Extrapolation ex = {};
    ex.loader(8, 6579100, 1000, 6578100);
    cout<<ex.extrapolator(6578503.7310264464, 0.24891772244761756, 3.5650866011713163)[0]<<endl;
     */
    //Triangulation tr = {};
    //tr.map(8, 100, r0+1000000, r0+1001000, 1000);
    //ParallelTriangulation pTr = {};
    //pTr.parMap(2, 100, r0, r0+8000, 1000);
    
    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    Extrapolation ex = {};
    ex.loader(8, 6579100, 1000, 6578100);
    
    std::uniform_real_distribution<> r_gen(6578100, 6579100);
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(0, 2*pi);
    std::mt19937_64 gen((int)time(0));
    double timeC = 0, timeE = 0;
    clock_t tC, tE;
    ofstream infileV("random extrapolated.txt");
    ofstream infileC("random counted.txt");
    cout<<"Counting potentials"<<endl;
    for(int i=0; i<10000; i++){
        double r = r_gen(gen);
        double theta = th_gen(gen);
        double fi = fi_gen(gen);
        
        tE = clock();
        vector<double> eex = ex.extrapolator(r, theta, fi);
        tE = clock() - tE;
        timeE += (double)tE / CLOCKS_PER_SEC;
        infileV<<eex[0]<<" "<<eex[1]<<" "<<eex[2]<<" "<<eex[3]<<" ";
        
        tC = clock();
        double U_pc = pc.potential(r, theta, fi);
        double aR = pc.accR(r, theta, fi);
        double aTh = pc.accTh(r, theta, fi);
        double aFi = pc.accFi(r, theta, fi);
        tC = clock() - tC;
        timeC += (double)tC / CLOCKS_PER_SEC;
        infileC<<U_pc<<" "<<aR<<" "<<aTh<<" "<<aFi<<" ";
    }
    cout<<timeE/timeC<<endl;
    
    /*
    clock_t t, te;
    t = clock();
    cout<<ex.extrapolator(6378900, 0.354, -3.018)[3]<<endl;
    t = clock() - t;
    cout<<(double)t / CLOCKS_PER_SEC<<endl;
    te = clock();
    cout<<pc.accFi(   6378900, 0.354, -3.018)<<endl;
    te = clock() - te;
    cout<<(double)te / CLOCKS_PER_SEC<<endl;
     */
    /*
    for(int i=0; i<ex.vert_arr[0].size(); i++){
        cout<<ex.vert_arr[0][i].index<<" "<<ex.vert_arr[0][i].r<<" "<<ex.vert_arr[0][i].theta<<" "<<ex.vert_arr[0][i].fi<<" "<<ex.vert_arr[0][i].U<<endl;
    }
     
    for(int i=0; i<ex.tr_arr.size(); i++){
        cout<<ex.tr_arr[i].index<<" fthr: "<<ex.tr_arr[i].fatherInd<<" "<<ex.tr_arr[i].childInd[0]<<" "<<ex.tr_arr[i].childInd[1]<<" "<<ex.tr_arr[i].childInd[2]<<" "<<ex.tr_arr[i].childInd[3]<<" "<<ex.tr_arr[i].V[0]<<" "<<ex.tr_arr[i].V[1]<<
        " "<<ex.tr_arr[i].V[2]<<endl;
    }
     */
}
