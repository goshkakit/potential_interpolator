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
    Triangulation tr = {};
    for(int deg = 5; deg<=8; deg++){
        ofstream infile("th on deg "+to_string(deg)+".txt");
        tr.map(deg, 10, r0+1000000, r0+1001000, 1000);
        for(int i=0; i<tr.tr_arr.size(); i++){
            infile<<abs(tr.vert_arr.at(tr.tr_arr[i].V[0]).theta - tr.vert_arr.at(tr.tr_arr[i].V[1]).theta)<<" "<<
                    abs(tr.vert_arr.at(tr.tr_arr[i].V[1]).theta - tr.vert_arr.at(tr.tr_arr[i].V[2]).theta)<<" "<<
                    abs(tr.vert_arr.at(tr.tr_arr[i].V[2]).theta - tr.vert_arr.at(tr.tr_arr[i].V[0]).theta)<<" ";
        }
    }
     */
    /*
    Triangulation tr = {};
    for(int i=3; i<=93; i+=10){
        tr.map(5, i, r0+1000000, r0+1001000, 1000);
    }
    
    PotentialCounter pc = {};
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    ofstream infile("dataset_"+to_string(r0+1000000)+".txt");
    std::uniform_real_distribution<> r_gen(r0+1000000, r0+1001000);
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(0, 2*pi);
    std::mt19937_64 gen((int)time(0));
    
    for(int i=0; i<=100000; i++){
        double r = r_gen(gen);
        double theta = th_gen(gen);
        double fi = fi_gen(gen);
        infile<<r<<" "<<theta<<" "<<fi<<endl;
    }
    
    for(int i=3; i<=93; i+=10){
        pc.length = i;
        string garm = to_string(i);
        Extrapolation ex = {};
        ex.loader(i, 5, r0+1001000, 1000, r0+1000000);
        double timeC = 0, timeE = 0;
        clock_t tC, tE;
        ofstream infileV("random extrapolated "+garm+".txt");
        ofstream infileC("random counted "+garm+".txt");
        string line;
        ifstream in("dataset_"+to_string(r0+1000000)+".txt");
        const char ch = '\n';
        char mass[200] = {};
        for(int i=0; i<100000; i++){
            memset(mass, 0, 200);
            in.getline(mass, 199, ch);
            double r, theta, fi;
            sscanf(mass, "%lf %lf %lf", &r, &theta, &fi);
            
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
    }
    */
    
    /*
    Triangulation tr = {};
    tr.map(4, 10, r0+20000000, r0+20001000, 1000);
    tr.map(5, 10, r0+20000000, r0+20001000, 1000);
    tr.map(6, 10, r0+20000000, r0+20001000, 1000);
    tr.map(7, 10, r0+20000000, r0+20001000, 1000);
    tr.map(8, 10, r0+20000000, r0+20001000, 1000);
    //ParallelTriangulation pTr = {};
    //pTr.parMap(2, 100, r0, r0+8000, 1000);
    
    PotentialCounter pc = {};
    pc.length = 10;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    double stR = r0+20000000;
    double enR = stR + 1000;
    ofstream infile("dataset_"+to_string(stR)+".txt");
    std::uniform_real_distribution<> r_gen(stR, enR);
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(0, 2*pi);
    std::mt19937_64 gen((int)time(0));
    
    for(int i=0; i<=100000; i++){
        double r = r_gen(gen);
        double theta = th_gen(gen);
        double fi = fi_gen(gen);
        infile<<r<<" "<<theta<<" "<<fi<<endl;
    }
    
    for(int degr=4; degr<=8; degr++){
        string degg = to_string(degr);
        
        Extrapolation ex = {};
        ex.loader(10, degr, enR, 1000, stR);
        double timeC = 0, timeE = 0;
        clock_t tC, tE;
        ofstream infileV("random extrapolated "+degg+".txt");
        ofstream infileC("random counted "+degg+".txt");
        string line;
        ifstream in("dataset_"+to_string(stR)+".txt");
        const char ch = '\n';
        char mass[200] = {};
        for(int i=0; i<100000; i++){
            memset(mass, 0, 200);
            in.getline(mass, 199, ch);
            double r, theta, fi;
            sscanf(mass, "%lf %lf %lf", &r, &theta, &fi);
            
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
    }
     */
    /*
    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    int stR = 6378100;
    int st = 2500000;
    int enR = 16378100;
     
    ofstream infile("dataset_"+to_string(stR)+".txt");
    std::uniform_real_distribution<> r_gen(stR, enR);
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(0, 2*pi);
    std::mt19937_64 gen((int)time(0));
    
    for(int i=0; i<=100000; i++){
        double r = r_gen(gen);
        double theta = th_gen(gen);
        double fi = fi_gen(gen);
        infile<<r<<" "<<theta<<" "<<fi<<endl;
    }
    
    Extrapolation ex = {};
    ex.loader(4, enR, st, stR);
    double timeC = 0, timeE = 0;
    clock_t tC, tE;
    ofstream infileV("random extrapolated " + to_string(st) + ".txt");
    ofstream infileC("random counted " + to_string(st) + ".txt");
    string line;
    ifstream in("dataset_"+to_string(stR)+".txt");
    const char ch = '\n';
    char mass[200] = {};
    for(int i=0; i<100000; i++){
        memset(mass, 0, 200);
        in.getline(mass, 199, ch);
        double r, theta, fi;
        sscanf(mass, "%lf %lf %lf", &r, &theta, &fi);
        
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
    cout<<timeC/100000<<" "<<timeE/timeC<<endl;
    
     */
    /*
    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    for(int degr=5; degr<=8; degr++){
        string degg = to_string(degr);
        
        Extrapolation ex = {};
        ex.loader(degr, 6579100, 1000, 6578100);
        
        std::uniform_real_distribution<> r_gen(6578100, 6579100);
        std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
        std::uniform_real_distribution<> fi_gen(0, 2*pi);
        std::mt19937_64 gen((int)time(0));
        double timeC = 0, timeE = 0;
        clock_t tC, tE;
        ofstream infileV("random extrapolated "+degg+".txt");
        ofstream infileC("random counted "+degg+".txt");
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
    }
     */
}
