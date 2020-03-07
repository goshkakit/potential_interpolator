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

int main(int argc, const char * argv[]) {
    //Triangulation tr = {};
    //tr.mesher(r0, 3);
    
    //for(int i=0; i<tr.vert_arr.size(); i++){
    //    cout<<tr.vert_arr[i].r<<" "<<tr.vert_arr[i].theta<<" "<<tr.vert_arr[i].fi<<" "<<tr.vert_arr[i].U<<endl;
    //}
    //tr.map(3, r0);
    //tr.map(6, 100, r0+4000, 100);
    //tr.zeroMeshCreator(7000000);
    //cout<<atan2(-1,-1)<<endl;
    /*
    for(int i=0; i<tr.tr_arr.size(); i++){
        cout<<tr.tr_arr[i].index<<" fthr: "<<tr.tr_arr[i].fatherInd<<" "<<tr.tr_arr[i].childInd[0]<<" "<<tr.tr_arr[i].childInd[1]<<" "<<tr.tr_arr[i].childInd[2]<<" "<<tr.tr_arr[i].childInd[3]<<endl;
    }
    cout<<"\n";
    for(int i=0; i<tr.vert_arr.size(); i++){
        cout<<tr.vert_arr[i].index<<" "<<tr.vert_arr[i].r<<" "<<tr.vert_arr[i].theta<<" "<<tr.vert_arr[i].fi<<" "<<endl;
    }
    */
    //cout<<tr.tr_arr.size()<<endl;
    //PotentialCounter pc = {};
    //pc.length = 100;
    //pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    //for(int r=r0; r<=r0+4000; r+=100){
    //for(int i=1; i<=5; i++){
    //    tr.map(6, r);
    //    tr.map(i, 8478100);
    //}
    //}
    
    Extrapolation ex = {};
    //"v_6_100_6381600_100.txt"
    //"Triangles_6.txt"
    
    ex.loader(6);
    
    cout<<ex.extrapolator(6378130, -0.0708, -3.0159)<<endl;
    //cout<<pc.potential(6378130, -0.0708, -3.0159)<<endl;
    /*
    for(int i=0; i<ex.vert_arr[0].size(); i++){
        cout<<ex.vert_arr[0][i].index<<" "<<ex.vert_arr[0][i].r<<" "<<ex.vert_arr[0][i].theta<<" "<<ex.vert_arr[0][i].fi<<" "<<ex.vert_arr[0][i].U<<endl;
    }
     
    for(int i=0; i<ex.tr_arr.size(); i++){
        cout<<ex.tr_arr[i].index<<" fthr: "<<ex.tr_arr[i].fatherInd<<" "<<ex.tr_arr[i].childInd[0]<<" "<<ex.tr_arr[i].childInd[1]<<" "<<ex.tr_arr[i].childInd[2]<<" "<<ex.tr_arr[i].childInd[3]<<" "<<ex.tr_arr[i].V[0]<<" "<<ex.tr_arr[i].V[1]<<
        " "<<ex.tr_arr[i].V[2]<<endl;
    }
     */
    /*
    double r = 6378100, theta, fi, rmax = 6387900;
    
    std::mt19937_64 gen(time(0));
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(-pi, pi);
    theta = th_gen(gen);
    fi = fi_gen(gen);
    
    ofstream infile("U_counted.txt");
    ofstream infileE("U_extrapolated.txt");
    int operations=0;
    time_t startC, endC;
    time(&startC);
    for(int i = r; i< rmax; i+= 1){
        operations++;
        double U_pc = pc.potential(i, theta, fi);
        infile<<U_pc<<" ";
    }
    time(&endC);
    double t_calc = difftime(endC, startC);
    time_t startE, endE;
    time(&startE);
    for(int i = r; i< rmax; i+= 1){
        double U_ex = ex.extrapolator(i, theta, fi);
        infileE<<U_ex<<" ";
    }
    time(&endE);
    double t_extr = difftime(endE, startE);
    
    cout<<"calc: "<<t_calc<<"; extr: "<<t_extr<<endl;
    double meanC = t_calc / operations;
    double meanE = t_extr / operations;
    cout<<operations<<endl;
    infile.close();
    infileE.close();
    
    cout<<"t_mean counted: "<<meanC<<"\n"<<"t_mean extrapolated: "<<meanE<<endl;
     */
}
