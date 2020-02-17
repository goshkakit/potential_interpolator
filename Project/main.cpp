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
#include "Squarer.hpp"

int main(int argc, const char * argv[]) {
    Triangulation tr = {};
    tr.mesher(r0, 1);
    
    //for(int i=0; i<tr.vert_arr.size(); i++){
    //    cout<<tr.vert_arr[i].r<<" "<<tr.vert_arr[i].theta<<" "<<tr.vert_arr[i].fi<<" "<<tr.vert_arr[i].U<<endl;
    //}
    //tr.map(3, r0);
    //tr.globalMapDot(1, 100, r0+100, 100);
    //tr.zeroMeshCreator(7000000);
    for(int i=0; i<tr.tr_arr.size(); i++){
        cout<<tr.tr_arr[i].index<<" fthr: "<<tr.tr_arr[i].fatherInd<<" "<<tr.tr_arr[i].V[0]<<" "<<tr.tr_arr[i].V[1]->index<<" "<<tr.tr_arr[i].V[2]->index<<" "<<endl;
    //    cout<<tr.tr_arr[i].index<<" "<<tr.tr_arr[i].V[0]->theta
    //                            <<" "<<tr.tr_arr[i].V[1]->theta
    //                            <<" "<<tr.tr_arr[i].V[2]->theta<<" "<<endl;
    }
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
    /*
    Extrapolation ex = {};
    
    double calc, t_calc=0, extr, t_extr=0;
    double r = 6378100, theta, fi;
    
    std::mt19937_64 gen(time(0));
    std::uniform_real_distribution<> th_gen(-pi/2, pi/2);
    std::uniform_real_distribution<> fi_gen(-pi, pi);
    theta = th_gen(gen);
    fi = fi_gen(gen);
    cout<<theta<<" "<<fi<<endl;
    
    ofstream infile(to_string(theta)+"_"+to_string(fi)+".txt");
    while(r<r0 + 4000){
        time_t start, end, e_start, e_end;
        time(&start);
        calc = pc.potential(r, theta, fi);
        time(&end);
        t_calc+=difftime(end, start);
        
        time(&e_start);
        extr = ex.extrapolator(6, r, theta, fi);
        time(&e_end);
        t_extr+=difftime(e_end,e_start);
        
        infile<<r<<" "<<calc<<" "<<extr<<"\n";
        cout<<"rad: "<<r<<endl;
        r+=10;
    }
    infile.close();
    cout<<"Calc time: "<<t_calc<<", extr time: "<<t_extr<<endl;
    */
    //Squarer sq = {};
    //cout<<sq.approxer(r0+2000, pi/31, pi/61)<<endl;
}
