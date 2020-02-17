//
//  Squarer.cpp
//  Project
//
//  Created by Георгий on 04.02.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "Squarer.hpp"

void Squarer::mesher(){
    /*
     r = [6378100, 6578100]
     theta = [-pi/2, pi/2]
     fi = [0, 2*pi]
     */
    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    vector<vector<vector<double> > > arr;
    
    ofstream infile("arr.txt");
    
    arr.resize(this->i_max);
    for (int i = 0; i < this->i_max; ++i) {
        arr[i].resize(this->j_max);
        for (int j = 0; j < this->j_max; ++j){
            arr[i][j].resize(this->k_max);
            for (int k = 0; k < this->k_max; ++k){
                arr[i][j][k] = pc.potential(r0 + i * this->r_step, -pi/2 + j * this->th_step, k * this->fi_step);
                cout<< i << " " << j << " " << k << " " << arr[i][j][k] << endl;
                infile << arr[i][j][k] << " ";
            }
        }
    }
    infile << endl;
    infile.close();
}

double Squarer::distance(double th1, double th2, double fi1, double fi2){
    return acos(sin(th1) * sin(th2) + cos(th1) * cos(th2) * cos(fi1 - fi2));
}

double Squarer::approxer(double r, double theta, double fi){
    double U=0;
    ifstream infile;
    infile.open("arr.txt");
    double num;
    while (infile >> num) {
        this->data.push_back(num);
    }
    
    int i, j, k, index;
    i = (int) ((r - r0)/this->r_step);
    j = (int) ((theta + pi/2)/this->th_step);
    k = (int) ((fi)/this->fi_step);
    index = i * 64800 + j * 360 + k;
    //if((r - r0)/this->r_step == i && (theta + pi/2)/this->th_step == j && (fi)/this->fi_step == k){
    //    U = this->data[index];
    //}
    //else{
    //    double U1, U2;
        
        //U = U2 * (r - i*this->r_step) / this->r_step + U1 * ((i+1)*this->r_step - r) / this->r_step;
    //}
    
    return U;
}
