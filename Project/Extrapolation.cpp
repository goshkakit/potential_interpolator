//
//  Extrapolation.cpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "Extrapolation.hpp"

void Vertice::ReadFromStr(const char* str){
    sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf", &index, &U, &r, &theta, &fi, &x, &y, &z);
}

void Tr::ReadFromStr(const char *str){
    sscanf(str, "%d %d %d %d %d %d %d %d %d", &index, &fatherInd, &childInd[0], &childInd[1], &childInd[2], &childInd[3],
                                              &V[0], &V[1], &V[2]);
}

void Extrapolation::stepSet(int step){
    this->stepRad = step;
}

void Extrapolation::radSet(double r){
    this->maxRad = r;
}

void Extrapolation::verticeLoader(string filename){
    string line;
    ifstream in(filename);
    const char ch = '\n';
    char mass[100] = {};
    for(int rad = 0; rad <= int((maxRad - r0)/stepRad); rad++) {
        vector<Vertice> v = {};
        this->vert_arr.push_back(v);
        for(int i = 0; i <= this->vertAm; i++){
            Vertice newVert;
            memset(mass, 0, 100);
            in.getline(mass, 99, ch);
            newVert.ReadFromStr(mass);
            this->vert_arr[rad].push_back(newVert);
        }
    }
}

void Extrapolation::triangleLoader(string filename){
    string line;
    ifstream in(filename);
    const char ch = '\n';
    char mass[90] = {};
    for(int i=0; i<=this->trAm; i++){
        Tr tr = {};
        memset(mass, 0, 90);
        in.getline(mass, 89, ch);
        tr.ReadFromStr(mass);
        this->tr_arr.push_back(tr);
    }
}

double Extrapolation::sphDist(double theta, double fi, int idx, int rIdx){
    return acos(sin(theta) * sin(this->vert_arr.at(rIdx).at(idx).theta) +
                cos(theta) * cos(this->vert_arr.at(rIdx).at(idx).theta) *
                cos(fi - this->vert_arr.at(rIdx).at(idx).fi));
}

void Extrapolation::loader(int deg){
    string degree = to_string(deg);
    this->radSet(6388100);
    this->stepSet(100);
    if(deg == 6){
        this->vertAm = 40961;
        this->trAm = 109219;
    }
    if(deg == 7){
        this->vertAm = 163841;
        this->trAm = 436899;
    }
    this->verticeLoader("v_" + degree + "_100_100.txt");
    this->triangleLoader("Triangles_" + degree + ".txt");
}

int Extrapolation::zeroSearcher(int rIdx, double x, double y, double z){
    for(int i=0; i<20; i++){
        if(isInTr(rIdx, i, x, y, z))
            return i;
    }
    return -1;
}

int Extrapolation::searcher(int fthrIdx, int rIdx, double x, double y, double z){
    Tr& t = this->tr_arr.at(fthrIdx);
    if(t.childInd[0] != -1){
        for(int j=0; j<4; j++){
            if(isInTr(rIdx, t.childInd[j], x, y, z))
                return searcher(t.childInd[j], rIdx, x, y, z);
        }
    }
    else
        return fthrIdx;
    return -1;
}

bool Extrapolation::isInTr(int rIdx, int idx, double x, double y, double z){
    //double A, B, C, lmbd, l, m, V1x, V1y, V1z, V2x, V2y, V2z;
    
    Tr& t = this->tr_arr[idx];
    Vertice& V0 = this->vert_arr.at(rIdx).at(t.V[0]);
    Vertice& V1 = this->vert_arr.at(rIdx).at(t.V[1]);
    Vertice& V2 = this->vert_arr.at(rIdx).at(t.V[2]);
    
    //cout<<x<<" "<<y<<" "<<z<<"\n";
    //cout<<V0.x<<" "<<V0.y<<" "<<V0.z<<"\n";
    //cout<<V1.x<<" "<<V1.y<<" "<<V1.z<<"\n";
    //cout<<V2.x<<" "<<V2.y<<" "<<V2.z<<"\n";
    
    double delta = V0.x * (V1.y * V2.z - V2.y * V1.z) - V1.x * (V0.y * V2.z - V2.y * V0.z) + V2.x * (V0.y * V1.z - V1.y * V0.z);
    double deltaX = x * (V1.y * V2.z - V2.y * V1.z) - V1.x * (y * V2.z - V2.y * z) + V2.x * (y * V1.z - V1.y * z);
    double deltaY = V0.x * (y * V2.z - V2.y * z) - x * (V0.y * V2.z - V2.y * V0.z) + V2.x * (V0.y * z - y * V0.z);
    double deltaZ = V0.x * (V1.y * z - y * V1.z) - V1.x * (V0.y * z - y * V0.z) + x * (V0.y * V1.z - V1.y * V0.z);
    
    double lmbd = delta / (deltaX + deltaY + deltaZ);
    double a = deltaX / (deltaX + deltaY + deltaZ);
    double b = deltaY / (deltaX + deltaY + deltaZ);
    double c = deltaZ / (deltaX + deltaY + deltaZ);
    
    if(a > 0 && b > 0 && c > 0 && lmbd > 0)
        return true;
    
    /*
    A = ((V1.y - V0.y) * (V2.z - V0.z) - (V2.y - V0.y) * (V1.z - V0.z));
    B = ((V1.x - V0.x) * (V2.z - V0.z) - (V2.x - V0.x) * (V1.z - V0.z));
    C = ((V1.x - V0.x) * (V2.y - V0.y) - (V2.x - V0.x) * (V1.y - V0.y));
    
    lmbd = (V0.x * A - V0.y * B + V0.z * C) / (x * A - y * B + z * C);
    
    x *= lmbd;
    y *= lmbd;
    z *= lmbd;
    
    V1x = V1.x - V0.x;
    V1y = V1.y - V0.y;
    V1z = V1.z - V0.z;
    V2x = V2.x - V0.x;
    V2y = V2.y - V0.y;
    V2z = V2.z - V0.z;
    x -= V0.x;
    y -= V0.y;
    z -= V0.z;
    
    l = (x * V1x + y * V1y + z * V1z) / (V1x * V1x + V1y * V1y + V1z * V1z);
    
    if(l >= 0 && l <= 1){
        m = (x * V2x + y * V2y + z * V2z) / (V2x * V2x + V2y * V2y + V2z * V2z);
        if(m >= 0 && (l + m <= 1)){
            return true;
        }
    }
     */
    return false;
}

double Extrapolation::layerCounter(double theta, double fi, int trIdx, int rIdx){
    double U = 0;
    Tr& t = this->tr_arr.at(trIdx);
    double delta = 1e-10;
    double d1 = sphDist(theta, fi, t.V[0], rIdx) + delta;
    double d2 = sphDist(theta, fi, t.V[1], rIdx) + delta;
    double d3 = sphDist(theta, fi, t.V[2], rIdx) + delta;
    U =(this->vert_arr.at(rIdx).at(t.V[0]).U / d1 +
        this->vert_arr.at(rIdx).at(t.V[1]).U / d2 +
        this->vert_arr.at(rIdx).at(t.V[2]).U / d3) / (1 / d1 + 1 / d2 + 1 / d3);

    return U;
}

double Extrapolation::counter(double r, double theta, double fi, int tr, int rIdx){
    //.....U3..U1..U..U2..U4..
    if(rIdx > 0){
        double U1, U2, U3, U4, U, w1, w2;
        U1 = layerCounter(theta, fi, tr, rIdx);
        U2 = layerCounter(theta, fi, tr, rIdx + 1);
        U3 = layerCounter(theta, fi, tr, rIdx - 1);
        U4 = layerCounter(theta, fi, tr, rIdx + 2);
        w1 = (r - (r0 + (stepRad * rIdx))) / stepRad;
        w2 = (r0 + (stepRad * (rIdx + 1)) - r) / stepRad;
        U = (U2 * w1 + U1 * w2) / (w1 + w2);
        
        double delta = 1e-10, d1, d2, d3, d4;
        d1 = 1 / (abs(U   -   U1) + delta);
        d2 = 1 / (abs(U2  -   U)  + delta);
        d3 = 1 / (abs(U   -   U3) + delta);
        d4 = 1 / (abs(U4  -   U)  + delta);
        
        U = (U4 / d4 + U3 / d3 + U2 / d2 + U1 / d1) / (1 / d1 + 1 / d2 + 1 / d3 + 1 / d4);
        return U;
    }
    else{
        double delta = 1e-10, U1, U2, U4, U, w1, w2;
        U1 = layerCounter(theta, fi, tr, rIdx);
        U2 = layerCounter(theta, fi, tr, rIdx + 1);
        U4 = layerCounter(theta, fi, tr, rIdx + 2);
        w1 = 1 / (abs(r - (r0 + (stepRad * rIdx))) + delta);
        w2 = 1 / (abs(r0 + (stepRad * (rIdx + 1)) - r) + delta);
        U = (U2 * w1 + U1 * w2) / (w1 + w2);
        
        double d1, d2, d4;
        d1 = 1 / (abs(U   -   U1) + delta);
        d2 = 1 / (abs(U2  -   U)  + delta);
        d4 = 1 / (abs(U4  -   U)  + delta);
        
        U = (U4 / d4 + U2 / d2 + U1 / d1) / (1 / d1 + 1 / d2 + 1 / d4);
        return U;
    }
}

double Extrapolation::extrapolator(double r, double theta, double fi){
    double U = 0;
    double x = r * sin(theta) * cos(fi);
    double y = r * sin(theta) * sin(fi);
    double z = r * cos(theta);
    int rIdx = int(round((r-r0)/stepRad));
    int zeroIdx = zeroSearcher(rIdx, x, y, z);
    int tr = searcher(zeroIdx, rIdx, x, y, z);
    
    U = counter(r, theta, fi, tr, rIdx);
    return U;
}
