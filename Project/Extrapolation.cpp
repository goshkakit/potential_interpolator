//
//  Extrapolation.cpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "Extrapolation.hpp"

void Vertice::ReadFromStr(const char* str){
    sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &index, &U, &r, &theta, &fi, &x, &y, &z, &accR, &accTh, &accFi);
}

void Tr::ReadFromStr(const char *str){
    sscanf(str, "%d %d %d %d %d %d %d %d %d", &index, &fatherInd, &childInd[0], &childInd[1], &childInd[2], &childInd[3],
                                              &V[0], &V[1], &V[2]);
}

void Extrapolation::stepSet(int step){
    this->stepRad = step;
}

void Extrapolation::radSet(int stR, int maxR){
    this->startRad = stR;
    this->maxRad = maxR;
}

void Extrapolation::verticeLoader(string filename){
    string line;
    ifstream in(filename);
    const char ch = '\n';
    char mass[120] = {};
    for(int rad = 0; rad <= int((maxRad - r0)/stepRad); rad++) {
        vector<Vertice> v = {};
        this->vert_arr.push_back(v);
        for(int i = 0; i <= this->vertAm; i++){
            Vertice newVert;
            memset(mass, 0, 120);
            in.getline(mass, 119, ch);
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

void Extrapolation::loader(int deg, int maxR, int step, int stR){
    string degree = to_string(deg);
    string st = to_string(step);
    this->radSet(stR, maxR);
    this->stepSet(step);
    if(deg == 5){
        this->vertAm = 10241;
        this->trAm = 27299;
    }
    if(deg == 6){
        this->vertAm = 40961;
        this->trAm = 109219;
    }
    if(deg == 7){
        this->vertAm = 163841;
        this->trAm = 436899;
    }
    if(deg == 8){
        this->vertAm = 655361;
        this->trAm = 1747619;
    }
    if(startRad == r0)
        this->verticeLoader("v_A" + degree + "_100_" + st +".txt");
    else
        this->verticeLoader("v_A" + degree + "_100_" + st+ "_" + to_string(startRad) + "-" + to_string(maxRad) + ".txt");
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
    int count = 0;
    if(t.childInd[0] != -1){
        for(int j=0; j<4; j++){
            if(isInTr(rIdx, t.childInd[j], x, y, z))
                return searcher(t.childInd[j], rIdx, x, y, z);
            else
                count++;
        }
    }
    else
        return fthrIdx;
    return fthrIdx;
}

bool Extrapolation::isInTr(int rIdx, int idx, double x, double y, double z){
    Tr& t = this->tr_arr[idx];
    Vertice& V0 = this->vert_arr.at(rIdx).at(t.V[0]);
    Vertice& V1 = this->vert_arr.at(rIdx).at(t.V[1]);
    Vertice& V2 = this->vert_arr.at(rIdx).at(t.V[2]);
    
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
    return false;
}

vector<double> Extrapolation::layerCounter(double theta, double fi, int trIdx, int rIdx){
    double U = 0;
    vector<double> res;
    Tr& t = this->tr_arr.at(trIdx);
    double delta = 1e-10;
    double d1 = sphDist(theta, fi, t.V[0], rIdx) + delta;
    double d2 = sphDist(theta, fi, t.V[1], rIdx) + delta;
    double d3 = sphDist(theta, fi, t.V[2], rIdx) + delta;
    U = (this->vert_arr.at(rIdx).at(t.V[0]).U / d1 +
         this->vert_arr.at(rIdx).at(t.V[1]).U / d2 +
         this->vert_arr.at(rIdx).at(t.V[2]).U / d3) / (1 / d1 + 1 / d2 + 1 / d3);
    res.push_back(U);
    
    double accR = (this->vert_arr.at(rIdx).at(t.V[0]).accR / d1 +
                   this->vert_arr.at(rIdx).at(t.V[1]).accR / d2 +
                   this->vert_arr.at(rIdx).at(t.V[2]).accR / d3) / (1 / d1 + 1 / d2 + 1 / d3);
    res.push_back(accR);
    double accTh = (this->vert_arr.at(rIdx).at(t.V[0]).accTh / d1 +
                    this->vert_arr.at(rIdx).at(t.V[1]).accTh / d2 +
                    this->vert_arr.at(rIdx).at(t.V[2]).accTh / d3) / (1 / d1 + 1 / d2 + 1 / d3);
    res.push_back(accTh);
    double accFi = (this->vert_arr.at(rIdx).at(t.V[0]).accFi / d1 +
                    this->vert_arr.at(rIdx).at(t.V[1]).accFi / d2 +
                    this->vert_arr.at(rIdx).at(t.V[2]).accFi / d3) / (1 / d1 + 1 / d2 + 1 / d3);
    res.push_back(accFi);
    return res;
}

vector<double> Extrapolation::counter(double r, double theta, double fi, int tr, int rIdx){
    vector<double> res, V1, V2;
    double delta = 1e-10;
    V1 = layerCounter(theta, fi, tr, rIdx);
    double U1 = V1[0];
    double accR1 = V1[1];
    double accTh1 = V1[2];
    double accFi1 = V1[3];
    
    if(r != maxRad){
        V2 = layerCounter(theta, fi, tr, rIdx + 1);
        double U2 = V2[0];
        double accR2 = V2[1];
        double accTh2 = V2[2];
        double accFi2 = V2[3];
        double w1 = 1 / (abs(r - (r0 + (stepRad * rIdx))) + delta);
        double w2 = 1 / (abs(r0 + (stepRad * (rIdx + 1)) - r) + delta);
        double U = (U2 * w1 + U1 * w2) / (w1 + w2);
        double accR = (accR2 * w1 + accR1 * w2) / (w1 + w2);
        double accTh = (accTh2 * w1 + accTh1 * w2) / (w1 + w2);
        double accFi = (accFi2 * w1 + accFi1 * w2) / (w1 + w2);
        res.push_back(U);
        res.push_back(accR);
        res.push_back(accTh);
        res.push_back(accFi);
    }else{
        res.push_back(U1);
        res.push_back(accR1);
        res.push_back(accTh1);
        res.push_back(accFi1);
    }
    return res;
    /*
    double d1, d2, d3, d4;
    d1 = 1 / (abs(U   -   U1) + delta);
    d2 = 1 / (abs(U2  -   U)  + delta);
    d4 = 1 / (abs(U4  -   U)  + delta);
    
    if(rIdx > 0){
        U3 = layerCounter(theta, fi, tr, rIdx - 1);
        d3 = 1 / (abs(U   -   U3) + delta);
        U = (U4 / d4 + U3 / d3 + U2 / d2 + U1 / d1) / (1 / d1 + 1 / d2 + 1 / d3 + 1 / d4);
    }
    else
        U = (U4 / d4 + U2 / d2 + U1 / d1) / (1 / d1 + 1 / d2 + 1 / d4);
     */
    
}

vector<double> Extrapolation::extrapolator(double r, double theta, double fi){
    vector<double> V;
    theta += pi/2;
    //For mapping uncomment:
    //fi -= pi;
    double x = r * sin(theta) * cos(fi);
    double y = r * sin(theta) * sin(fi);
    double z = r * cos(theta);
    int rIdx = int(round((r-startRad)/stepRad));
    int zeroIdx = zeroSearcher(rIdx, x, y, z);
    int tr = searcher(zeroIdx, rIdx, x, y, z);
    
    V = counter(r, theta, fi, tr, rIdx);
    return V;
}
