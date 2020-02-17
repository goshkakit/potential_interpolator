//
//  Extrapolation.cpp
//  Project
//
//  Created by Георгий on 04.01.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "Extrapolation.hpp"

void Vertice_data::ReadFromStr(const char* str)
{
    sscanf(str, "%lf %lf %lf %lf", &U, &r, &theta, &fi);
}

vector<Vertice_data> Extrapolation::LoadFromFile(string filename){
    string line;
    vector<Vertice_data> vect;
    ifstream in(filename);
    const char ch = '\n';
    char mass[90] = {};
    while (!in.eof()) {
        Vertice_data newVert;
        memset(mass, 0, 90);
        in.getline(mass, 89, ch);
        newVert.ReadFromStr(mass);
        vect.push_back(newVert);
    }
    return vect;
}
/*
vector<string> Extrapolation::name_searcher(double r){
    vector<string> name;
    
    
    int num = r/100000;
    if(num < 85){
        int res = (int)r%100000;
        if(res < 28100){
            name.push_back(to_string(num-1) + "78100_");
            name.push_back(to_string(num) + "28100_");
        }
        else if(res >=28100 && res < 78100){
            name.push_back(to_string(num) + "28100_");
            name.push_back(to_string(num) + "78100_");
        }
        else{
            if(num != 84){
                name.push_back(to_string(num) + "78100_");
                name.push_back(to_string(num+1) + "28100_");
            }
            else{
                name.push_back(to_string(num) + "78100_");
                name.push_back("8500000_");
            }
        }
    }
    else if (num >= 85 && num < 100){
        name.push_back(to_string(num) + "00000_");
        name.push_back(to_string(num+1) + "00000_");
    }
    else if (num >= 100 && num < 150){
        int snum = (int)num/10;
        int res = (int)num%10;
        if(res < 5){
            name.push_back(to_string(snum) + "000000_");
            name.push_back(to_string(snum) + "500000_");
        }
        else{
            name.push_back(to_string(snum) + "500000_");
            name.push_back(to_string(snum+1) + "000000_");
        }
    }
    else if(num >= 150 && num < 400){
        name.push_back(to_string(num) + "00000_");
        name.push_back(to_string(num+10) + "00000_");
    }
    //errors??
    
    return name;
}
*/
double Extrapolation::distance(double theta, double fi, Vertice_data V){
    return V.r * acos(sin(theta) * sin(V.theta) + cos(theta) * cos(V.theta) * cos(fi - V.fi));
}

double Extrapolation::triangle_counter(vector<Vertice_data> vect, double theta, double fi){
    double U=0;
    double dist1, dist2, dist3, tdist;
    Vertice_data V1, V2, V3, tV;
    
    dist1 = distance(theta, fi, vect[0]);
    V1 = vect[0];
    dist2 = distance(theta, fi, vect[1]);
    V2 = vect[1];
    if(dist1 > dist2){
        tdist = dist1;
        dist1 = dist2;
        dist2 = tdist;
        tV = V1;
        V1 = V2;
        V2 = tV;
    }
    dist3 = distance(theta, fi, vect[2]);
    V3 = vect[2];
    if(dist3 < dist1){
        tdist = dist1;
        dist1 = dist3;
        dist3 = tdist;
        tdist = dist2;
        dist2 = dist3;
        dist3 = tdist;
        tV = V1;
        V1 = V3;
        V3 = tV;
        tV = V2;
        V2 = V3;
        V3 = tV;
    }
    else if(dist3 > dist1 && dist3 < dist2){
        tdist = dist2;
        dist2 = dist3;
        dist3 = tdist;
        tV = V2;
        V2 = V3;
        V3 = tV;
    }
    
    for(int i=3; i<vect.size(); i++){
        double cur_dist = distance(theta, fi, vect[i]);
        if(cur_dist < dist1){
            tdist = dist1;
            dist1 = cur_dist;
            dist3 = tdist;
            tdist = dist2;
            dist2 = dist3;
            dist3 = tdist;
            tV = V1;
            V1 = vect[i];
            V3 = tV;
            tV = V2;
            V2 = V3;
            V3 = tV;
        }
        else if (cur_dist >= dist1 && cur_dist < dist2){
            tdist = dist2;
            dist2 = cur_dist;
            dist3 = tdist;
            tV = V2;
            V2 = vect[i];
            V3 = tV;
        }
        else if (cur_dist >= dist2 && cur_dist < dist3){
            dist3 = cur_dist;
            V3 = vect[i];
        }
    }
    if(dist1 == 0){
        U = V1.U;
    }
    else{
        U = (V1.U/(dist1 * dist1 * dist1) + V2.U/(dist2 * dist2 * dist2) +
             V3.U/(dist3 * dist3 * dist3)) / (1/(dist1 * dist1 * dist1) +
                                              1/(dist2 * dist2 * dist2) +
                                              1/(dist3 * dist3 * dist3));
    }
    //U = (dist1 * V1.U + dist2 * V2.U + dist3 * V3.U) / (dist1 + dist2 + dist3);
    return U;
}

double Extrapolation::extrapolator(int deg, double r, double theta, double fi){
    double U=0;
    
    /*
    vector<string> name = this->name_searcher(r);
    vector<vertice_data> lower_r = this->LoadFromFile(name[0] + to_string(deg) + ".txt");
    vector<vertice_data> higher_r = this->LoadFromFile(name[1] + to_string(deg) + ".txt");
    
    U = (this->triangle_counter(lower_r, theta, fi) * (r - lower_r[0].r) +
         this->triangle_counter(higher_r, theta, fi) * (higher_r[0].r - r)) /
        (higher_r[0].r - lower_r[0].r);
     */
    return U;
}

