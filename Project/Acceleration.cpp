//
//  Acceleration.cpp
//  Project
//
//  Created by Георгий on 08.03.2020.
//  Copyright © 2020 Георгий. All rights reserved.
//

#include "Acceleration.hpp"

void Acceleration::dSet(double dR, double dTh, double dFi){
    this->dR = dR;
    this->dTh = dTh;
    this->dFi = dFi;
}

void Acceleration::countLoader(){
    this->pc = PotentialCounter();
    this->pc.length = 100;
    this->pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
}

void Acceleration::extrLoader(){
    this->ex = Extrapolation();
    this->ex.loader(6, 6389100);
}

double Acceleration::U0(double r){
    return f_c * M_c / r;
}

double Acceleration::U2(double r, double theta){
    return -(f_c * M_c * J2 * r0 * r0 * (3* sin(theta) * sin(theta) - 1)) / (2 * r * r * r);
}

double Acceleration::accR_C(double r, double theta, double fi){
    double Up = U0(r + dR) + U2(r + dR, theta) + this->pc.potential(r + dR, theta, fi);
    double Ul = U0(r)      + U2(r, theta)      + this->pc.potential(r,      theta, fi);
    return (Up - Ul) / dR;
}

double Acceleration::accR_E(double r, double theta, double fi){
    double Up = U0(r + dR) + U2(r + dR, theta) + this->ex.extrapolator(r + dR, theta, fi);
    double Ul = U0(r)      + U2(r, theta)      + this->ex.extrapolator(r,      theta, fi);
    return (Up - Ul) / dR;
}

double Acceleration::accTh_C(double r, double theta, double fi){
    double Up = U2(r, theta + dTh) + this->pc.potential(r, theta + dTh, fi);
    double Ul = U2(r, theta      ) + this->pc.potential(r, theta      , fi);
    return (Up - Ul) / (r*dTh);
}

double Acceleration::accTh_E(double r, double theta, double fi){
    double Up = U2(r, theta + dTh) + this->ex.extrapolator(r, theta + dTh, fi);
    double Ul = U2(r, theta      ) + this->ex.extrapolator(r, theta      , fi);
    return (Up - Ul) / (r*dTh);
}

double Acceleration::accFi_C(double r, double theta, double fi){
    double Up = this->pc.potential(r, theta, fi + dFi);
    double Ul = this->pc.potential(r, theta, fi      );
    return (Up - Ul) / (r*dFi);
}

double Acceleration::accFi_E(double r, double theta, double fi){
    double Up = this->ex.extrapolator(r, theta, fi + dFi);
    double Ul = this->ex.extrapolator(r, theta, fi);
    return (Up - Ul) / (r*dFi);
}

double Acceleration::g_C(double r, double theta, double fi){
    double aR = accR_C(r, theta, fi);
    double aTh = accTh_C(r, theta, fi);
    double aFi = accFi_C(r, theta, fi);
    return sqrt(aR*aR + aTh*aTh + aFi*aFi);
}

double Acceleration::g_E(double r, double theta, double fi){
    double aR = accR_E(r, theta, fi);
    double aTh = accTh_E(r, theta, fi);
    double aFi = accFi_E(r, theta, fi);
    return sqrt(aR*aR + aTh*aTh + aFi*aFi);
}
