//
//  Triangulation.cpp
//  Project
//
//  Created by Георгий on 09.12.2019.
//  Copyright © 2019 Георгий. All rights reserved.
//

#include "Triangulation.hpp"

void Triangulation::setRad(double r)
{
    this->r=r;
}

void VerticeWithCompared::toSph()
{
    this->r = sqrt(this->x * this->x +
                   this->y * this->y +
                   this->z * this->z);
    
    this->theta = acos(this->z/sqrt(this->x * this->x +
                                    this->y * this->y +
                                    this->z * this->z));
    if(this->x == 0 && this->y == 0){
        this->fi = 0;
    }
    else if(this->x == 0 && this->y < 0){
        this->fi = 3*pi/2;
    }
    else if(this->x == 0 && this->y > 0){
        this->fi = pi/2;
    }
    else if(this->x > 0 && this->y == 0){
        this->fi = 0;
    }
    else if(this->x < 0 && this->y == 0){
        this->fi = pi;
    }
    else if(this->x > 0 && this->y > 0){
        this->fi = atan(this->y/this->x);
    }
    else if(this->x < 0){
        this->fi = pi+atan(this->y/this->x);
    }
    else if(this->x > 0 && this->y < 0){
        this->fi = 2*pi+atan(this->y/this->x);
    }
}

bool VerticeWithCompared::operator<(const VerticeWithCompared &rht) const{
    return Triangulation::distCounter(*this, *compared) < Triangulation::distCounter(rht, *compared);
}

Triangle::Triangle(){
    index = 0;
    fatherInd = 0;
    V.clear();
}

Triangle::Triangle(int index, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3){
    fatherInd = 0;
    this->index = index;
    V.push_back(V1);
    V.push_back(V2);
    V.push_back(V3);
}

bool Triangle::operator==(const Triangle &tr) const{
    bool eq = false;
    for(int i=0; i<this->V.size(); i++){
        for(int j=0; j<tr.V.size(); j++){
            if(this->V[i]->index == tr.V[j]->index){
                eq = true;
                break;
            }
            else{
                eq = false;
            }
        }
        if(eq == false){
            break;
        }
    }
    return eq;
}

double Triangulation::distCounter(VerticeWithCompared V1, VerticeWithCompared V2){
    double dist = sqrt((V1.x - V2.x) * (V1.x - V2.x) +
                       (V1.y - V2.y) * (V1.y - V2.y) +
                       (V1.z - V2.z) * (V1.z - V2.z));
    return dist;
}

void Triangulation::zero_triangles(){
    this->tr_arr.clear();
    //1-4
    for(int i=0; i<4; i++){
        Triangle tr = Triangle(i, &vert_arr[0], &vert_arr[i+1], &vert_arr[i+2]);
        this->tr_arr.push_back(tr);
    }
    //5
    Triangle tr5 = Triangle(4, &vert_arr[0], &vert_arr[5], &vert_arr[1]);
    this->tr_arr.push_back(tr5);
    //6,8,10,12
    for(int i=1; i<5; i++){
        Triangle tr = Triangle(i+4, &vert_arr[i], &vert_arr[i+1], &vert_arr[i+5]);
        this->tr_arr.push_back(tr);
    }
    //7,9,11,13
    for(int i=2; i<6; i++){
        Triangle tr = Triangle(i+7, &vert_arr[i], &vert_arr[i+4], &vert_arr[i+5]);
        this->tr_arr.push_back(tr);
    }
    //14-15
    for(int i = 5; i<=6; i++){
        Triangle tr = Triangle(i+8, &vert_arr[i], &vert_arr[10], &vert_arr[1]);
        this->tr_arr.push_back(tr);
    }
    //16-19
    for(int i=6; i<10; i++){
        Triangle tr = Triangle(i+9, &vert_arr[11], &vert_arr[i], &vert_arr[i+1]);
        this->tr_arr.push_back(tr);
       }
    //20
    Triangle tr20 = Triangle(19, &vert_arr[11], &vert_arr[10], &vert_arr[6]);
    this->tr_arr.push_back(tr20);
    
    this->globalIndex=19;
}

void Triangulation::zeroMeshCreator(double r)
{
    this->setRad(r);
    const float H_ANGLE = pi / 180 * 72;
    const float V_ANGLE = atanf(1.0f / 2);
    this->vert_arr.clear();
    this->vert_arr.resize(12);
    float z, xy;
    float hAngle1 = -pi / 2 - H_ANGLE / 2;
    float hAngle2 = -pi / 2;
    
    this->vert_arr[0].x = 0;
    this->vert_arr[0].y = 0;
    this->vert_arr[0].z = this->r;
    this->vert_arr[0].toSph();
    this->vert_arr[0].index = 0;
    
    for(int i = 1; i <= 5; ++i)
    {
        z = this->r * sin(V_ANGLE);
        xy = this->r * cos(V_ANGLE);
        
        this->vert_arr[i].x = xy * cosf(hAngle1);
        
        this->vert_arr[i+5].x = xy * cosf(hAngle2);
        
        this->vert_arr[i].y = xy * sinf(hAngle1);
        
        this->vert_arr[i+5].y = xy * sinf(hAngle2);
        
        this->vert_arr[i].z = z;
        
        this->vert_arr[i+5].z = -z;
        
        this->vert_arr[i].toSph();
        this->vert_arr[i+5].toSph();
        this->vert_arr[i].index = i;
        this->vert_arr[i+5].index = i+5;
        
        hAngle1 += H_ANGLE;
        hAngle2 += H_ANGLE;
    }
    
    this->vert_arr[11].x = 0;
    this->vert_arr[11].y = 0;
    this->vert_arr[11].z = -1.0*this->r;
    this->vert_arr[11].toSph();
    this->vert_arr[11].index = 11;
    this->vIndex=12;
    
    //PotentialCounter pc = {};
    //pc.length = this->polynomDeg;
    //pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    //for(int j=0; j<12; j++){
    //    this->vert_arr[j].U = pc.potential(this->vert_arr[j].r, this->vert_arr[j].theta - pi/2, this->vert_arr[j].fi - pi);
    //}
    
    //this->zero_triangles();
}

vector<VerticeWithCompared> Triangulation::closestVerticesFinder(VerticeWithCompared V, vector<VerticeWithCompared> tmpVertices, int mode){
    
    for (auto &v: tmpVertices)
        v.compared = &V;
    vector<VerticeWithCompared> closestVertices;
    
    int verticesNum = 0;
    if(mode == 1){
        verticesNum = 2;
    }
    else{
        for(int i = 0; i < 12; i++){
            
            if(V.x == tmpVertices[i].x &&
               V.y == tmpVertices[i].y &&
               V.z == tmpVertices[i].z){
                verticesNum = 5;
                break;
            }
            else{
                verticesNum = 6;
            }
        }
    }
    
    sort(tmpVertices.begin(),tmpVertices.end());
    
    for(int k = 1; k <= verticesNum; k++){
        closestVertices.push_back(tmpVertices[k]);
    }
    
    return closestVertices;
}

VerticeWithCompared Triangulation::verticeCreator(VerticeWithCompared V1, VerticeWithCompared V2){
    VerticeWithCompared newVertice;
    
    newVertice.x = (V1.x + V2.x)/2.0;
    newVertice.y = (V1.y + V2.y)/2.0;
    newVertice.z = (V1.z + V2.z)/2.0;
    
    double norm = this->r / sqrt(newVertice.x * newVertice.x +
                                 newVertice.y * newVertice.y +
                                 newVertice.z * newVertice.z);
    newVertice.x *= norm;
    newVertice.y *= norm;
    newVertice.z *= norm;
    
    newVertice.toSph();
    
    return newVertice;
}

bool Triangulation::vertDetector(VerticeWithCompared V){
    for(int i = 0; i < this->vert_arr.size(); i++){
        if(V.x == this->vert_arr[i].x &&
           V.y == this->vert_arr[i].y &&
           V.z == this->vert_arr[i].z){
            //this->localVIndex = i;
            return false;
        }
    }
    return true;
}

void Triangulation::triangledetector(Triangle* tr, vector<Triangle> trArr){
    for(int k=0; k<trArr.size(); k++){
        if(*tr == trArr[k]){
            this->localTrIndex = k;
            vector<int> indexes = {tr->V[0]->index, tr->V[1]->index, tr->V[2]->index};
            vector<VerticeWithCompared*> tmpV = trArr[k].V;
            trArr[k].V.clear();
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    if(tmpV[j]->index == indexes[i]){
                        trArr[k].V.push_back(tmpV[j]);
                        break;
                    }
                }
            }
            break;
        }
    }
}

vector<Triangle> Triangulation::trCreator(Triangle* fthrTr, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3){
    vector<Triangle> newTr;
    Triangle tr1 = Triangle(this->globalIndex + 1, fthrTr->V[0], V1, V2);
    tr1.fatherInd = fthrTr->index;
    fthrTr->childInd.push_back(tr1.index);
    newTr.push_back(tr1);
    
    Triangle tr2 = Triangle(this->globalIndex + 2, V1, fthrTr->V[1], V3);
    tr2.fatherInd = fthrTr->index;
    fthrTr->childInd.push_back(tr2.index);
    newTr.push_back(tr2);
    
    Triangle tr3 = Triangle(this->globalIndex + 3, V2, V3, fthrTr->V[2]);
    tr3.fatherInd = fthrTr->index;
    fthrTr->childInd.push_back(tr3.index);
    newTr.push_back(tr3);
    
    Triangle tr4 = Triangle(this->globalIndex + 4, V1, V2, V3);
    tr4.fatherInd = fthrTr->index;
    fthrTr->childInd.push_back(tr4.index);
    newTr.push_back(tr4);
    
    fthrTr->isDone = true;
    this->globalIndex = this->globalIndex + 4;
    return newTr;
}
/*
void Triangulation::mesher(double r, int degree){
    this->zeroMeshCreator(r);
    for(int iter = 1; iter <= degree; iter++)
    {
        vector<VerticeWithCompared> tmpVertices = this->vert_arr;
        //vector<Triangle> tmpTriangles = this->tr_arr;
        
        for(int j = 0; j < tmpVertices.size(); j++){
        vector<VerticeWithCompared> closestVertices = closestVerticesFinder(tmpVertices[j], tmpVertices, 0);
            for(int l = 0; l < closestVertices.size(); l++){
                vector<VerticeWithCompared> thirdVertice = closestVerticesFinder(closestVertices[l], closestVertices, 1);
                for(int i = 0; i < thirdVertice.size(); i++){
                    Triangle tmpTr = Triangle(0, &tmpVertices[j], &closestVertices[l], &thirdVertice[i]);
                    
                    VerticeWithCompared V1 = verticeCreator(tmpVertices[j], closestVertices[l]);
                    VerticeWithCompared* V1Ind;
                    VerticeWithCompared V2 = verticeCreator(tmpVertices[j], thirdVertice[i]);
                    VerticeWithCompared* V2Ind;
                    VerticeWithCompared V3 = verticeCreator(closestVertices[l], thirdVertice[i]);
                    VerticeWithCompared* V3Ind;
                    
                    if(vertDetector(V1)){
                        V1.index = this->vIndex;
                        this->vIndex += 1;
                        this->vert_arr.push_back(V1);
                        V1Ind = &this->vert_arr[V1.index];
                    }
                    else{
                        V1Ind = &this->vert_arr[localVIndex];
                    }
                    
                    if(vertDetector(V2)){
                        V2.index = this->vIndex;
                        this->vIndex += 1;
                        this->vert_arr.push_back(V2);
                        V2Ind = &this->vert_arr[V2.index];
                    }
                    else{
                        V2Ind = &this->vert_arr[localVIndex];
                    }
                    
                    if(vertDetector(V3)){
                        V3.index = this->vIndex;
                        this->vIndex += 1;
                        this->vert_arr.push_back(V3);
                        V3Ind = &this->vert_arr[V3.index];
                    }
                    else{
                        V3Ind = &this->vert_arr[localVIndex];
                    }
                    
                    //Triangle* fthrTr;
                    triangledetector(&tmpTr, this->tr_arr);
                    //if(fthrTr->V.size() != 3) { cout<<"no such triangle\n"; }
                    
                    if(this->tr_arr[localTrIndex].isDone == false){
                        vector<Triangle> newTriangles = trCreator(&this->tr_arr[localTrIndex], V1Ind, V2Ind, V3Ind);
                        //tmpTriangles[this->localTrIndex] = this->tr_arr[localTrIndex];
                        //this->tr_arr[this->localTrIndex] = fthrTr;
                        for(int h=0; h<4; h++){
                            this->tr_arr.push_back(newTriangles[h]);
                        }
                    }
                }
            }
        }
    }
}
*/




void Triangulation::mesherDot(double r, int degree){
    //PotentialCounter pc = {};
    //pc.length = polynomDeg;
    //pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    this->zeroMeshCreator(r);
    for(int iter = 1; iter <= degree; iter++)
    {
        vector<VerticeWithCompared> tmpVertices = this->vert_arr;
        for(int j = 0; j < tmpVertices.size(); j++){
            vector<VerticeWithCompared> closestVertices = closestVerticesFinder(tmpVertices[j], tmpVertices, 0);
            for(int l = 0; l < closestVertices.size(); l++){
                vector<VerticeWithCompared> thirdVertice = closestVerticesFinder(closestVertices[l], closestVertices, 1);
                for(int i = 0; i < thirdVertice.size(); i++){
                    VerticeWithCompared V1 = verticeCreator(tmpVertices[j], closestVertices[l]);
                    VerticeWithCompared V2 = verticeCreator(tmpVertices[j], thirdVertice[i]);
                    VerticeWithCompared V3 = verticeCreator(closestVertices[l], thirdVertice[i]);
                    if(vertDetector(V1)){
                        //V1.U = pc.potential(V1.r, V1.theta - pi/2, V1.fi - pi);
                        this->vert_arr.push_back(V1);
                    }
                    if(vertDetector(V2)){
                        //V2.U = pc.potential(V2.r, V2.theta - pi/2, V2.fi - pi);
                        this->vert_arr.push_back(V2);
                    }
                    if(vertDetector(V3)){
                        //V3.U = pc.potential(V3.r, V3.theta - pi/2, V3.fi - pi);
                        this->vert_arr.push_back(V3);
                    }
                }
            }
        }
    }
}

/*
void Triangulation::map(int degree, double r){
    PotentialCounter pc = {};
    pc.length = 100;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    time_t start, end;
    time(&start);
    cout<<"work on map with degree "<<degree<<" and r = "<<r<<" started"<<endl;
    this->mesher(r, degree);
    this->vert_arr.clear();
    cout<<"mesh is ready, counting potentials"<<endl;
    string rad = to_string(int(r));
    string deg = to_string(degree);
    ofstream infile(rad+"_"+deg+".txt");
    for (int j = 0; j < this->tr_arr.size(); j++)
    {
        for(int i=0; i<3; i++){
            this->tr_arr[j].V[i]->U = pc.potential(this->tr_arr[j].V[i]->r, this->tr_arr[j].V[i]->theta - pi/2,
                                                   this->tr_arr[j].V[i]->fi - pi);
            infile<<this->tr_arr[j].V[i]->U<<" "<<this->tr_arr[j].V[i]->r<<" "<<this->tr_arr[j].V[i]->theta - pi/2<<" "<<
                    this->tr_arr[j].V[i]->fi - pi <<"\n";
        }
    }
    infile << endl;
    infile.close();
    time(&end);
    printf("Completed in %f min\n", difftime(end, start)/60.0);
}
*/
void Triangulation::globalMapDot(int degree, int polynomDegree, int maxRad, int step){
    PotentialCounter pc = {};
    pc.length = polynomDegree;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    string deg  = to_string(degree);
    string polD = to_string(polynomDegree);
    string maxR = to_string(maxRad);
    string st   = to_string(step);
    ofstream infile(deg+"_"+polD+"_"+maxR+"_"+st+".txt");
    cout<<"Map in progress:"<<endl;
    time_t start, end;
    time(&start);
    
    for(int r=r0; r<=maxRad; r+= step){
        time_t oneMapSt, oneMapF;
        time(&oneMapSt);
        this->mesherDot(r, degree);
        cout<<"Counting potentials on rad "<<r<<endl;
        for(int i=0; i<this->vert_arr.size(); i++){
            this->vert_arr[i].U = pc.potential(this->vert_arr[i].r, this->vert_arr[i].theta - pi/2, this->vert_arr[i].fi - pi);
            infile<<this->vert_arr[i].U<<" "<<
                    this->vert_arr[i].r<<" "<<
                    this->vert_arr[i].theta - pi/2<<" "<<
                    this->vert_arr[i].fi - pi<<"\n";
        }
        time(&oneMapF);
        cout<<"Map with rad "<<r<<" is done in "<<difftime(oneMapF, oneMapSt)/60.0<<" min"<<endl;
    }
    infile.close();
    time(&end);
    cout<<"Global map is done in "<<difftime(end, start)/60.0<<" min"<<endl;
}

