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
    if(this->x == 0 && this->y == 0)
        this->fi = 0;
    else
        this->fi = atan2(this->y, this->x);
      /*
        if(this->x == 0 && this->y < 0){
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
       */
}

bool VerticeWithCompared::operator<(const VerticeWithCompared &rht) const{
    return Triangulation::distCounter(*this, *compared) < Triangulation::distCounter(rht, *compared);
}

Triangle::Triangle(){
    index = 0;
    fatherInd = 0;
    V.clear();
}

Triangle::Triangle(int index, VerticeWithCompared* V1, VerticeWithCompared* V2, VerticeWithCompared* V3, int fatherIndex){
    fatherInd = fatherIndex;
    this->index = index;
    V.push_back(V1->index);
    V.push_back(V2->index);
    V.push_back(V3->index);
}

Triangle::Triangle(int index, int V1, int V2,int V3, int fatherIndex){
    fatherInd = fatherIndex;
    this->index = index;
    V.push_back(V1);
    V.push_back(V2);
    V.push_back(V3);
}

bool Triangle::operator==(const Triangle &tr) const{
    bool eq = false;
    for(int i=0; i<this->V.size(); i++){
        for(int j=0; j<tr.V.size(); j++){
            if(this->V.at(i) == tr.V.at(j)){
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
        Triangle tr = Triangle(i, &vert_arr.at(0), &vert_arr.at(i+1), &vert_arr.at(i+2));
        this->tr_arr.push_back(tr);
    }
    //5
    Triangle tr5 = Triangle(4, &vert_arr.at(0), &vert_arr.at(5), &vert_arr.at(1));
    this->tr_arr.push_back(tr5);
    //6,8,10,12
    for(int i=1; i<5; i++){
        Triangle tr = Triangle(i+4, &vert_arr.at(i), &vert_arr.at(i+1), &vert_arr.at(i+5));
        this->tr_arr.push_back(tr);
    }
    //7,9,11,13
    for(int i=2; i<6; i++){
        Triangle tr = Triangle(i+7, &vert_arr.at(i), &vert_arr.at(i+4), &vert_arr.at(i+5));
        this->tr_arr.push_back(tr);
    }
    //14-15
    for(int i = 5; i<=6; i++){
        Triangle tr = Triangle(i+8, &vert_arr.at(i), &vert_arr.at(10), &vert_arr.at(1));
        this->tr_arr.push_back(tr);
    }
    //16-19
    for(int i=6; i<10; i++){
        Triangle tr = Triangle(i+9, &vert_arr.at(11), &vert_arr.at(i), &vert_arr.at(i+1));
        this->tr_arr.push_back(tr);
       }
    //20
    Triangle tr20 = Triangle(19, &vert_arr.at(11), &vert_arr.at(10), &vert_arr.at(6));
    this->tr_arr.push_back(tr20);
    
    this->globalIndex=19;
}

void Triangulation::zeroMeshCreator(double r)
{
    this->setRad(r);
    const float H_ANGLE = pi / 180 * 72;
    const float V_ANGLE = atan(1.0 / 2);
    this->vert_arr.clear();
    VerticeWithCompared tre = {};
    this->vert_arr.resize(12, tre);
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
        
        this->vert_arr[i].x = xy * cos(hAngle1);
        
        this->vert_arr[i+5].x = xy * cos(hAngle2);
        
        this->vert_arr[i].y = xy * sin(hAngle1);
        
        this->vert_arr[i+5].y = xy * sin(hAngle2);
        
        this->vert_arr[i].z = z;
        
        this->vert_arr[i+5].z = -z;
        
        this->vert_arr.at(i).toSph();
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
    this->zero_triangles();
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
            if(V.x == tmpVertices.at(i).x &&
               V.y == tmpVertices.at(i).y &&
               V.z == tmpVertices.at(i).z){
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
        closestVertices.push_back(tmpVertices.at(k));
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
        if(this->vert_arr.at(i).x == V.x &&
           this->vert_arr.at(i).y == V.y &&
           this->vert_arr.at(i).z == V.z ){
            this->localVIndex = i;
            return false;
        }
    }
    return true;
}

void Triangulation::mesher(double r, int degree){
    this->zeroMeshCreator(r);
    vector<vector<Triangle>> all_triangles;
    all_triangles.push_back(tr_arr);
    tr_arr.clear();
    for(int iter = 1; iter <= degree; iter++)
    {
        all_triangles.push_back(vector<Triangle>());
        for (int tr_num = 0; tr_num < all_triangles.at(iter-1).size(); tr_num++)
        {
            Triangle& t = all_triangles.at(iter-1).at(tr_num);
            VerticeWithCompared Vs[3] = {};
            int Vptrs[3] = {};
            for (int idx = 0; idx < 3; idx ++) {
                Vs[idx] = verticeCreator(this->vert_arr.at(t.V.at(idx % 3)), this->vert_arr.at(t.V.at((idx + 1) % 3)));
                if(vertDetector(Vs[idx])){
                    Vs[idx].index = this->vIndex;
                    this->vIndex += 1;
                    this->vert_arr.push_back(Vs[idx]);
                    Vptrs[idx] = Vs[idx].index;
                }
                else{
                    Vptrs[idx] = localVIndex;
                }
            }
            
            Triangle tr1 = Triangle(this->globalIndex + 1, t.V[0], Vptrs[0], Vptrs[2], t.index);
           
            Triangle tr2 = Triangle(this->globalIndex + 2, t.V[1], Vptrs[0], Vptrs[1], t.index);
         
            Triangle tr3 = Triangle(this->globalIndex + 3, t.V[2], Vptrs[1], Vptrs[2], t.index);
            
            Triangle tr4 = Triangle(this->globalIndex + 4, Vptrs[0], Vptrs[1], Vptrs[2], t.index);
            
            t.childInd.clear();
            t.childInd.push_back(tr1.index);
            t.childInd.push_back(tr2.index);
            t.childInd.push_back(tr3.index);
            t.childInd.push_back(tr4.index);
            all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr1);
            all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr2);
            all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr3);
            all_triangles.at(iter).insert(all_triangles.at(iter).end(), tr4);
            this->globalIndex = this->globalIndex + 4;
        }
        tr_arr.insert(tr_arr.end(), all_triangles.at(iter-1).begin(), all_triangles.at(iter-1).end());
        if(iter == degree){
            tr_arr.insert(tr_arr.end(), all_triangles.at(iter).begin(), all_triangles.at(iter).end());
        }
    }
}

void Triangulation::mesherDot(double r, int degree){
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
                        this->vert_arr.push_back(V1);
                    }
                    if(vertDetector(V2)){
                        this->vert_arr.push_back(V2);
                    }
                    if(vertDetector(V3)){
                        this->vert_arr.push_back(V3);
                    }
                }
            }
        }
    }
}

void Triangulation::map(int degree, int polynomDegree, int maxRad, int step){
    PotentialCounter pc = {};
    pc.length = polynomDegree;
    pc.LoadFromFile("/Users/georgij/С++/Project/Project/data.txt");
    
    string deg  = to_string(degree);
    string polD = to_string(polynomDegree);
    string st   = to_string(step);
    ofstream infileV;
    infileV.open("v_"+deg+"_"+polD+"_"+st+".txt", std::ios_base::app);
    ofstream infile("Triangles_" + deg + ".txt");
    cout<<"Map in progress:"<<endl;
    time_t start, end;
    time(&start);
    
    cout<<"Counting triangles with degree "<<degree<<endl;
    this->mesher(r0, degree);
    for (int j = 0; j < this->tr_arr.size(); j++)
    {
        infile<<this->tr_arr[j].index<<" "<<this->tr_arr[j].fatherInd<<" ";
        for(int k=0; k<4; k++){
            infile<<this->tr_arr[j].childInd[k]<<" ";
        }
        for(int i=0; i<3; i++){
            infile<<this->vert_arr.at(this->tr_arr[j].V[i]).index<<" ";
        }
        infile<<"\n";
    }
    
    infile<<endl;
    infile.close();
    
    cout<<"Meshing..."<<endl;
     
    for(int r=r0; r<=maxRad; r+= step){
        time_t oneMapSt, oneMapF;
        time(&oneMapSt);
        this->mesher(r, degree);
        cout<<"Counting potentials on rad "<<r<<endl;

        for(int i=0; i<this->vert_arr.size(); i++){
            this->vert_arr[i].U = pc.potential(this->vert_arr[i].r, this->vert_arr[i].theta, this->vert_arr[i].fi);
            infileV<<this->vert_arr[i].index <<" "<<
                     this->vert_arr[i].U     <<" "<<
                     this->vert_arr[i].r     <<" "<<
                     this->vert_arr[i].theta <<" "<<
                     this->vert_arr[i].fi    <<" "<<
                     this->vert_arr[i].x     <<" "<<
                     this->vert_arr[i].y     <<" "<<
                     this->vert_arr[i].z     <<"\n";
            
        }
        
        time(&oneMapF);
        cout<<"Map with rad "<<r<<" is done in "<<difftime(oneMapF, oneMapSt)/60.0<<" min"<<endl;
    }
    infileV.close();
    time(&end);
    cout<<"Global map is done in "<<difftime(end, start)/60.0<<" min"<<endl;
}

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
            this->vert_arr[i].U = pc.potential(this->vert_arr[i].r, this->vert_arr[i].theta, this->vert_arr[i].fi);
            infile<<this->vert_arr[i].index <<" "<<
                    this->vert_arr[i].U     <<" "<<
                    this->vert_arr[i].r     <<" "<<
                    this->vert_arr[i].theta <<" "<<
                    this->vert_arr[i].fi    <<"\n";
        }
        time(&oneMapF);
        cout<<"Map with rad "<<r<<" is done in "<<difftime(oneMapF, oneMapSt)/60.0<<" min"<<endl;
    }
    infile.close();
    time(&end);
    cout<<"Global map is done in "<<difftime(end, start)/60.0<<" min"<<endl;
}

