#include <iostream>
#include <functional>
#include "PotentialMesh.h"

#define N 128
#define K 1
#define WG 0.9
#define TOL 1e-9

#define DELTAXY 0.01
#define LPOTENTIAL 0.1
#define RPOTENTIAL -0.1
#define UPOTENTIAL 0
#define DPOTENTIAL 0

#define X1 0.4
#define X2 0.88
#define Y1 0.64
#define Y2 0.64
#define Q 80
#define R  0.051

using namespace std;

int main() {
    BorderPotential borderPotential = BorderPotential(LPOTENTIAL, RPOTENTIAL, UPOTENTIAL, DPOTENTIAL);
    PotentialMesh potentialMesh = PotentialMesh(N, K, WG, TOL, DELTAXY, DELTAXY, &borderPotential, [](double x, double y)->double {
        x= x * DELTAXY;
        y= y * DELTAXY;
        if((x - X1)*(x - X1)  + (y - Y1)*(y - Y1) <= R*R)
            return Q;
        if((x - X2)*(x - X2)  + (y - Y2)*(y - Y2) <= R*R)
            return -Q;
        return 0;
    });
    potentialMesh.RelaxationG();
    //potentialMesh.PrintPotentialMesh();


    return 0;
}