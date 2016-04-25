//
// Created by Kamil on 24.04.2016.
//
#ifndef IMN2_MESH_H
#define IMN2_MESH_H
#include <iostream>
#include <functional>
#include <vector>
#include <cmath>


using namespace std;

class BorderPotential {
private:
    double _potenatial[4] = {0,0,0,0};
public:
    BorderPotential(double lPotential, double rPotential,  double uPotential, double dPotential) {
        SetPotential(lPotential, rPotential,  uPotential, dPotential);
    }
    void SetPotential(double lPotential, double rPotential,  double uPotential, double dPotential) {
        _potenatial[0] = lPotential;
        _potenatial[1] = rPotential;
        _potenatial[2] = uPotential;
        _potenatial[3] = dPotential;
    }
    double GetLPotential() {
        return _potenatial[0];
    }
    double GetRPotential() {
        return _potenatial[1];
    }
    double GetUPotential() {
        return _potenatial[2];
    }
    double GetDPotential() {
        return _potenatial[3];
    }

};

class PotentialMesh {
private:
    int _n, _k, _iter;
    double _dx, _dy, _wg, _TOL, _SNow, _SPrev, _STOL;
    BorderPotential * _borderPotential;
    function < double( int x, int y) > _densityFn;
    vector<vector<double > >  _densityMesh;
    vector<vector<double > >  _densityMeshPrevious;
    void FillBorder() {
        for(int i=1;i<_n-1;i++)
        {
            _densityMeshPrevious[i][0] = _borderPotential->GetLPotential();
            _densityMeshPrevious[i][_n-1] = _borderPotential->GetRPotential();
        }
        for(int i=0;i<_n;i++)
        {
            _densityMeshPrevious[0][i] = _borderPotential->GetUPotential();
            _densityMeshPrevious[_n-1][i] = _borderPotential->GetDPotential();
        }
    }
public:
    PotentialMesh(int n, int k, double wg, double TOL,  double dx, double dy, BorderPotential * borderPotential,  function < double( int x, int y) > densityFn)  {
        InitMesh(n, k, wg, TOL, dx, dy, borderPotential, densityFn);
    }
    void InitMesh(int n, int k, double wg, double TOL, double dx, double dy, BorderPotential * borderPotential,  function < double( int x, int y) > densityFn) {
        _iter = 0;
        _n = n;
        _k = k;
        _wg = wg;
        _TOL = TOL;
        _STOL = _TOL;
        _dx = dx;
        _dy = dy;
        _borderPotential = borderPotential;
        _densityFn = densityFn;
        _densityMesh.resize(n);
        _densityMeshPrevious.resize(n);
        for(auto &i: _densityMesh) {
            i.resize(n);
        }
        for(auto &i: _densityMeshPrevious) {
            i.resize(n);
        }
        FillBorder();
    }

    void ResetState() {
        _iter = 0;
        for (auto vect: _densityMesh) {
            fill(vect.begin(), vect.end(), 0);
        }
        for (auto vect: _densityMeshPrevious) {
            fill(vect.begin(), vect.end(), 0);
        }
        FillBorder();
    }
    void SetWG(double wg){
        _wg = wg;
    }
    void PrintPotentialMesh() {
        for (auto y: _densityMesh) {
            for (auto x: y) {
                cout << x << ", ";
            }
            cout << "\b\b\n";
        }
    }
    void RelaxationG() {
        while (_STOL >= _TOL) {
            for (int j = _k; j < _n - _k; j++) {
                for (int i = _k; i < _n - _k; i++) {
                    _densityMesh[j][i] = (_densityMeshPrevious[j][i + _k]
                                          + _densityMeshPrevious[j][i - _k]
                                          + _densityMeshPrevious[j + _k][i]
                                          + _densityMeshPrevious[j - _k][i]
                                          + (_k * _dx * _k * _dx) * _densityFn(i, j))
                                         / (4);
                }
            }
            for (int j = _k; j < _n - _k; j++) {
                for (int i = _k; i < _n - _k; i++) {
                    _densityMesh[j][i] = (1 - _wg) * _densityMeshPrevious[j][i] + _wg * _densityMesh[j][i];
                }
            }
            for (int i = _k; i < _n - _k; i++) {
                for (int j = _k; j < _n - _k; j++) {
                  _SNow = pow((_k * _dx), 2)
                          * (0.5 * (
                                    pow(((_densityMesh[j][i+_k] - _densityMesh[j][i])/(_k*_dx)), 2)
                                  + pow(((_densityMesh[j+_k][i] - _densityMesh[j][i])/(_k*_dx)), 2)
                                   ) - _densityFn(i, j)*_densityMesh[j][i]);
                }
            }
            if(_iter > 0) {
                _STOL = abs((_SNow - _SPrev)/_SPrev);
            }
            cout << _iter <<  " "  << _SNow << "\n";
            _SPrev = _SNow;
            _densityMeshPrevious = _densityMesh;
            _iter++;
        }


    }


protected:

};

#endif //IMN2_MESH_H
