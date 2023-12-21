#pragma once
#include <vector>
#include <string>
using namespace std;

#define pid pair<int, double>
#define PointSet vector<vector<double>>

PointSet readFile(const string& fname);

void writeFile(const PointSet& p, const string& fname, const vector<int>& indexW, const vector<int>& indexU);

vector<double> f1(const vector<double>& w, const vector<double>& u, const double& tau);

vector<double> f2(const vector<double>& w, const vector<double>& u, const double& tau);

vector<double> f3(const vector<double>& w, const vector<double>& u, const double& tau);

PointSet computeNewSet(const PointSet& W, const PointSet& U, const double& tau, vector<int>& indexW, vector<int>& indexU);

int Round(const double& x);

vector<vector<int>> hyperCube(const int& n_n, const int& n);

PointSet reduceSet(const PointSet& W, const double& h, const vector<vector<int>>& dx, vector<int>& indexW, vector<int>& indexU);

bool notInPhaseConstraint(const vector<vector<double>>& pCon, const double& px, const double& py);

double sqr_dist(const vector<double>& a, const vector<double>& b);