#include <fstream>
#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include "speed_challenge.h"


PointSet readFile(const string& fname) {

	ifstream ins(fname.c_str());

	PointSet points;
	if (ins.is_open()) {
		int n, m;
		ins >> n >> m;
		points = PointSet(n, vector<double>(m, 0));

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				ins >> points[i][j];
			}
		}
	}
	else {
		cout << "Open file with error." << endl;
	}

	ins.close();
	return points;
}

void writeFile(const PointSet& p, const string& fname, const vector<int>& indexW, const vector<int>& indexU) {

	ofstream outfile(fname.c_str());
	for (int i = 0; i < p.size(); ++i) {
		for (int j = 0; j < p[0].size(); ++j) {
			outfile << p[i][j] << " ";
		}
		outfile << indexW[i] << " " << indexU[i];
		outfile << "\n";
	}
	outfile.close();
}

vector<double> f1(const vector<double>& w, const vector<double>& u, const double& tau) {

	vector<double> result(w.size());
	result[0] = w[0] + tau * cos(w[2]);
	result[1] = w[1] + tau * sin(w[2]);
	result[2] = w[2] + tau * u[0];
	return result;
}

vector<double> f2(const vector<double>& w, const vector<double>& u, const double& tau) {

	vector<double> result(w.size());
	result[0] = w[0] + tau * w[1];
	result[1] = w[1] + tau * cos(w[4]);
	result[2] = w[2] + tau * w[3];
	result[3] = w[3] + tau * sin(w[4]);
	result[4] = w[4] + tau * u[0];
	return result;
}

vector<double> f3(const vector<double>& w, const vector<double>& u, const double& tau) {
	double m = 1;
	double g = 9.8;
	double L = 1;
	vector<double> result(w.size());
	result[0] = w[0] + tau * w[1];
	result[1] = w[1] + tau * (-m*g*sin(w[0])+u[0]) / (m * L);
	return result;
}


PointSet computeNewSet(const PointSet& W, const PointSet& U, const double& tau, vector<int>& indexW, vector<int>& indexU) {

	int n = W.size(), m = U.size();
	PointSet result(n * m);
	indexW.resize(n * m);
	indexU.resize(n * m);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			result[j + i * m] = f2(W[i], U[j], tau);
			indexW[j + i * m] = i;
			indexU[j + i * m] = j;
		}
	}

	return result;
}


int Round(const double& x) {

	int id = int(x);
	return x < 0 ? id - 1 : id;
}


vector<vector<int>> hyperCube(const int& n_n, const int& n) {
	vector<vector<int>> result(n_n, vector<int>(n));
	for (int i = 0; i < n_n; i++) {
		int tmp = i;
		for (int k = n - 1; k >= 0; k--) {
			int tmp_i = tmp;
			result[i][k] = tmp % 2;
			tmp = tmp_i / 2;
		}
	}
	return result;
}

PointSet reduceSet(const PointSet& W, const double& h, const vector<vector<int>>& dx, vector<int>& indexW, vector<int>& indexU) {
	PointSet resultSet;
	map<vector<int>, pid> netSet;
	vector<int> min_net_point;
	int n = W[0].size();
	int n_n = pow(2, n);
	vector<int> cur_net_point;
	pid value;
	for (int i = 0; i < W.size(); ++i) {
		min_net_point.clear();

		for (double c : W[i]) {
			min_net_point.push_back(Round(c / h));
		}

		for (int k = 0; k < n_n; ++k) {
			cur_net_point.clear();
			value.second = 0;
			for (int j = 0; j < n; ++j) {
				cur_net_point.push_back(min_net_point[j] + dx[k][j]);
				value.second += abs(cur_net_point.back() * h - W[i][j]);
			}
			value.first = i;

			map<vector<int>, pid>::iterator it = netSet.find(cur_net_point);
			if (it != netSet.end()) {

				if ((it->second).second > value.second) {

					(it->second).first = value.first;
					(it->second).second = value.second;
				}
			}
			else {
				netSet[cur_net_point] = value;
			}
		}
	}

	set<int> index;
	vector<int> resultIndexW, resultIndexU;
	for (auto elem : netSet) {
		index.insert((elem.second).first);
	}
	for (int elem : index) {
		resultSet.push_back(W[elem]);
		resultIndexW.push_back(indexW[elem]);
		resultIndexU.push_back(indexU[elem]);
	}
	indexW = resultIndexW;
	indexU = resultIndexU;
	return resultSet;
}

bool notInPhaseConstraint(const vector<vector<double>>& pCon, const double& px, const double& py) {
	bool flag = false;
	int n = pCon.size();
	double D = 0;
	for (int i = 0; i < n - 1; i++)
	{
		D = (px - pCon[i][0]) * (pCon[i + 1][1] - pCon[i][1]) - (py - pCon[i][1]) * (pCon[i + 1][0] - pCon[i][0]);
		if (D >= 0)
			flag = true;
	}
	D = (px - pCon[n - 1][0]) * (pCon[0][1] - pCon[n - 1][1]) - (py - pCon[n - 1][1]) * (pCon[0][0] - pCon[n - 1][0]);
	if (D >= 0)
		flag = true;
	return flag;
}

double sqr_dist(const vector<double>& a, const vector<double>& b) {
	return pow(a[0] - b[0], 2) + pow(a[2] - b[2], 2);
}