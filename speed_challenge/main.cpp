#include "speed_challenge.h"
#include <fstream>
#include <iostream>
#include <ctime>
#include <set>

vector<int> indexW;
vector<int> indexU;

int main()
{
	PointSet W = readFile("W0.txt");
	vector<double> start = W[0];
	PointSet U = readFile("U.txt");
	vector<double> target = readFile("target.txt")[0];
	cout.precision(5);
	double h = 0.1;
	int border = 64000;
	int n = W[0].size();
	int n_n = pow(2, n);
	vector<vector<int>> dx = hyperCube(n_n, n);
	double tau = 0.05;

	int t = clock();
	int index_closest = 0;
	int k;
	double min = 100000;
	for (k = 1; k <= 100; ++k) {
		set<int> index;
		vector<int> resultIndexW, resultIndexU;
		PointSet resultSet;
		W = computeNewSet(W, U, tau, indexW, indexU);
		for (int o = 0; o < W.size(); ++o) {
			//if (notInPhaseConstraint({ {3.0, 0}, {3.5, 0}, {3.5, 0.5}, {3.0, 0.5} }, W[o][0], W[o][1]) &&
				//notInPhaseConstraint({ {3.9, 0.9}, {4.5, 0.9}, {4.5, 1.7}, {3.9, 1.7} }, W[o][0], W[o][1]) &&
				//notInPhaseConstraint({ {5.1, 0.5}, {5.5, 0.5}, {5.5, 1.5}, {5.1, 1.5} }, W[o][0], W[o][1]))
			if (notInPhaseConstraint({ {3.3, -0.3}, {3.5, -0.3}, {3.5, 0.7}, {3.3, 0.7} }, W[o][0], W[o][2]) &&
				notInPhaseConstraint({ {6, -1.8}, {6.2, -1.8}, {6.2, -0.2}, {6, -0.2} }, W[o][0], W[o][2]) &&
				notInPhaseConstraint({ {8.5, 0.5}, {8.5, 0.5}, {8.7, 1.3}, {8.7, 1.3} }, W[o][0], W[o][2]))
				index.insert(o);
		}
		for (int elem : index) {
			resultSet.push_back(W[elem]);
			resultIndexW.push_back(indexW[elem]);
			resultIndexU.push_back(indexU[elem]);
		}
		indexW = resultIndexW;
		indexU = resultIndexU;
		W = resultSet;
		cout << k << "\t" << W.size() << "\t";

		double sqr_d = 0;
		bool isTargetHit = false;
		for (int i = 0; i < W.size(); ++i) {
			sqr_d = sqr_dist(W[i], target);
			if (sqr_d < min)
				min = sqr_d;
			if (sqr_d < 0.007) {
				index_closest = i;
				isTargetHit = true;
				break;
			}

		}
		cout << min << "    ";
		if (isTargetHit) break;

		if (W.size() > border) {
			W = reduceSet(W, h, dx, indexW, indexU);
		}
		cout << W.size() << endl;
		writeFile(W, "W/W" + to_string(k) + ".txt", indexW, indexU);
	}
	cout << "\nTime=" << clock() - t << endl;
	writeFile(W, "result1.txt", indexW, indexU);
	cout << "Finish" << endl;
	cout << index_closest << ":";
	for (int c = 0; c < n; c++) {
		cout << " " << W[index_closest][c];
	}
	cout << " " << indexW[index_closest] << " " << indexU[index_closest] << endl;
	PointSet optimalControl;
	PointSet optimalTrajectory;
	optimalTrajectory.push_back(W[index_closest]);
	optimalControl.push_back(U[indexU[index_closest]]);
	int indexPrevPoint = indexW[index_closest];
	int indexControl;

	for (int m = k - 1; m > 0; --m) {

		vector<double> cur;
		cur.resize(n);
		ifstream ins(("W/W" + to_string(m) + ".txt").c_str());
		for (int i = 0; i < indexPrevPoint; i++) {
			ins.ignore(100, '\n');
		}
		for(int z = 0; z < n; z++)
			ins >> cur[z];
		ins >> indexPrevPoint;
		ins >> indexControl;
		cout << m << ":";
		for (int c = 0; c < n; c++)
			cout << " " << cur[c];
		cout << " " << indexPrevPoint << " " << indexControl << endl;
		optimalTrajectory.push_back(cur);
		optimalControl.push_back(U[indexControl]);
		ins.close();
	}
	cout << optimalTrajectory.size() << " " << optimalControl.size() << endl;

	optimalTrajectory.push_back(start);
	reverse(optimalTrajectory.begin(), optimalTrajectory.end());
	optimalTrajectory.push_back(target);

	ofstream traectoryFile("trajectory.txt");
	for (int i = 0; i < optimalTrajectory.size(); ++i) {
		for (int c = 0; c < n - 1; c++)
			traectoryFile << optimalTrajectory[i][c] << " ";
		traectoryFile << optimalTrajectory[i][n - 1] << "\n";
	}
	traectoryFile.close();

	reverse(optimalControl.begin(), optimalControl.end());

	ofstream controlFile("control1.txt");
	for (int i = 0; i < optimalControl.size(); ++i) {
		controlFile << i + 1 << " " << optimalControl[i][0] << "\n";
	}
	controlFile.close();

	return 0;
}