#define _USE_MATH_DEFINES

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>

using namespace std;

// значение функции f(x) = ctg(x) + x^2
double f_x(double x){
	return (1 / tan(x)) + pow(x, 2);
}

vector<double> gauss1(vector<vector<double>> A, vector<double> b) {
	for (int i = 0; i < A.size(); i++) {
		if (A[i][i] != 0) {
			double ved = A[i][i];
			for (int j = i; j < A.size(); j++) {
				if (A[i][j] != 0) {
					A[i][j] /= ved;
				}
			}
			if (b[i] != 0) {
				b[i] /= ved;
			}

			vector<double> str(A.size());
			double free_b = b[i];
			for (int j = 0; j < A.size(); j++) {
				str[j] = A[i][j];
			}
			for (int j = i + 1; j < A.size(); j++) {
				if (A[j][i] != 0) {
					vector<double> d_str = str;
					free_b = b[i];
					for (int k = 0; k < A.size(); k++) {
						d_str[k] *= A[j][i];
					}
					free_b *= A[j][i];
					for (int k = 0; k < A.size(); k++) {
						A[j][k] -= d_str[k];
					}
					b[j] -= free_b;
				}
			}
		}
		else {
			for (int j = i + 1; j < A.size(); j++) {
				if (A[j][i] != 0) {
					vector<double>Ai;
					vector<double>Aj;
					double bi = b[i];
					double bj = b[j];
					for (int k = 0; k < A.size(); k++) {
						Ai.push_back(A[i][k]);
						Aj.push_back(A[j][k]);
					}
					for (int k = 0; k < A.size(); k++) {
						A[i][k] = Aj[k];
						A[j][k] = Ai[k];
					}
					b[i] = bj;
					b[j] = bi;


					break;
				}
			}
			double ved = A[i][i];
			for (int j = i; j < A.size(); j++) {
				if (A[i][j] != 0) {
					A[i][j] /= ved;
				}
			}
			if (b[i] != 0) {
				b[i] /= ved;
			}

			vector<double> str(A.size());
			double free_b = b[i];
			for (int j = 0; j < A.size(); j++) {
				str[j] = A[i][j];
			}
			for (int j = i + 1; j < A.size(); j++) {
				if (A[j][i] != 0) {
					vector<double> d_str = str;
					free_b = b[i];
					for (int k = 0; k < A.size(); k++) {
						d_str[k] *= A[j][i];
					}
					free_b *= A[j][i];
					for (int k = 0; k < A.size(); k++) {
						A[j][k] -= d_str[k];
					}
					b[j] -= free_b;
				}
			}

		}
	}

	vector<double> x(A.size());
	x[x.size() - 1] = b[b.size() - 1];
	for (int i = x.size() - 2; i >= 0; i--) {
		double res = 0;
		for (int j = i + 1; j < x.size(); j++) {
			res += A[i][j] * x[j];
		}
		x[i] = b[i] - res;
	}
	return x;
}


// интер.полином Лагранжа с равностоящими узлами
vector<pair<double, vector<vector<double>>>> LagrRavn(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
		funcInNodes[i] = f_x(nodes[i]);
	}

	vector<pair<double, vector<vector<double>>>> res;
	for (int i = 0; i < n + 1; i++) {
		double znam = 1;
		vector<vector<double>> work;
		for (int j = 0; j < n + 1; j++) {
			if (i != j) {
				znam *= (nodes[i] - nodes[j]);
				work.push_back({ 1, -nodes[j] });
			}
		}
		double mn = funcInNodes[i] / znam;
		pair<double, vector<vector<double>>> q = make_pair(mn, work);
		res.push_back(q);
	}
	return res;
}
// интерполяционный многочлен Лагранжа с оптимальными узлами
vector<pair<double, vector<vector<double>>>> LagrOpt(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < n + 1; i++) {
		funcInNodes[i] = f_x(nodes[i]);
	}

	vector<pair<double, vector<vector<double>>>> res;
	for (int i = 0; i < n + 1; i++) {
		double znam = 1;
		vector<vector<double>> work;
		for (int j = 0; j < n + 1; j++) {
			if (i != j) {
				znam *= (nodes[i] - nodes[j]);
				work.push_back({ 1, -nodes[j] });
			}
		}
		double mn = funcInNodes[i] / znam;
		pair<double, vector<vector<double>>> q = make_pair(mn, work);
		res.push_back(q);
	}
	return res;

}
// максимальное отклонение Лагранжа
double OtklonLn(vector<pair<double, vector<vector<double>>>> nums, int m, double a, double b) {
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = a + i * ((b - a) / m);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}

	vector<double> LnInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		double res = 0;
		for (int j = 0; j < nums.size(); j++) {
			double pr_res = 1.0;
			for (int k = 0; k < nums[j].second.size(); k++) {
				pr_res *= checkNodes[i] + nums[j].second[k][1];
			}
			pr_res *= nums[j].first;
			res += pr_res;
		}
		LnInCheckNodes[i] = res;
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - LnInCheckNodes[i]) > maxRes) {
			maxRes = abs(funcInCheckNodes[i] - LnInCheckNodes[i]);
		}
	}
	return maxRes;
}
double OtklonLopt(vector<pair<double, vector<vector<double>>>> nums, int m, double a, double b) {
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (m + 1))) + 3.2);
	}
	for (int i = 0; i < m + 1; i++) {
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	
	vector<double> LoptInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		double res = 0;
		for (int j = 0; j < nums.size(); j++) {
			double pr_res = 1.0;
			for (int k = 0; k < nums[j].second.size(); k++) {
				pr_res *= checkNodes[i] + nums[j].second[k][1];
			}
			pr_res *= nums[j].first;
			res += pr_res;
		}
		LoptInCheckNodes[i] = res;
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - LoptInCheckNodes[i]) > maxRes) {
			maxRes = abs(funcInCheckNodes[i] - LoptInCheckNodes[i]);
		}
	}
	return maxRes;
}




vector<pair<double, vector<vector<double>>>> NewtonRavn(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
		funcInNodes[i] = f_x(nodes[i]);
	}
	// формируем матрицу для разделенных разностей
	
	vector<vector<double>> razn;
	for (int i = 0; i < n+1; i++) {
		vector<double> q(n+1);
		razn.push_back(q);
	}

	// заполняем матрицу разделённых разностей
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j < n+1 - i; j++) {
			if (i == 0) {
				razn[j][i] = funcInNodes[j];
			}
			else {
				razn[j][i] = (razn[j + 1][i - 1] - razn[j][i - 1]) / (nodes[i + j] - nodes[j]);
			}
		}
	}
	vector<pair<double, vector<vector<double>>>> res;
	for (int i = 0; i < n + 1; i++) {
		if (i == 0) {
			double sv = razn[i][i];
			vector<vector<double>> q = { {} };
			res.push_back(make_pair(sv, q));
		}
		else {
			double sv = razn[0][i];
			vector<vector<double>> q = {};
			for (int j = 0; j < i; j++) {
				q.push_back({ 1, -nodes[j] });
			}
			res.push_back(make_pair(sv, q));
		}
	}
	return res;
}
vector<pair<double, vector<vector<double>>>> NewtonOpt(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
		funcInNodes[i] = f_x(nodes[i]);
	}
	// формируем матрицу для разделенных разностей
	vector<vector<double>> razn;
	for (int i = 0; i < n + 1; i++) {
		vector<double> q(n + 1);
		razn.push_back(q);
	}
	// заполняем матрицу разделённых разностей
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1 - i; j++) {
			if (i == 0) {
				razn[j][i] = funcInNodes[j];
			}
			else {
				razn[j][i] = (razn[j + 1][i - 1] - razn[j][i - 1]) / (nodes[i + j] - nodes[j]);
			}
		}
	}

	vector<pair<double, vector<vector<double>>>> res;
	for (int i = 0; i < n + 1; i++) {
		if (i == 0) {
			double sv = razn[i][i];
			vector<vector<double>> q = { {} };
			res.push_back(make_pair(sv, q));
		}
		else {
			double sv = razn[0][i];
			vector<vector<double>> q = {};
			for (int j = 0; j < i; j++) {
				q.push_back({ 1, -nodes[j] });
			}
			res.push_back(make_pair(sv, q));
		}
	}
	return res;
}
double OtklonNn(vector<pair<double, vector<vector<double>>>> nums, int m, double a, double b) {
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = a + i * ((b - a) / m);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}

	vector<double> NnInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		double res = 0;
		for (int j = 0; j < nums.size(); j++) {
			if (j == 0) {
				res += nums[j].first;
			}
			else {
				double pr_res = nums[j].first;
				for (int k = 0; k < nums[j].second.size(); k++) {
					pr_res *= (checkNodes[i] + nums[j].second[k][1]);
				}
				res += pr_res;
			}
		}
		NnInCheckNodes[i] = res;
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - NnInCheckNodes[i]) > maxRes) {
			maxRes = abs(funcInCheckNodes[i] - NnInCheckNodes[i]);
		}
	}
	return maxRes;
}
double OtklonNopt(vector<pair<double, vector<vector<double>>>> nums, int m, double a, double b) {
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (m + 1))) + 3.2);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	vector<double> NoptInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		double res = 0;
		for (int j = 0; j < nums.size(); j++) {
			if (j == 0) {
				res += nums[j].first;
			}
			else {
				double pr_res = nums[j].first;
				for (int k = 0; k < nums[j].second.size(); k++) {
					pr_res *= (checkNodes[i] + nums[j].second[k][1]);
				}
				res += pr_res;
			}
		}
		NoptInCheckNodes[i] = res;
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - NoptInCheckNodes[i]) > maxRes) {
			maxRes = abs(funcInCheckNodes[i] - NoptInCheckNodes[i]);
		}
	}
	return maxRes;
}




// линейный сплайн S(1,0) с равностоящими узлами
vector<vector<double>> S10(double a, double b, int n) {
	// узлы и значение функции в узлах
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
		funcInNodes[i] = f_x(nodes[i]);
	}

	// lk
	vector<vector<double>> lk;
	for (int i = 0; i < n; i++) {
		vector<double> a(2);
		lk.push_back(a);
	}

	for (int i = 0; i < lk.size(); i++) {
		lk[i][0] = (double)(funcInNodes[i + 1] - funcInNodes[i]) / (nodes[i + 1] - nodes[i]);
		lk[i][1] = funcInNodes[i + 1] - nodes[i + 1] * lk[i][0];
	}

	return lk;
}
// линейный сплайн S(1,0) с оптимальными узлами
vector<vector<double>> S10_opt(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < n + 1; i++) {
		funcInNodes[i] = f_x(nodes[i]);
	}
	// lk
	vector<vector<double>> lk;
	for (int i = 0; i < n; i++) {
		vector<double> a(2);
		lk.push_back(a);
	}

	for (int i = 0; i < lk.size(); i++) {
		lk[i][0] = (double)(funcInNodes[i + 1] - funcInNodes[i]) / (nodes[i + 1] - nodes[i]);
		lk[i][1] = funcInNodes[i + 1] - nodes[i + 1] * lk[i][0];
	}

	return lk;
}
double OtklonS10(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = a + i * ((b - a) / m);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
	}
	vector<double> S10InCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		int y=-1;
		for (int j = 0; j < n; j++) {
			if (checkNodes[i] >= nodes[j] && checkNodes[i] <= nodes[j + 1]) {
				y = j;
				break;
			}
		}
		if (y == -1) {
			y = n - 1;
		}
		S10InCheckNodes[i] = nums[y][0] * checkNodes[i] + nums[y][1];
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(S10InCheckNodes[i] - funcInCheckNodes[i]) > maxRes) {
			maxRes = abs(S10InCheckNodes[i] - funcInCheckNodes[i]);
		}
	}
	return maxRes;
}
double OtklonS10opt(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (m + 1))) + 3.2);
	}
	sort(checkNodes.begin(), checkNodes.end());
	for (int i = 0; i < m + 1; i++) {
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}

	vector<double> S10optInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		int y;
		if (checkNodes[i] <= nodes[0]) {
			y = 0;
		}
		else {
			if (checkNodes[i] >= nodes[nodes.size() - 1]) {
				y = nums.size() - 1;
			}
			else {
				for (int j = 0; j < n; j++) {
					if (checkNodes[i] >= nodes[j] && checkNodes[i] <= nodes[j + 1]) {
						y = j;
						break;
					}
				}
			}
		}
		S10optInCheckNodes[i] = nums[y][0] * checkNodes[i] + nums[y][1];
	}
	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - S10optInCheckNodes[i]) > maxRes)
			maxRes = abs(funcInCheckNodes[i] - S10optInCheckNodes[i]);
	}
	return maxRes;
}




vector<vector<double>> S21(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
		funcInNodes[i] = f_x(nodes[i]);
	}
	vector<vector<double>> lk;
	for (int i = 0; i < n; i++) {
		vector<double> a(3);
		lk.push_back(a);
	}

	// матрица коэф-тов
	vector<vector<double>> A;
	for (int i = 0; i < n * 3; i++) {
		vector<double> a(3 * lk.size());
		A.push_back(a);
	}
	// вектор свободных членов
	vector<double> freeB(3 * n);

	int node = 0;
	int col = 0;

	for (int i = 0; i < (nodes.size() - 2) * 2 + 2; i++) {
		A[i][col] = pow(nodes[node], 2);
		A[i][col + 1] = nodes[node];
		A[i][col + 2] = 1;
		freeB[i] = funcInNodes[node];
		if (i % 2 == 0) {
			node++;
		}
		else {

			col += 3;
		}
	}
	int row = (nodes.size() - 2) * 2 + 2;
	col = 0;
	for (int i = 0; i < nodes.size() - 2; i++) {
		A[row][col] = 2 * nodes[i + 1];
		A[row][col + 1] = 1;
		A[row][col + 3] = -2 * nodes[i + 1];
		A[row][col + 4] = -1;
		row++;
		col += 3;
	}
	A[A.size() - 1][A.size() - 3] = 2 * nodes[nodes.size() - 1];
	A[A.size() - 1][A.size() - 2] = 1;
	
	// находим коэф-ты пр помощи решения СЛАУ методом Гаусса
	vector<double> koef = gauss1(A, freeB);

	int ind = 0;
	for (int i = 0; i < lk.size(); i++) {
		lk[i][0] = koef[ind];
		lk[i][1] = koef[ind + 1];
		lk[i][2] = koef[ind + 2];
		ind += 3;
	}
	return lk;
}
vector<vector<double>> S21_opt(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < n + 1; i++) {
		funcInNodes[i] = f_x(nodes[i]);
	}
	vector<vector<double>> lk;
	for (int i = 0; i < n; i++) {
		vector<double> a(3);
		lk.push_back(a);
	}
	// матрица коэф-тов
	vector<vector<double>> A;
	for (int i = 0; i < n * 3; i++) {
		vector<double> a(3 * lk.size());
		A.push_back(a);
	}
	// вектор свободных членов
	vector<double> freeB(3 * n);

	int node = 0;
	int col = 0;

	for (int i = 0; i < (nodes.size() - 2) * 2 + 2; i++) {
		A[i][col] = pow(nodes[node], 2);
		A[i][col + 1] = nodes[node];
		A[i][col + 2] = 1;
		freeB[i] = funcInNodes[node];
		if (i % 2 == 0) {
			node++;
		}
		else {

			col += 3;
		}
	}
	int row = (nodes.size() - 2) * 2 + 2;
	col = 0;
	for (int i = 0; i < nodes.size() - 2; i++) {
		A[row][col] = 2 * nodes[i + 1];
		A[row][col + 1] = 1;
		A[row][col + 3] = -2 * nodes[i + 1];
		A[row][col + 4] = -1;
		row++;
		col += 3;
	}
	A[A.size() - 1][A.size() - 3] = 2 * nodes[nodes.size() - 1];
	A[A.size() - 1][A.size() - 2] = 1;

	// находим коэф-ты пр помощи решения СЛАУ методом Гаусса
	vector<double> koef = gauss1(A, freeB);

	int ind = 0;
	for (int i = 0; i < lk.size(); i++) {
		lk[i][0] = koef[ind];
		lk[i][1] = koef[ind + 1];
		lk[i][2] = koef[ind + 2];
		ind += 3;
	}
	return lk;
}
double OtklonS21(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = a + i * ((b - a) / m);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
	}
	vector<double> S21InCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		int y = -1;
		for (int j = 0; j < n; j++) {
			if (checkNodes[i] >= nodes[j] && checkNodes[i] <= nodes[j + 1]) {
				y = j;
				break;
			}
		}
		if (y == -1) {
			y = n - 1;
		}
		S21InCheckNodes[i] = nums[y][0] * pow(checkNodes[i], 2) + nums[y][1] * checkNodes[i] + nums[y][2];
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(S21InCheckNodes[i] - funcInCheckNodes[i]) > maxRes) {
			maxRes = abs(S21InCheckNodes[i] - funcInCheckNodes[i]);
		}
	}
	return maxRes;
}
double OtklonS21opt(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (m + 1))) + 3.2);
	}
	sort(checkNodes.begin(), checkNodes.end());
	for (int i = 0; i < m + 1; i++) {
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	vector<double> S21optInCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		int y;
		if (checkNodes[i] <= nodes[0]) {
			y = 0;
		}
		else {
			if (checkNodes[i] >= nodes[nodes.size() - 1]) {
				y = nums.size() - 1;
			}
			else {
				for (int j = 0; j < n; j++) {
					if (checkNodes[i] >= nodes[j] && checkNodes[i] <= nodes[j + 1]) {
						y = j;
						break;
					}
				}
			}
		}
		S21optInCheckNodes[i] = nums[y][0] * pow(checkNodes[i], 2) + nums[y][1] * checkNodes[i] + nums[y][2];
	}
	for (int i = 0; i < m + 1; i++) {
		if (abs(funcInCheckNodes[i] - S21optInCheckNodes[i]) > maxRes)
			maxRes = abs(funcInCheckNodes[i] - S21optInCheckNodes[i]);
	}
	return maxRes;
}




vector<vector<double>> S32(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
		funcInNodes[i] = f_x(nodes[i]);
	}
	double h = nodes[1] - nodes[0];
	vector<vector<double>> H;
	for (int i = 0; i < n - 1; i++) {
		vector<double> a(n - 1);
		H.push_back(a);
	}
	for (int i = 0; i < n-1; i++) {
		if (i == 0) {
			H[i][i] = 4 * h;
			H[i][i + 1] = h;
		}
		else {
			if (i == n - 2) {
				H[i][i] = 4 * h;
				H[i][i - 1] = h;
			}
			else {
				H[i][i - 1] = h;
				H[i][i] = 4 * h;
				H[i][i + 1] = h;
			}
		}
	}
	vector<double> gamma(n-1);
	for (int i = 1; i < n; i++) {
		gamma[i - 1] = 6 * ((funcInNodes[i + 1] - 2 * funcInNodes[i] + funcInNodes[i - 1]) / h);
	}
	vector<double> vt_pr = gauss1(H, gamma);
	vector<double> second_der(n+1);
	for (int i = 0; i < n+1; i++) {
		if (i == 0 || i == n) {
			second_der[i] = 0;
		}
		else {
			second_der[i] = vt_pr[i - 1];
		}
	}

	vector<double> first_der(n);
	for (int i = 0; i < n; i++) {
		first_der[i] = (funcInNodes[i + 1] - funcInNodes[i]) / h - second_der[i + 1] * h / 6 - second_der[i] * h / 3;
	}
	vector<vector<double>> koef;
	for (int i = 0; i < n; i++) {
		vector<double> a(4);
		a[0] = (second_der[i + 1] - second_der[i]) / (6 * h);
		a[1] = second_der[i] / 2;
		a[2] = first_der[i];
		a[3] = funcInNodes[i];
		koef.push_back(a);
	}
	return koef;
}
vector<vector<double>>S32_opt(double a, double b, int n) {
	vector<double> nodes(n + 1);
	vector<double> funcInNodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < n + 1; i++) {
		funcInNodes[i] = f_x(nodes[i]);
	}
	vector<double> h(n);
	for (int i = 0; i < n; i++) {
		h[i] = nodes[i + 1] - nodes[i];
	}

	vector<vector<double>> H;
	for (int i = 0; i < n-1; i++) {
		vector<double> a(n-1);
		H.push_back(a);
	}
	for (int i = 0; i < n-1; i++) {
		if (i == 0) {
			H[i][i] = 2 * (h[i] + h[i + 1]);
			H[i][i + 1] = h[i + 1];
		}
		else {
			if (i == n-2) {
				H[i][i] = 2 * (h[i] + h[i + 1]);
				H[i][i - 1] = h[i];
			}
			else {
				H[i][i - 1] = h[i];
				H[i][i] = 2 * (h[i] + h[i + 1]);
				H[i][i + 1] = h[i + 1];
			}
		}
	}
	vector<double> gamma(n-1);
	for (int i = 1; i < n; i++) {
		gamma[i - 1] = 6 * ((funcInNodes[i + 1] - funcInNodes[i]) / h[i] - (funcInNodes[i] - funcInNodes[i - 1]) / h[i - 1]);
	}
	vector<double> vt_pr = gauss1(H, gamma);
	vector<double> second_der(n+1);
	for (int i = 0; i < n+1; i++) {
		if (i == 0 || i == n) {
			second_der[i] = 0;
		}
		else {
			second_der[i] = vt_pr[i - 1];
		}
	}
	vector<double> first_der(n);
	for (int i = 0; i < n; i++) {
		first_der[i] = (funcInNodes[i + 1] - funcInNodes[i]) / h[i] - (second_der[i + 1] * h[i]) / 6 - (second_der[i] * h[i]) / 3;
	}
	vector<vector<double>> koef;
	for (int i = 0; i < n; i++) {
		vector<double> a(4);
		a[0] = (second_der[i + 1] - second_der[i]) / (6 * h[i]);
		a[1] = second_der[i] / 2;
		a[2] = first_der[i];
		a[3] = funcInNodes[i];
		koef.push_back(a);
	}
	return koef;
}
double OtklonS32(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = a + i * ((b - a) / m);
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = a + i * ((b - a) / n);
	}
	vector<double> S32InCheckNodes(m + 1);
	for (int i = 0; i < m + 1; i++) {
		int y = -1;
		for (int j = 0; j < n; j++) {
			if (checkNodes[i] >= nodes[j] && checkNodes[i] <= nodes[j + 1]) {
				y = j;
				break;
			}
		}
		if (y == -1) {
			y = n - 1;
		}
		for (int r = 0; r < nums[y].size(); r++) {
			S32InCheckNodes[i] += nums[y][r] * pow(checkNodes[i] - nodes[y], 4 - 1 - r);
		}
	}

	for (int i = 0; i < m + 1; i++) {
		if (abs(S32InCheckNodes[i] - funcInCheckNodes[i]) > maxRes) {
			maxRes = abs(S32InCheckNodes[i] - funcInCheckNodes[i]);
		}
	}
	return maxRes;
}
double OtklonS32opt(vector<vector<double>> nums, int m, double a, double b) {
	int n = nums.size();
	double maxRes = -1;
	vector<double> checkNodes(m + 1);
	vector<double> funcInCheckNodes(m + 1);
	vector<double> nodes(n + 1);
	for (int i = 0; i < n + 1; i++) {
		nodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (n + 1))) + 3.2);
	}
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < m + 1; i++) {
		checkNodes[i] = 0.5 * (2.8 * cos((M_PI * (2 * i + 1)) / (2 * (m + 1))) + 3.2);
	}
	sort(checkNodes.begin(), checkNodes.end());
	for (int i = 0; i < m + 1; i++) {
		funcInCheckNodes[i] = f_x(checkNodes[i]);
	}

	for (int i = 0; i < checkNodes.size(); i++) {
		int j;
		if (checkNodes[i] < nodes[0]) {
			j = 0;
		}
		else {
			if (checkNodes[i] > nodes[n - 1]) {
				j = nums.size() - 1;
			}
			else {
				if (checkNodes[i] >= nodes[0] && checkNodes[i] <= nodes[n - 1]) {
					for (int y = 0; y < nodes.size() - 1; y++) {
						if (checkNodes[i] >= nodes[y] && checkNodes[i] <= nodes[y + 1]) {
							j = y;
						}
					}
				}
			}
		}
		double f = f_x(checkNodes[i]);
		double p = 0;
		for (int r = 0; r < nums[j].size(); r++) {
			p += nums[j][r] * pow(checkNodes[i] - nodes[j], 4 - 1 - r);
		}
		if (abs(f - p) > maxRes) {
			maxRes = abs(f - p);
		}
	}
	return maxRes;

}
int main() {
	double a = 0.2;
	double b = 3;

	string method;
	cout << "Enter method " << '\n';
	cin >> method;

	if (method == "LagrRavn") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<pair<double, vector<vector<double>>>> res = LagrRavn(a, b, n);
		for (int i = 0; i < n + 1; i++) {
			cout << res[i].first << "   ";
			for (int j = 0; j < res[i].second.size(); j++) {
				cout << res[i].second[j][0] << " " << res[i].second[j][1] << "   ";
			}
			cout << '\n';
		}
	}
	if (method == "LagrOpt") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<pair<double, vector<vector<double>>>> res = LagrOpt(a, b, n);
		for (int i = 0; i < n + 1; i++) {
			cout << res[i].first << "   ";
			for (int j = 0; j < res[i].second.size(); j++) {
				cout << res[i].second[j][0] << " " << res[i].second[j][1] << "   ";
			}
			cout << '\n';
		}
	}
	if (method == "OtklonL") {
		vector<int>N(8);
		vector<int>M(8);
		vector<double> Ln(8);
		vector<double> Lopt(8);
		for (int i = 0; i < 8; i++) {
			int n, m;
			cout << "Enter the amount of nodes (minimum 2) and the amount of CheckNodes" << '\n';
			cin >> n >> m;
			N[i] = n;
			M[i] = m;
			Ln[i] = OtklonLn(LagrRavn(a, b, n), m, a, b);
			Lopt[i] = OtklonLopt(LagrOpt(a, b, n), m, a, b);
		}
		for (int i = 0; i < 8; i++) {
			cout << N[i] << " " << M[i] << " " << Ln[i] << " " << Lopt[i] << '\n';
		}
	}

	if (method == "NewRavn") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<pair<double, vector<vector<double>>>> res = NewtonRavn(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			cout << res[i].first << "   ";
			for (int j = 0; j < res[i].second.size(); j++) {
				if (res[i].second[j].size() == 0) {
					continue;
				}
				else {
					cout << res[i].second[j][0] << " " << res[i].second[j][1] << "   ";
				}
			}
			cout << '\n';
		}
	}
	if (method == "NewOpt") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<pair<double, vector<vector<double>>>> res = NewtonOpt(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			cout << res[i].first << "   ";
			for (int j = 0; j < res[i].second.size(); j++) {
				if (res[i].second[j].size() == 0) {
					continue;
				}
				else {
					cout << res[i].second[j][0] << " " << res[i].second[j][1] << "   ";
				}
			}
			cout << '\n';
		}
	}
	if (method == "OtklonN") {
		vector<int>N(7);
		vector<int>M(7);
		vector<double> Nn(7);
		vector<double> Nopt(7);
		for (int i = 0; i < 7; i++) {
			int n, m;
			cout << "Enter the amount of nodes (minimum 2) and the amount of CheckNodes" << '\n';
			cin >> n >> m;
			N[i] = n;
			M[i] = m;
			Nn[i] = OtklonNn(NewtonRavn(a, b, n), m, a, b);
			Nopt[i] = OtklonNopt(NewtonOpt(a, b, n), m, a, b);
		}
		for (int i = 0; i < 7; i++) {
			cout << N[i] << " " << M[i] << " " << Nn[i] << " " << Nopt[i] << '\n';
		}
	}

	if (method == "S10") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S10(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "S10_opt") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S10_opt(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "OtklonS10") {
		vector<int> amountN(4);
		vector<int> amountM(4);
		vector<double> resMaxS10(4);
		vector<double> resMaxS10opt(4);
		for (int i = 0; i < 4; i++) {
			int n, m;
			cout << "Enter the amount of nodes (minimum 2) and the amount of CheckNodes" << '\n';
			cin >> n >> m;
			amountN[i] = n;
			amountM[i] = m;
			resMaxS10[i] = OtklonS10(S10(a, b, n),m, a, b);
			resMaxS10opt[i] = OtklonS10opt(S10_opt(a, b, n), m, a, b);
		}
		for (int i = 0; i < 4; i++) {
			cout << amountN[i] << " " << amountM[i] << " " << resMaxS10[i] << " " << resMaxS10opt[i] << '\n';
		}
	}

	if (method == "S21") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S21(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "S21_opt") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S21_opt(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "OtklonS21") {
		vector<int> amountN(4);
		vector<int> amountM(4);
		vector<double> resMaxS21(4);
		vector<double> resMaxS21opt(4);
		for (int i = 0; i < 4; i++) {
			int n, m;
			cout << "Enter the amount of nodes (minimum 2) and the amount of CheckNodes" << '\n';
			cin >> n >> m;
			amountN[i] = n;
			amountM[i] = m;
			resMaxS21[i] = OtklonS21(S21(a, b, n), m, a, b);
			resMaxS21opt[i] = OtklonS21opt(S21_opt(a, b, n), m, a, b);
		}
		for (int i = 0; i < 4; i++) {
			cout << amountN[i] << " " << amountM[i] << " " << resMaxS21[i] << " " << resMaxS21opt[i] << '\n';
		}
	}

	if (method == "S32") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S32(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "S32_opt") {
		int n;
		cout << "Enter the amount of nodes (minimum 2)" << '\n';
		cin >> n;
		vector<vector<double>> res = S32_opt(a, b, n);
		for (int i = 0; i < res.size(); i++) {
			for (int j = 0; j < res[i].size(); j++) {
				cout << res[i][j] << " ";
			}
			cout << '\n';
		}
	}
	if (method == "OtklonS32") {
		vector<int> amountN(4);
		vector<int> amountM(4);
		vector<double> resMaxS32(4);
		vector<double> resMaxS32opt(4);
		for (int i = 0; i < 1; i++) {
			int n, m;
			cout << "Enter the amount of nodes (minimum 2) and the amount of CheckNodes" << '\n';
			cin >> n >> m;
			amountN[i] = n;
			amountM[i] = m;
			resMaxS32[i] = OtklonS32(S32(a, b, n), m, a, b);
			resMaxS32opt[i] = OtklonS32opt(S32_opt(a, b, n), m, a, b);
		}
		for (int i = 0; i < 4; i++) {
			cout << amountN[i] << " " << amountM[i] << " " << resMaxS32[i] << " " << resMaxS32opt[i] << '\n';
		}
	}
}