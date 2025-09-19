#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <random>
#include <set>
#include <iomanip>

using namespace std;


double f_x(double x) {
	return x * log(pow(x + 2, 2));
}

double polynom(vector<double> a, double x) {
	double res = 0;
	for (int i = 0; i < a.size(); i++) {
		res += a[i] * pow(x, i);
	}
	return res;
}

// метод Гаусса для решения СЛАУ
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

double generate_real_random(double mini, double maxi) {
	mt19937 rng(random_device{}());
	uniform_real_distribution<double> distribution(mini, maxi);
	double random_number = distribution(rng);
	return random_number;
}

pair<vector<double>, vector<vector<double>>> generate_nodes(int m, int m_f, double eps) {
	vector<double> x(m+1);
	x[0] = -1;
	x[m] = 1;
	for (int i = 1; i < m; i++) {
		x[i] = x[i - 1] + double(2) / double(m);
	}

	vector<vector<double>> f;
	for (int i = 0; i < m + 1; i++) {
		double f_xi = f_x(x[i]);
		double down = f_xi - eps;
		double up = f_xi + eps;

		vector<double> temp;
		while (temp.size() != m_f) {
			double temp_zn = generate_real_random(down, up);
			int k = 0;
			for (int i = 0; i < temp.size(); i++) {
				if (temp[i] == temp_zn) {
					k++;
					break;
				}
			}
			if (k == 0) {
				temp.push_back(temp_zn);
			}
		}
		sort(temp.begin(), temp.end());
		f.push_back(temp);
	}
	return make_pair(x, f);
}

vector<vector<double>> transpon(vector<vector<double>> A) {
	int n = A.size();
	int m = A[0].size();
	vector<vector<double>> res(m, vector<double>(n));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			res[i][j] = A[j][i];
		}
	}
	return res;
}

vector<vector<double>> matrix_mult(vector<vector<double>> m1, vector<vector<double>> m2) {
	int n = m1.size();
	vector<vector<double>> res(n, vector<double>(n));
	for (int i = 0; i < m1.size(); i++) {
		for (int j = 0; j < m2[0].size(); j++) {
			for (int k = 0; k < m2.size(); k++) {
				res[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return res;
}

vector<double> matrix_vector_mult(vector<vector<double>> m, vector<double> v) {
	vector<double> res(m.size());
	for (int i = 0; i < m.size(); i++) {
		for (int j = 0; j < v.size(); j++) {
			res[i] += m[i][j] * v[j];
		}
	}
	return res;
}

double summa_mist(vector<double> x, vector<double> a) {
	double res = 0;
	for (int i = 0; i < x.size(); i++) { // а тут m - кол-во иксов или общее кол-во 3*m
		res += pow(f_x(x[i]) - polynom(a, x[i]), 2);
	}
	return res;
}

int main() {
	setlocale(LC_ALL, "Russian");

	double eps = pow(10, -3);

	int n;
	cout << "Степень полинома " << '\n';
	cin >> n;

	int m;
	cout << "Количество узлов" << '\n';
	cin >> m;

	int m_f;
	cout << "Количество значений y в каждом узле" << '\n';
	cin >> m_f;

	// НОРМАЛЬНЫЕ УРАВНЕНИЯ
	pair<vector<double>, vector<vector<double>>> temp = generate_nodes(m, m_f, eps);
	vector<double> x = temp.first;  // генерируем x
	vector<vector<double>> y = temp.second;  // значения y в точках х

	// формируем f, E
	//Et * E * a = Et * f
	vector<double> f_ur; // вектор свободных членов
	vector<vector<double>> E; // матрица Е
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < m_f; j++) {
			// добавляем в f_ur y[i][j]
			f_ur.push_back(y[i][j]);
			// формируем вектор temp_e
			vector<double> temp_e;
			for (int k = 0; k < n + 1; k++) {
				temp_e.push_back(pow(x[i], k));
			}
			E.push_back(temp_e);
		}
	}

	// решаем систему и находим коэф-ты
	vector<vector<double>> E_t = transpon(E);
	vector<vector<double>> E_tE = matrix_mult(E_t, E);
	vector<double> free_b = matrix_vector_mult(E_t, f_ur);
	vector<double> a = gauss1(E_tE, free_b);



	cout << "Нормальные уравнения" << '\n';
	// a0 - свободный член, a1 - при первой степени ...
	for (int i = 0; i < a.size(); i++) {
		cout << a[i] << " ";
	}
	cout << '\n';


	// ОРТОГОНАЛЬНЫЕ УРАВНЕНИЯ
	cout << "Ортогональные уравнения" << '\n';

	// делаем удобный вид для x, y
	vector<double> x_temp;
	vector<double> y_temp;
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < m_f; j++) {
			x_temp.push_back(x[i]);
			y_temp.push_back(y[i][j]);
		}
	}

	// все полиномы q
	vector<vector<double>> q;
	q.push_back({ 1 });

	// q1 считаем
	double s1 = 0;
	for (int i = 0; i < x_temp.size(); i++) {
		s1 += x_temp[i];
	}
	q.push_back({ -s1 / x_temp.size(), 1 });

	// считаем оставшиеся q 
	for (int j = 1; j < n; j++) {
		vector<double> q_curr(j + 2);
		double alpha_v = 0;
		double alpha_n = 0;
		double beta_v = 0;
		double beta_n = 0;
		for (int i = 0; i < x_temp.size(); i++) {
			alpha_v += x_temp[i] * pow(polynom(q[j], x_temp[i]), 2);
			alpha_n += pow(polynom(q[j], x_temp[i]), 2);
			beta_v += x_temp[i] * polynom(q[j], x_temp[i]) * polynom(q[j - 1], x_temp[i]);
			beta_n += pow(polynom(q[j-1], x_temp[i]), 2);
		}
		double alpha = alpha_v/alpha_n;
		double beta = beta_v/beta_n;
		for (int i = 0; i < q[j].size(); i++) {
			q_curr[i + 1] += q[j][i];
		}
		for (int i = 0; i < q[j].size(); i++) {
			q_curr[i] -= alpha * q[j][i];
		}
		for (int i = 0; i < q[j - 1].size(); i++) {
			q_curr[i] -= beta * q[j - 1][i];
		}
		q.push_back(q_curr);
	}

	// найти a0, ..., a_n
	vector<double> a_koef(n + 1);
	for (int k = 0; k <= n; k++) {
		double a_v = 0;
		double a_n = 0;
		for (int i = 0; i < x_temp.size(); i++) {
			// числительи знаменатель
			a_v += y_temp[i] * polynom(q[k], x_temp[i]);
			a_n += pow(polynom(q[k], x_temp[i]), 2);
		}
		a_koef[k] = a_v / a_n;
	}

	// итоговый ортогональный полином - sum a_k * q_k
	vector<double> q_res(n + 1);
	for (int i = 0; i < a_koef.size(); i++) {
		for (int j = 0; j < q[i].size(); j++) {
			q_res[j] += a_koef[i] * q[i][j];
		}
	}
	for (int i = 0; i < q_res.size(); i++) {
		cout << q_res[i] << " ";
	}
	cout << '\n';
	//СЧИТАЕМ СУММУ КВАДРАТОВ ОШИБОК
	cout << "Нормальные уравнения " << summa_mist(x_temp, a) << '\n';
	cout << "Ортогональные уравнения " << summa_mist(x_temp, q_res);
}