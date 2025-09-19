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

// определния минора
vector<vector<double>> minor(vector<vector<double>> matrix, int row, int col) {
	int n = matrix.size();
	vector<vector<double>> res;
	for (int i = 0; i < n - 1; i++) {
		vector<double> temp(n - 1);
		res.push_back(temp);
	}
	if (row == 0) {
		if (col == 0) {
			for (int i = 1; i < n; i++) {
				for (int j = 1; j < n; j++) {
					res[i - 1][j - 1] = matrix[i][j];
				}
			}
		}
		else {
			if (col == n - 1) {
				for (int i = 1; i < n; i++) {
					for (int j = 0; j < n - 1; j++) {
						res[i - 1][j] = matrix[i][j];
					}
				}
			}
			else {
				for (int i = 1; i < n; i++) {
					for (int j = 0; j < col; j++) {
						res[i - 1][j] = matrix[i][j];
					}
				}
				for (int i = 1; i < n; i++) {
					for (int j = col + 1; j < n; j++) {
						res[i - 1][j - 1] = matrix[i][j];
					}
				}
			}
		}
	}
	else {
		if (row == n - 1) {
			if (col == 0) {
				for (int i = 0; i < n-1; i++) {
					for (int j = 1; j < n; j++) {
						res[i][j - 1] = matrix[i][j];
					}
				}
			}
			else {
				if (col == n - 1) {
					for (int i = 0; i < n - 1; i++) {
						for (int j = 0; j < n - 1; j++) {
							res[i][j] = matrix[i][j];
						}
					}
				}
				else {
					for (int i = 0; i < n-1; i++) {
						for (int j = 0; j < col; j++) {
							res[i][j] = matrix[i][j];
						}
					}
					for (int i = 0; i < n-1; i++) {
						for (int j = col + 1; j < n; j++) {
							res[i][j - 1] = matrix[i][j];
						}
					}
				}
			}
		}
		else {
			if (col == 0) {
				for (int i = 0; i < row; i++) {
					for (int j = 1; j < n;j++) {
						res[i][j - 1] = matrix[i][j];
					}
				}
				for (int i = row + 1; i < n; i++) {
					for (int j = 1; j < n; j++) {
						res[i - 1][j - 1] = matrix[i][j];
					}
				}
			}
			else {
				if (col == n - 1) {
					for (int i = 0; i < row; i++) {
						for (int j = 0; j < n-1; j++) {
							res[i][j] = matrix[i][j];
						}
					}
					for (int i = row + 1; i < n; i++) {
						for (int j = 0; j < n - 1; j++) {
							res[i - 1][j] = matrix[i][j];
						}
					}
				}
				else {
					for (int i = 0; i < row; i++) {
						for (int j = 0; j < col; j++) {
							res[i][j] = matrix[i][j];
						}
					}
					for (int i = 0; i < row; i++) {
						for (int j = col + 1; j < n; j++) {
							res[i][j - 1] = matrix[i][j];
						}
					}
					for (int i = row + 1; i < n; i++) {
						for (int j = 0; j < col; j++) {
							res[i - 1][j] = matrix[i][j];
						}
					}
					for (int i = row + 1; i < n; i++) {
						for (int j = col + 1; j < n; j++) {
							res[i - 1][j - 1] = matrix[i][j];
						}
					}
				}
			}
		}
	}
	return res;
}

// функция для нахождения нормы вектора
double norm_vector(vector<double> y) {
	double summa = 0;
	for (int i = 0; i < y.size(); i++) {
		summa += y[i] * y[i];
	}
	return pow(summa, 0.5);
}

// перемножение матриц
vector<vector<double>> matrix_mult(vector<vector<double>> m1, vector<vector<double>> m2) {
	int n = m1.size();
	vector<vector<double>> res;
	for (int i = 0; i < n; i++) {
		vector<double> temp(n);
		res.push_back(temp);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				res[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return res;
}

// QR разложение матрицы
pair<vector<vector<double>>, vector<vector<double>>> QR_razl(vector<vector<double>> A) {
	// короче как-то не так получается матрица Q (посчитать ручками)
	int n = A.size();
	vector<vector<vector<double>>> Q;
	vector<vector<vector<double>>> R;

	// как бы нулевой шаг
	vector<vector<double>> Q0;
	vector<vector<double>> R0 = A;
	for (int i = 0; i < n; i++) {
		vector<double> temp(n);
		temp[i] = 1;
		Q0.push_back(temp);
	}
	Q.push_back(Q0);
	R.push_back(R0);

	// начинаем нормальные шаги
	for (int i = 0; i < n-1; i++) {
		// создаем векторы y, z
		vector<double> y(n - i);
		vector<double> z(n - i);
		z[0] = 1;
		for (int j = i; j < n; j++) {
			y[j-i] = R[R.size() - 1][j][i];
		}

		// считаем альфа 
		double alfa = norm_vector(y);

		// считаем вектор, норму котрого в r
		vector<double> t(n - i);
		for (int j = 0; j < n - i; j++) {
			t[j] = y[j] - alfa * z[j];
		}
		double r = norm_vector(t);

		// вектор w - должен быть единичным
		vector<double> w(n - i);
		for (int j = 0; j < n - i; j++) {
			w[j] = t[j] / r;
		}

		vector<vector<double>> Q_temp;
		vector<vector<double>> R_temp;
		if (i == 0) {
			// на первом шаге - Q, R
			for (int j = 0; j < n; j++) {
				vector<double> tempor(n);
				tempor[j] = 1;
				Q_temp.push_back(tempor);
			}
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					Q_temp[j][k] -= 2 * w[j] * w[k];
				}
			}
			R_temp = matrix_mult(Q_temp, A);
			Q.push_back(Q_temp);
			R.push_back(R_temp);
		}
		else {
			vector<vector<double>> Qt_t;
			for (int j = 0; j < n - i; j++) {
				vector<double> t(n - i);
				t[j] = 1;
				Qt_t.push_back(t);
			}
			for (int j = 0; j < n - i; j++) {
				for (int k = 0; k < n - i; k++) {
					Qt_t[j][k] -= 2 * w[j]*w[k];
				}
			}
			vector<vector<double>> R_t;
			for (int j = 0; j < R[R.size() - 1].size()-1; j++) {
				vector<double> t(R[R.size() - 1].size() - 1);
				R_t.push_back(t);
			}
			for (int j = i; j < R[R.size() - 1].size(); j++) {
				for (int k = i; k < R[R.size() - 1].size(); k++) {
					R_t[j - i][k - i] = R[R.size() - 1][j][k];
				}
			}
			vector<vector<double>> Rt_t = matrix_mult(Qt_t, R_t);
			
			vector<vector<double>> R_temp;
			vector<vector<double>> Q_temp;

			// заполним сначала Q
			for (int j = 0; j < n; j++) {
				vector<double> t(n);
				Q_temp.push_back(t);
				R_temp.push_back(t);
			}
			for (int j = 0; j < i; j++) {
				Q_temp[j][j] = 1;
			}
			for (int j = i; j < n; j++) {
				for (int k = i; k < n; k++) {
					Q_temp[j][k] = Qt_t[j - i][k - i];
				}
			}

			// теперь заполняем R
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < n; k++) {
					R_temp[j][k] = R[R.size() - 1][j][k];
					R_temp[k][j] = R[R.size() - 1][k][j];
				}
			}
			for (int j = i; j < n; j++) {
				for (int k = i; k < n; k++) {
					R_temp[j][k] = Rt_t[j - i][k - i];
				}
			}
			R.push_back(R_temp);
			Q.push_back(Q_temp);
		}
	}

	// итоговая матрица Q - Q_itog
	vector<vector<double>> Q_T = matrix_mult(Q[1], Q[0]);
	for (int i = 2; i < Q.size();i++) {
		vector<vector<double>> Q_temp = matrix_mult(Q[i], Q_T);
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				Q_T[j][k] = Q_temp[j][k];
			}
		}
	}
	// транспонируем матрицу
	vector<vector<double>> Q_itog;
	for (int i = 0; i < n; i++) {
		vector<double> w(n);
		Q_itog.push_back(w);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q_itog[i][j] = Q_T[j][i];
		}
	}

	vector<vector<double>> R_itog = R[R.size() - 1];
	return make_pair(Q_itog, R_itog);
}

// функция для нахождения определителя
int determinant(vector<vector<double>> matrix) {
	if (matrix.size() == 1) {
		return matrix[0][0];
	}
	else {
		if (matrix.size() == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		}
		else {
			int det = 0;
			for (int col = 0; col < matrix.size(); col++) {
				det += pow(-1, col) * matrix[0][col] * determinant(minor(matrix, 0, col));
			}
			return det;
		}
	}
}

// нахождение обратной матрицы
vector<vector<double>> inverse_matrix(vector<vector<double>> matrix) {
	int n = matrix.size();
	vector<vector<double>> res;
	for (int i = 0; i < n; i++) {
		vector<double> temp(n);
		res.push_back(temp);
	}

	double obr_det = (double)1 / (double)determinant(matrix);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int alg_dop = determinant(minor(matrix, i, j));
			res[j][i] = pow(-1, i + j) * alg_dop * obr_det;
		}
	}
	return res;
}

// генерация случайных чисел от 1 до 20
int generate_random1(void) {
	random_device rd;
	mt19937 gen(rd());
	int min = 1;
	int max = 20;
	uniform_int_distribution<> dis(min, max);
	return dis(gen);
}

// генерация случайных чисел от 1 до 5
int generate_random2(void) {
	random_device rd;
	mt19937 gen(rd());
	int min = 1;
	int max = 5;
	uniform_int_distribution<> dis(min, max);
	return dis(gen);
}

// перемножение матрицы и вектора
vector<double> matrix_vector_mult(vector<vector<double>> m, vector<double> v) {
	int n = v.size();
	vector<double> res(n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res[i] += m[i][j] * v[j];
		}
	}
	return res;
}

// функция для генерирования матрицы А через Lambda, C
vector<vector<double>> generate_A(int n) {
	// generation matrix A
	vector<vector<double>> Lambda;
	for (int i = 0; i < n; i++) {
		vector<double> temp(n);
		temp[i] = generate_random1();
		Lambda.push_back(temp);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << Lambda[i][j] << " ";
		}
		cout << '\n';
	}
	cout << '\n';

	// generating C, det != 0
	vector<vector<double>> C;
	for (int i = 0; i < n; i++) {
		vector<double> temp(n);
		C.push_back(temp);
	}

	while (determinant(C) == 0) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = generate_random1();
			}
		}
	}
	vector<vector<double>> inv_C = inverse_matrix(C);

	vector<vector<double>> A = matrix_mult(matrix_mult(inv_C, Lambda), C);
	return A;
}

// степенной метод: возращает макс.по модулю с.ч. и соотв.с.в.
pair<double, vector<double>> step_method(vector<vector<double>> A) {
	// размерность матрицы, ограничения
	int n = A.size();
	double delta = pow(10, -8);
	double rtol = pow(10, -6);

	vector<vector<double>> y;
	vector<vector<double>> z;
	vector<vector<double>> lambda;
	vector<double> s_vector;
	double s_numb = 0;

	// генерируем y0
	vector<double> y0(n);
	for (int i = 0; i < n; i++) {
		y0[i] = generate_random2();
	}
	y.push_back(y0);

	// нормируем и получаем z0
	double y0_norm = norm_vector(y0);
	vector<double> z0(n);
	for (int i = 0; i < n; i++) {
		z0[i] = y0[i] / y0_norm;
	}
	z.push_back(z0);

	// итерации
	for (int i = 0; i < 100; i++) {
		// получили следующий вектор y - перемножаем A и последний полученный z 
		vector<double> y_next = matrix_vector_mult(A, z[z.size() - 1]);
		y.push_back(y_next);

		// нормируем и получаем z_next
		vector<double> z_next(n);
		double y_next_norm = norm_vector(y_next);
		for (int j = 0; j < n; j++) {
			z_next[j] = y_next[j] / y_next_norm;
		}
		z.push_back(z_next);

		// получаем лямбда_i: либо делим y_i (k) / z_i (k-1), либо 0-ая координата
		vector<double> l_next(n);
		for (int j = 0; j < n; j++) {
			if (abs(z[z.size() - 2][j]) > delta) {
				l_next[j] = y[y.size() - 1][j] / z[z.size() - 2][j];
			}
			else {
				l_next[j] = 0;
			}
		}

		lambda.push_back(l_next);

		// начиная со 2 итерации, начинаем проверять лямбды
		if (i >= 1) {
			// вектор разности лямбд
			vector<double> razn_l(n);
			for (int j = 0; j < n; j++) {
				razn_l[j] = lambda[lambda.size() - 1][j] - lambda[lambda.size() - 2][j];
			}

			// нормируем вектор разности
			double norm_lambda = norm_vector(razn_l);
			// если условие выполнено, то последниее z - собственный вектор
			if (norm_lambda <= rtol * max(norm_vector(lambda[lambda.size() - 1]), norm_vector(lambda[lambda.size() - 2]))) {
				for (int j = 0; j < n; j++) {
					s_vector = z[z.size() - 1];
				}
				break;
			}
		}
	}
	// считаем мощность мн-ва I, суммируем лямбда_i (k)
	int I = 0;
	double summa = 0;
	for (int i = 0; i < n; i++) {
		if (abs(lambda[lambda.size() - 1][i]) > delta) {
			I++;
			summa += lambda[lambda.size() - 1][i];
		}
	}

	// присваиваем значение в собственное число и собственный вектор
	s_numb = summa / I;
	s_vector = z[z.size() - 1];
	return make_pair(s_numb, s_vector);
}

// QR алгоритм для нахождения всех с.ч.
vector<double> qr(vector<vector<double>> A) {
	// ограничения eps, размерность матрицы, инициализируем вектор собственных чисел
	double eps = pow(10, -8);
	int n = A.size();
	vector<double> sobs_num;

	// переносим А в work
	vector<vector<double>> work;
	for (int i = 0; i < n; i++) {
		work.push_back(A[i]);
	}

	// строим матрицу Хессенберга
	for (int i = 0; i < n - 2; i++) {
		// формируем si
		double summa = 0;
		double si;
		for (int j = i + 1; j < n; j++) {
			summa += pow(work[j][i], 2);
		}
		if (work[i + 1][i] >= 0) {
			si = pow(summa, 0.5);
		}
		else {
			si = -pow(summa, 0.5);
		}

		// формируем мю i-ые
		double mui = 1 / (double)pow(2 * si * (si - work[i + 1][i]), 0.5);

		// вектор v
		vector<double> v(n);
		v[i + 1] = mui * (work[i + 1][i] - si);
		for (int j = i + 2; j < n; j++) {
			v[j] = mui * work[j][i];
		}

		// формируем матрицу Hi
		vector<vector<double>> Hi;
		for (int j = 0; j < n; j++) {
			vector<double> temp(n);
			temp[j] = 1;
			Hi.push_back(temp);
		}

		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				Hi[j][k] -= 2 * v[j] * v[k];
			}
		}
		// промежуточный этап матрицы Хессенберга - перемножение Hi, work(A), Hi
		vector<vector<double>> Hessen = matrix_mult(matrix_mult(Hi, work), Hi);
		// переносим Hessen в work
		work = Hessen;
	}

	// вектор векторов итераций - добавляем work (матрица Хессенберга)
	vector<vector<vector<double>>> iterations;
	iterations.push_back(work);

	// начинаются итерации
	for (int q = 0; q < 100; q++) {
		// определяем матрицу, на основе которой будем строить новую
		vector<vector<double>> B = iterations[iterations.size() - 1];

		// инициализируем новую и заполняем
		vector<vector<double>> C; // B(k) - b(nn)*E = C
		for (int i = 0; i < B.size(); i++) {
			vector<double> p = B[i];
			p[i] -= B[B.size()-1][B.size() - 1];
			C.push_back(p);
		}

		// матрицу С раскладываем в QR
		vector<vector<double>> Q_C = QR_razl(C).first;
		vector<vector<double>> R_C = QR_razl(C).second;

		// строим матрицу D = R_C * Q_C + b(nn) * E = B(k+1)
		vector<vector<double>> D = matrix_mult(R_C, Q_C);
		for (int i = 0; i < D.size(); i++) {
			D[i][i] += B[B.size() - 1][B.size() - 1];
		}
		
		// следим за поддиагональным элементом, пока он не стал меньше eps
		// когда стал, добавляем диаг.эл-т в вектор собст.чисел и понижаем размерность
		if (abs(D[D.size() - 1][D.size() - 2]) < eps) {
			if (D.size() > 2) {
				sobs_num.push_back(D[D.size() - 1][D.size() - 1]);
				vector<vector<double>> D_temp;
				for (int i = 0; i < D.size() - 1; i++) {
					vector<double> l(D.size() - 1);
					D_temp.push_back(l);
				}
				for (int i = 0; i < D_temp.size(); i++) {
					for (int j = 0; j < D_temp.size(); j++) {
						D_temp[i][j] = D[i][j];
					}
				}
				iterations.push_back(D_temp);
			}
			else {
				sobs_num.push_back(D[1][1]);
				sobs_num.push_back(D[0][0]);
			
				break;
			}
		}
		else {
			iterations.push_back(D);
		}
	}
	return sobs_num;
}

// обратный степенной метод для определения всех с.ч. и собственных векторов
pair<vector<double>, vector<vector<double>>> obr_stepen(vector<vector<double>> A) {
	// размерность и ограничения
	int n = A.size();
	double delta = pow(10, -8);
	double rtol = pow(10, -6);

	// векторы собстыенных чисел и собственных векторов
	vector<double> res_snum(n);
	vector<vector<double>> res_svector;

	// собственные числа (нужно, чтобы получить оценки сигма)
	vector<double> s_nums = qr(A);

	for (int q = 0; q < n; q++) {
		double curr_s = s_nums[q];

		// последовательность y, z, sigma для этого с.ч.
		vector<vector<double>> y;
		vector<vector<double>> z;
		vector<vector<double>> mu;
		vector<double> sigma;

		// выбираем начальный сдвиг для curr_s
		double s;
		while (true) {
			s = curr_s + pow(generate_random2(), -2);
			int k = 0;
			for (int j = 0; j < n; j++) {
				if (q != j) {
					if (abs(curr_s - s) < abs(s_nums[j] - s)) {
						k++;
					}
				}
			}
			if (k == n - 1) {
				break;
			}
		}
		sigma.push_back(s);

		// теперь генерируем y0 и z0
		vector<double> y0(n);
		for (int i = 0; i < n; i++) {
			y0[i] = generate_random2();
		}
		y.push_back(y0);

		vector<double> z0(n);
		double norm_y0 = norm_vector(y0);
		for (int i = 0; i < n; i++) {
			z0[i] = y0[i] / norm_y0;
		}
		z.push_back(z0);

		// начинаем итерации
		for (int w = 0; w < 100; w++) {
			// выбираем current sigma
			double curr_sigma = sigma[sigma.size() - 1];

			// матрица A с минус sigma * E
			vector<vector<double>> curr_A;
			for (int i = 0; i < n; i++) {
				curr_A.push_back(A[i]);
				curr_A[curr_A.size() - 1][i] -= curr_sigma;
			}

			// получаем y(k) решаем систему Гауссом
			vector<double> y_curr = gauss1(curr_A, z[z.size() - 1]);
			y.push_back(y_curr);

			// получаем z
			double norm_y = norm_vector(y_curr);
			vector<double> z_curr(n);
			for (int j = 0; j < n; j++) {
				z_curr[j] = y_curr[j] / norm_y;
			}
			z.push_back(z_curr);

			// выбираем y, z для формирования mu
			vector<double> work_y = y[y.size() - 1];
			vector<double> work_z = z[z.size() - 2];
			vector<double> curr_mu(n);
			for (int j = 0; j < n; j++) {
				if (abs(work_y[j]) > delta) {
					curr_mu[j] = work_z[j] / work_y[j];
				}
				else {
					curr_mu[j] = 0;
				}
			}
			mu.push_back(curr_mu);

			// новый сдвиг
			double new_s = curr_sigma;
			double summa = 0;
			int I = 0;
			for (int i = 0; i < n; i++) {
				if (abs(curr_mu[i]) > delta) {
					I++;
					summa += curr_mu[i];
				}
			}
			new_s += summa / I;

			// проверяем сходимость и получаем собс.число через сигму и собственный вектор - последний z
			if (abs(curr_s - new_s) < pow(10, -6)) {
				res_snum[q] = new_s;
				res_svector.push_back(z[z.size() - 1]);
				break;
			}

		/*	if (w >= 1) {

				vector<double> check_mu(n);
				for (int i = 0; i < n; i++) {
					check_mu[i] = mu[mu.size() - 1][i] - mu[mu.size() - 2][i];
				}
				double norm_check = norm_vector(check_mu);
				if (norm_check <= rtol * max(norm_vector(mu[mu.size() - 1]), norm_vector(mu[mu.size() - 2]))) {
					res_snum[q] = sigma[sigma.size() - 1];
					res_svector.push_back(z[z.size() - 1]);
					break;
				}
			}*/
		}
	}

	return make_pair(res_snum, res_svector);
}

void main() {
	int n;
	cout << "Matrix dimension " << '\n';
	cin >> n;
	
	vector<vector<double>> A;
	// generating matrix A
	while (true) {
		A = generate_A(n);
		if (qr(A).size() == n) {
			vector<double> all_nums = qr(A);
			int smth = 0;
			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {
					if (to_string(all_nums[i]) == to_string(all_nums[j])) {
						smth++;
						break;
					}
				}
			}
			if (smth == 0) {
				break;
			}
		}
	}

	// printing matrix A
	cout << "matrix A" << '\n';
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << '\n';
	}

	//string method;
	/*cin >> method;*/

	//if (method == "stepen") {
	//	pair<double, vector<double>> res = step_method(A);
	//	cout << "s_number " << res.first << '\n';
	//	cout << "s_vector" << " ";
	//	for (int i = 0; i < n; i++) {
	//		cout << res.second[i] << " ";
	//	}
	//	cout << '\n';
	//}
	/*if (method == "qr") {
		vector<double> res = qr(A);
		for (int i = 0; i < res.size(); i++) {
			cout << setprecision(5) << res[i] << " ";
		}
	}*/
	/*if (method == "obr_stepen") {
		pair<vector<double>, vector<vector<double>>> obr = obr_stepen(A);
		for (int i = 0; i < n; i++) {
			cout << "s_num " << obr.first[i] << '\n';
			cout << "s_vector ";
			for (int j = 0; j < n; j++) {
				cout << obr.second[i][j] << " ";
			}
			cout << '\n';
		}
	}*/


	// printing results of all methods
	cout <<'\n' << "stepen_method" << '\n';
	pair<double, vector<double>> res = step_method(A);
	cout << "s_number " << res.first << '\n';
	cout << "s_vector" << " ";
	for (int i = 0; i < res.second.size(); i++) {
		cout << res.second[i] << " ";
	}
	cout << '\n' << '\n';

	cout << "qr" << '\n';
	vector<double> res1 = qr(A);
	for (int i = 0; i < res1.size(); i++) {
		cout << res1[i] << " ";
	}
	cout << '\n' << '\n';

	cout << "obr_stepen" << '\n';
	pair<vector<double>, vector<vector<double>>> obr = obr_stepen(A);
	for (int i = 0; i < n; i++) {
		cout << "s_num " << obr.first[i] << '\n';
		cout << "s_vector ";
		for (int j = 0; j < n; j++) {
			cout << obr.second[i][j] << " ";
		}
		cout << '\n';
	}
} 