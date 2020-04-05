#include <iostream>
#include <vector>

using namespace std;

vector<double>& operator+=(vector<double>& left, const vector<double>& right) {  //сумма векторов
	for (int i = 0; i < right.size(); i++)left[i] = left[i] + right[i];
	return left;
}

vector<double>& operator-=(vector<double>& left, const vector<double>& right) {  //разность вектор
	for (int i = 0; i < right.size(); i++)left[i] = left[i] - right[i];
	return left;
}

double scalar_product(vector<double> a, vector<double>b, int n){
	double c = 0.0;
	for (int i = 0; i < n; i++)
		c += a[i] * b[i];
	return c;
}
vector<double>  vec_proj(vector<double> a, vector<double> b, int n)
{
	double k = scalar_product(a, b,n) / scalar_product(b, b,n);
	for (int i = 0; i < n; i++) b[i] *= k;
	return b;
}
void mult(vector <vector <double>> A, vector <vector <double>> B,	//matrix*matrix
	vector <vector <double>> &R, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}

void mult(vector <vector <double>> A, vector <double> b,     //matrix*vector
	vector <double> &R, int n)
{
	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
			R[i] += A[i][k] * b[k];
}

void QR(vector <vector <double>> A, vector <vector <double>> &Q,
	vector <vector <double>> &R, int n)
{
	vector <vector <double>> Q_t(n);//Q transp
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)Q_t[j].push_back(0);
	}
	vector<double> proj_sum(n); // текущая сумма векторов проекций
	for (int i = 0; i < n; i++)proj_sum[i] = 0;
	
	for (int i = 0; i < n; i++) {	//транспонирование, так как нужно  разбиение по столбцам
		for (int j = 0; j < n; j++) {
			Q_t[i][j] = A[j][i];
			Q[i][j] = A[j][i];
		}
	}


	for (int i = 1; i < n; i++) {//вычисление ортогональных векторов		
		for (int k = 0; k < i; k++) {
				proj_sum += vec_proj(Q[i], Q_t[k],n);
		}
		Q_t[i] -= proj_sum; 
		proj_sum -= proj_sum;//proj_sum = 0;
	}

	
	//нормализация:
	double tmp; 
	for (int i = 0; i < n; i++) {
		tmp = scalar_product(Q_t[i], Q_t[i], n);
		for (int j = 0; j < n; j++) {		
			Q_t[i][j] /= sqrt(tmp);
		}
	}

	for (int i = 0; i < n; i++) {	// получение Q:
		for (int j = 0; j < n; j++) {
			Q[i][j] = Q_t[j][i];
		}
	}

	//R=Q_t*A
	mult(Q_t, A, R, n);

}

void slae_solution(vector <vector <double>> &Q,
	vector <vector <double>> &R, vector<double> &b, vector<double> &x, int n) {
	
	vector<double> y(n);
	vector <vector <double>> Res(n);//Rx=Q_t*b;
	vector <vector <double>> Q_t(n);//Q transp
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q_t[j].push_back(0);
			Res[j].push_back(0);
		}
	}
	for (int i = 0; i < n; i++) {	// получение Q transp :
		for (int j = 0; j < n; j++) {
			Q_t[i][j] = Q[j][i];
		}
	}
	
	mult(Q_t, b, y, n);//b to matrix and y to vector
	
	//обратный ход(Rx=Q_t*b):
	for (int i = 0; i < n - 1; i++)x[i] = 0;
	x[n - 1] = y[n - 1] / R[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= R[i][j] * x[j];
		}
		x[i] = (x[i] + y[i]) / R[i][i];
	}
}

void show(vector <vector <double>> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << A[i][j] << "\t";
		}
		cout << endl;
	}
}

int main()
{
	const double n = 4;
	double det = 1.0;
	vector <vector <double>> A(n), Q(n), R(n), Res(n),TMP(n);
	vector <double> b(n);
	vector <double> y(n), x(n);
	for (int i = 0; i < n; i++)
	{
		b[i] = rand() % 20;
		for (int j = 0; j < n; j++)
		{
			A[i].push_back(rand() % 20);
			Q[i].push_back(0);
			R[i].push_back(0);
			Res[i].push_back(0);
		}
	}
	QR(A, Q, R, n);

	cout << "b:" << endl;
	for (int i = 0; i < n; i++)cout << b[i]<<" ";
	cout << "Fisrt matrix" << endl;
	show(A, n);
	cout << "R matrix" << endl;
	show(R, n);
	cout << "Q matrix" << endl;
	show(Q, n);
	mult(Q, R, Res, n);
	cout << "Q*R matrix" << endl;
	show(Res, n);
	slae_solution(Q, R, b, x, n);
	cout << "slae solution:" << endl;
	for (int i = 0; i < n; i++)cout << x[i] << " ";
	
	system("pause");
	return 0;
}