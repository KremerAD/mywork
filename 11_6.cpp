#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#pragma warning(disable:4996)
using namespace std;

#define EPS 1e-14
#define PI 3.14159265358979323846
#define N_x 4 //number of partitions
#define N_y 4
#define a_x 1. //interval for axis x
#define b_x 4.
#define a_y 1. //interval for axis y
#define b_y 4.
#define n_f 3 //number of basis func

double f(double x, double y) {
    //return y + 5.;
    //return exp(-(x * x + y * y));
    //return pow(y, 2) * pow(1. - y, 2);
    return 1. / 3. * pow(y, 3) - 5. / 2. * pow(y, 2) + 4. * y;
}

void initial(vector < vector<double> >& matrix_real, vector < vector<double> >& Ind, double h_x, double h_y);
void approx(double h_x, vector < vector<double> > matrix_real, vector < vector<double> > Ind, vector < vector<double> >& matrix_approx);
void splain(double h_y, vector < vector<double> > matrix_real, vector < vector<double> > Ind, vector < vector<double> >& matrix_splain);
void write_in(ofstream& f, vector < vector<double> > matrix);
void renew(vector < vector<double> > matrix_splain, vector < vector<double> > matrix_approx, vector < vector<double> >& matrix_renew);
void error(vector < vector<double> > matrix_real, vector < vector<double> > matrix_renew, vector < vector<double> >& matrix_error);
void trans(vector < vector<double> > matrix, vector < vector<double> >& matrix_t); // transpose matrix
void multi(vector < vector<double> > a, vector < vector<double> > b, vector < vector<double> >& ab); //multiply matrix
void gauss(vector < vector<double> > a, vector<double>& ans); //method Gauss
double fumat(int i, int j, vector <double> x);

int main()
{
    ofstream fout1("out.txt");
    double h_x = (b_x - a_x) / (N_x - 1.), h_y = (b_y - a_y) / (N_y - 1.);
    vector < vector<double> > matrix_real, matrix_splain, matrix_approx, matrix_renew, matrix_error, Ind;
    matrix_real.resize(N_y);
    matrix_splain.resize(N_y);
    matrix_approx.resize(N_y);
    matrix_renew.resize(N_y);
    matrix_error.resize(N_y);
    Ind.resize(N_y);

    initial(matrix_real, Ind, h_x, h_y);
    fout1 << "Matrix of indicators:" << endl;
    write_in(fout1, Ind);
    fout1 << "Matrix real: " << endl;
    write_in(fout1, matrix_real);

    splain(h_y, matrix_real, Ind, matrix_splain);
    fout1 << "Matrix after splain:" << endl;
    write_in(fout1, matrix_splain);

    approx(h_x, matrix_real, Ind, matrix_approx);
    fout1 << "Matrix after approx:" << endl;
    write_in(fout1, matrix_approx);

    renew(matrix_splain, matrix_approx, matrix_renew);
    fout1 << "Matrix renew:" << endl;
    write_in(fout1, matrix_renew);

    error(matrix_real, matrix_renew, matrix_error);
    fout1 << "Matrix of errors:" << endl;
    write_in(fout1, matrix_error);

    fout1.close();
    return 0;
}

void initial(vector < vector<double> >& matrix_real, vector < vector<double> >& Ind, double h_x, double h_y) {
    //int k = 0;
    for (int i = 0; i < N_y; i++) {
        for (int j = 0; j < N_x; j++) {
            matrix_real[i].push_back(f(a_x + j * h_x, a_y + i * h_y));
            if (i == 2) {
                if ((j == 0) || ((j == 3))) {
                    Ind[i].push_back(2);
                }
                else {
                    Ind[i].push_back(0);
                }
            }
            if (i == 0) {
                Ind[i].push_back(2);
            }
            if ((i != 0) && (i != 2)) {
                Ind[i].push_back(2);
            }
            //Ind[i].push_back(rand() % 3);
            

            //k++;

        }
    }
}

double fumat(int i, int j, vector <double> x) {
    return pow(x[i], j);
}

void gauss(vector < vector<double> > a, vector<double>& ans) {
    int n = (int)a.size();
    int m = (int)a[0].size() - 1;

    vector<int> where(m, -1);
    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i)
            if (abs(a[i][col]) > abs(a[sel][col]))
                sel = i;
        if (abs(a[sel][col]) < EPS)
            continue;
        for (int i = col; i <= m; ++i)
            swap(a[sel][i], a[row][i]);
        where[col] = row;

        for (int i = 0; i < n; ++i)
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (int j = col; j <= m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }
    ans.assign(m, 0);
    for (int i = 0; i < m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
}

void write_in(ofstream& f, vector < vector<double> > matrix) {
    int i_max = matrix.size();
    int j_max = matrix[0].size();
    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            f << left << setw(13) << matrix[i][j] << " ";
        }
        f << endl;
    }
    f << endl;
}

void splain(double h_y, vector < vector<double> > matrix_real, vector < vector<double> > Ind, vector < vector<double> >& matrix_splain) {
    for (int j = 0; j < N_x; j++) {
        int N = 0;
        vector <double> d, x, y, M;
        for (int i = 0; i < N_y; i++) {
            if (fabs(Ind[i][j]) < 0.5) {
                matrix_splain[i].push_back(EPS);
            }
            else {
                matrix_splain[i].push_back(matrix_real[i][j]);
                y.push_back(matrix_real[i][j]);
                x.push_back(a_y + i * h_y);
                N++;
            }
        }
        vector < vector<double> > c;
        c.resize(N - 2);

        for (int n = 0; n < N - 2; n++) {
            for (int m = 0; m < N - 2; m++) {
                if (n == m) {
                    c[n].push_back((x[n + 2] - x[n]) / 3.);
                }
                else if (m == n + 1) {
                    c[n].push_back((x[m + 1] - x[m]) / 6.);
                }
                else if (n == m + 1) {
                    c[n].push_back((x[n + 1] - x[n]) / 6.);
                }
                else {
                    c[n].push_back(0.);
                }
            }
        }

        for (int i = 0; i < N - 2; i++)
            c[i].push_back((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - (y[i + 1] - y[i]) / (x[i + 1] - x[i]));
        cout << c[0][1] << endl;
        gauss(c, d);
        M.push_back(0.);
        for (int i = 0; i < N - 2; ++i) {
            M.push_back(d[i]);
        }
        M.push_back(0.);
        int temp = 0;
        for (int i = 0; i < N_y; i++) {
            if (Ind[i][j] < 0.5) {
                if (temp == N) temp -= 1;
                matrix_splain[i][j] = M[temp - 1] * pow(x[temp] - (a_y + i * h_y), 3) / (6.0 * (x[temp] - x[temp - 1])) +
                    M[temp] * pow((a_y + i * h_y) - x[temp - 1], 3) / (6.0 * (x[temp] - x[temp - 1])) +
                    (y[temp - 1] - (M[temp - 1] * pow(x[temp] - x[temp - 1], 2)) / 6.0) * (x[temp] - (a_y + i * h_y)) / (x[temp] - x[temp - 1]) +
                    (y[temp] - (M[temp] * pow(x[temp] - x[temp - 1], 2)) / 6.0) * (a_y + i * h_y - x[temp - 1]) / (x[temp] - x[temp - 1]);
            }
            else {
                temp++;
            }
        }
    }
}

void approx(double h_x, vector < vector<double> > matrix_real, vector < vector<double> > Ind, vector < vector<double> >& matrix_approx) {
    int N;
    for (int s = 0; s < N_y; s++) {
        N = 0;
        vector <double> x, x2, y, c;
        for (int j = 0; j < N_x; j++) {
            if (fabs(Ind[s][j]) < 0.5) {
                matrix_approx[s].push_back(0.);
                x2.push_back(a_x + j * h_x);
            }
            else {
                matrix_approx[s].push_back(matrix_real[s][j]);
                y.push_back(matrix_real[s][j]);
                x.push_back(a_x + j * h_x);
                N++;
            }
        }
        vector < vector<double> > a, at, ata, atb, matrix_g;
        a.resize(N);
        at.resize(n_f);
        ata.resize(n_f);
        atb.resize(n_f);
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < n_f; m++) {
                a[n].push_back(fumat(n, m, x));
            }
        }
        trans(a, at);
        multi(at, a, ata);
        double temp = 0;
        for (int n = 0; n < n_f; n++) {
            temp = 0;
            vector < vector <double> > vtemp, ytemp;
            vtemp.resize(N);
            ytemp.resize(N);
            for (int m = 0; m < N; m++) {
                ytemp[m].push_back(y[m]);
            }
            multi(at, ytemp, vtemp);
            atb[n].push_back(vtemp[n][0]);
        }
        for (int n = 0; n < n_f; n++) {
            ata[n].push_back(atb[n][0]);
        }
        gauss(ata, c);
        int temp0 = 0;
        for (int j = 0; j < N_x; j++) {
            if (Ind[s][j] < 0.5) {
                for (int k = 0; k < n_f; k++) {
                    matrix_approx[s][j] += c[k] * fumat(temp0, k, x2);
                }
                temp0++;
            }
        }
    }
}

void renew(vector < vector<double> > matrix_splain, vector < vector<double> > matrix_approx, vector < vector<double> >& matrix_renew) {
    for (int i = 0; i < N_y; i++) {
        for (int j = 0; j < N_x; j++) {
            matrix_renew[i].push_back((matrix_splain[i][j] + matrix_approx[i][j]) / 2);
        }
    }
}

void error(vector < vector<double> > matrix_real, vector < vector<double> > matrix_renew, vector < vector<double> >& matrix_error) {
    for (int i = 0; i < N_y; i++) {
        for (int j = 0; j < N_x; j++) {
            matrix_error[i].push_back(fabs(matrix_real[i][j] - matrix_renew[i][j]));
        }
    }
}

void trans(vector < vector<double> > matrix, vector < vector<double> >& matrix_t) {
    int n = (int)matrix.size();
    int m = (int)matrix[0].size();
    matrix_t.resize(m);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            matrix_t[j].push_back(matrix[i][j]);
        }
    }
}

void multi(vector < vector<double> > a, vector < vector<double> > b, vector < vector<double> >& ab) {
    int n_1 = (int)a.size();
    int m_1 = (int)a[0].size();
    int n_2 = (int)b.size();
    int m_2 = (int)b[0].size();
    ab.resize(n_1);
    for (int i = 0; i < n_1; i++) {
        for (int j = 0; j < m_2; j++) {
            double temp = 0;
            for (int k = 0; k < m_1; k++) {
                temp += a[i][k] * b[k][j];
            }
            ab[i].push_back(temp);
        }
    }
}