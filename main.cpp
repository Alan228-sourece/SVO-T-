#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double dostresh(double x) {
    return x * x * x + cos(x);
}

double f(double x) {
    return -6 * x + cos(x);
}

vector<double> mgrid(int n) {
    vector<double> x(n + 1);
    double h = 1.0 / n;
    for (int i = 0; i <= n; ++i) {
        x[i] = i * h;
    }
    return x;
}

void fillM(vector<vector<double>>& A, vector<double>& F, const vector<double>& x, double a, double b, int n) {
    double h = 1.0 / n;
    for (int i = 0; i < n - 1; ++i) {
        if (i == 0) {
            A[i][i] = 2.0 / (h * h);
            A[i][i + 1] = -1.0 / (h * h);
            F[i] = f(x[i + 1]) + a / (h * h);
        }
        else if (i == n - 2) {
            A[i][i - 1] = -1.0 / (h * h);
            A[i][i] = 2.0 / (h * h);
            F[i] = f(x[i + 1]) + b / (h * h);
        }
        else {
            A[i][i - 1] = -1.0 / (h * h);
            A[i][i] = 2.0 / (h * h);
            A[i][i + 1] = -1.0 / (h * h);
            F[i] = f(x[i + 1]);
        }
    }
}

vector<double> progonka(const vector<vector<double>>& A, const vector<double>& F, int n) {
    vector<double> alpha(n - 1, 0.0);
    vector<double> beta(n - 1, 0.0);
    vector<double> y(n + 1, 0.0);

    alpha[0] = A[0][1] / A[0][0];
    beta[0] = F[0] / A[0][0];

    for (int i = 1; i < n - 1; ++i) {
        double TihiyDen = A[i][i] - A[i][i - 1] * alpha[i - 1];
        alpha[i] = 0;
        if ((i < n - 2)) {
            alpha[i] = A[i][i + 1] / TihiyDen;
        }
        beta[i] = (F[i] - A[i][i - 1] * beta[i - 1]) / TihiyDen;
    }

    y[n - 1] = beta[n - 2];
    for (int i = n - 2; i > 0; --i) {
        y[i] = beta[i - 1] - alpha[i - 1] * y[i + 1];
    }

    return y;
}

void errors(const vector<double>& y, const vector<double>& x, int n, double& cerr, double& lerr) {
    cerr = 0.0;
    lerr = 0.0;

    for (int i = 0; i <= n; ++i) {
        double error = fabs(y[i] - dostresh(x[i]));
        cerr = max(cerr, error);
        lerr += error * error;
    }

    lerr = sqrt(lerr / (n + 1));
}
void solve(vector<double>& cv, vector<double>& lv, int n) {
    double a = dostresh(0.0);
    double b = dostresh(1.0);
    vector<double> x = mgrid(n);
    vector<vector<double>> A(n - 1, vector<double>(n - 1, 0.0));
    vector<double> F(n - 1, 0.0);
    fillM(A, F, x, a, b, n);
    vector<double> y = progonka(A, F, n);
    y[0] = a;
    y[n] = b;
    double cerr, lerr;
    errors(y, x, n, cerr, lerr);
    cv.push_back(cerr);
    lv.push_back(lerr);
    //cout << "C- ошибка: " << cerr << endl;
    //cout << "L2- ошибка: " << lerr << endl;

}
int main() {
    int n;
    vector<double> cv;
    vector<double> lv;
    for (n = 16; n <= 30000; n *= 8)
    {
        solve(cv, lv, n);
    }
    for (n = 16; n <= 30000; n *= 8)
    {
        cout << n << " ";
    }
    cout << endl;
    for (n = 0; n < cv.size(); n++)
    {
        cout << cv[n] << " ";
    }
    cout << endl;
    for (n = 0; n < lv.size(); n++)
    {
        cout << lv[n] << " ";
    }
    return 0;
}
