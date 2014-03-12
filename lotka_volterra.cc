#include <iostream>
#include <fstream>
using namespace std;

// Prey
#define X_init 1.0
// Predator
#define Y_init 1.0
// Step size
#define H 0.01
// Parameters
#define Alpha 5.2
#define Beta 3.4
#define Gamma 2.1
#define Delta 1.4
// Start time
#define T_start 0.0
// Finish time
#define T_finish 10

// Prey
double fx(double x, double y) { return (Alpha * x - Beta * x * y); }

// Predator
double fy(double x, double y) { return (-Gamma * y + Delta * x * y); }

// Runge-Kutta function
void runge_kutta(double *x, double *y, double h, double (*fx)(double, double),
                 double (*fy)(double, double)) {
  double dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4;

  dx1 = fx(*x, *y);
  dy1 = fy(*x, *y);
  dx2 = fx(*x + (h / 2.0) * dx1, *y + (h / 2.0) * dy1);
  dy2 = fy(*x + (h / 2.0) * dx1, *y + (h / 2.0) * dy1);
  dx3 = fx(*x + (h / 2.0) * dx2, *y + (h / 2.0) * dy2);
  dy3 = fy(*x + (h / 2.0) * dx2, *y + (h / 2.0) * dy2);
  dx4 = fx(*x + h * dx3, *y + h * dy3);
  dy4 = fy(*x + h * dx3, *y + h * dy3);

  *x += h * (dx1 + 2.0 * dx2 + 2.0 * dx3 + dx4) / 6.0;
  *y += h * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
}

int main() {
  double x, y;
  double i = T_start;
  x = X_init;
  y = Y_init;

  ofstream out("lotka_volterra.dat");

  if (!out) {
    cout << "Cannot open file.\n";
  }

  do {
    out << i << " " << x << " " << y << endl;
    runge_kutta(&x, &y, H, fx, fy);

    i += H;

  } while (i <= T_finish);

  out.close();

  return 0;
}
