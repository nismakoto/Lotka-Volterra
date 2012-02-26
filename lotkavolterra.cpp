#include <iostream>
#include <fstream>
using namespace std;

// 被食者の個体数の初期値
#define X_init 1.0
// 捕食者の個体数の初期値
#define Y_init 1.0
// 刻み幅
#define H 0.01
// ロトカ・ボルテラ方程式の各パラメータ
#define Alpha 5.2
#define Beta 3.4
#define Gamma 2.1
#define Delta 1.4
// 初期時間
#define T_start 0.0
// 終了時間
#define T_finish 10

// 被食者の微分方程式
double fx(double x, double y)
{
        return(Alpha*x-Beta*x*y);
}

// 捕食者の微分方程式
double fy(double x, double y)
{
        return(-Gamma*y+Delta*x*y);
}

// ルンゲ・クッタ法
void runge_kutta (double *x, double *y, double h, double (*fx)(double, double), double (*fy)(double, double))
{
        double dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4;

        dx1 = fx(*x, *y);
        dy1 = fy(*x, *y);
        dx2 = fx(*x+(h/2.0)*dx1, *y+(h/2.0)*dy1);
        dy2 = fy(*x+(h/2.0)*dx1, *y+(h/2.0)*dy1);
        dx3 = fx(*x+(h/2.0)*dx2, *y+(h/2.0)*dy2);
        dy3 = fy(*x+(h/2.0)*dx2, *y+(h/2.0)*dy2);
        dx4 = fx(*x+h*dx3, *y+h*dy3);
        dy4 = fy(*x+h*dx3, *y+h*dy3);

        *x += h * (dx1 + 2.0*dx2 + 2.0*dx3 + dx4) / 6.0;
        *y += h * (dy1 + 2.0*dy2 + 2.0*dy3 + dy4) / 6.0;
}

int main()
{
        double x, y;
        // 初期値の代入
        double i = T_start;
        x = X_init;
        y = Y_init;

        ofstream out("lotkavolterra.dat");

        // エラー処理
        if(!out)
        {
                cout << "Cannot open file.\n";
        }

        do{
                // ファイルに解を出力
                out << i << " " << x << " " << y << endl;

                // ルンゲ・クッタ法を計算する関数の引数に被食者と捕食者の関数及び変数の
                // ポインタを渡す
                runge_kutta(&x, &y, H, fx, fy);

                i += H;

        } while(i <= T_finish);

        // ファイルを閉じる
        out.close();

        return 0;
}
