#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

double Trapec_M(double(*f)(double), double a, double b, double E);
double Simpson_M(double(*f)(double), double a, double b, double E);
double SKF(double(*f)(double, double), double a, double b, double c, double d, double E);
double integral_1(double x)
{
	return (1 + x * x) / (1 + x * x * x);
}

double integral_2(double x, double y)
{
	return x * x / (1 + y * y);
}

int main()
{
	setlocale(LC_ALL, "Rus");
	double(*P[1])(double x) = { integral_1 };
	double(*D[1])(double x, double y) = { integral_2 };
	double ab[1][2];
	double abcd[1][4];
	ab[0][0] = 3;
	ab[0][1] = 4.254;
	abcd[0][0] = 0;
	abcd[0][1] = 4.0;
	abcd[0][2] = 1.0;
	abcd[0][3] = 2.0;
	
	cout << "\n Функция (1+x^2)/(1+x^3): \n";
	cout << "Метод трапеции при epsilon = 1e-4: " << endl;
	double integral = Trapec_M(*(P[0]), ab[0][0], ab[0][1], 1e-4);
	cout << "Метод трапеции при epsilon = 1e-5: " << endl;
	integral = Trapec_M(*(P[0]), ab[0][0], ab[0][1], 1e-5);
	cout << "Метод Симпсона при epsilon = 1e-4: " << endl;
	integral = Simpson_M(*(P[0]), ab[0][0], ab[0][1], 1e-4);
	cout << "Метод Симпсона при epsilon = 1e-5: " << endl;
	integral = Simpson_M(*(P[0]), ab[0][0], ab[0][1], 1e-5);
	cout << "\n Функция x^2/(1 + y^2): \n";
	cout << "Кубатурная формула Симпсона при E = 1e-6: " << endl;
	integral = SKF(*(D[0]), abcd[0][0], abcd[0][1], abcd[0][2], abcd[0][3], 1e-6);

	system("pause");
	return 0;
}

double Trapec_M(double (*f)(double), double a, double b, double E)
{

	int n = 2;
	double x = a;
	double h = (b - a) / n;
	double I_h = (f(a) + f(b)) / 2;
	int k = 0;
	for (int i = 1; i < n; i++)
	{
		x += h;
		I_h += f(x);
	}
	I_h *= h;
	k++;
	cout << "k= " << k << setw(10) << "h= " << h << setw(10) << "I= " << I_h << endl;

	n = 4;
	x = a;
	h = (b - a) / n;
	double I_h2 = (f(a) + f(b)) / 2;
	for (int i = 1; i < n; i++)
	{
		x += h;
		I_h2 += f(x);
	}
	I_h2 *= h;
	k++;
	cout << "k= " << k << setw(10) << "h= " << h << setw(15) << "I= " << I_h2 << endl;
	while (abs(I_h2 - I_h) > 3 * E)
	{
		I_h = I_h2;
		x = a;
		n *= 2;
		h = (b - a) / n;
		I_h2 = (f(a) + f(b)) / 2;
		for (int i = 1; i < n; i++)
		{
			x += h;
			I_h2 += f(x);
		}
		I_h2 *= h;
		k++;
		cout << "k= " << k << setw(10) << "h= " << h << setw(15) << "I= " << I_h2 << endl;
	}
	return I_h2;
}
double Simpson_M(double(*f)(double), double a, double b, double E)
{

	int n = 1;
	double x = a;
	double h = (b - a) / (2 * n);
	double I_h = f(x);
	int k = 0;
	for (int i = 1; i < 2 * n; i++)
	{
		x += h;
		if (i % 2 == 0)
		{
			I_h += 2 * f(x);
		}
		else
		{
			I_h += 4 * f(x);
		}
	}
	x += h;
	I_h += f(x);
	I_h *= h / 3;
	k++;
	cout << "k= " << k << setw(10) << "h= " << h << setw(10) << "I= " << I_h << endl;

	n = 2;
	x = a;
	h = (b - a) / (2 * n);
	double I_h2 = f(x);
	for (int i = 1; i < 2 * n; i++)
	{
		x += h;
		if (i % 2 == 0)
		{
			I_h2 += 2 * f(x);
		}
		else
		{
			I_h2 += 4 * f(x);
		}
	}
	x += h;
	I_h2 += f(x);
	I_h2 *= h / 3;
	k++;
	cout << "k= " << k << setw(10) << "h= " << h << setw(10) << "I= " << I_h2 << endl;
	while (abs(I_h2 - I_h) > 15 * E)
	{
		I_h = I_h2;
		n *= 2;
		x = a;
		h = (b - a) / (2 * n);
		I_h2 = f(x);
		for (int i = 1; i < 2 * n; i++)
		{
			x += h;
			if (i % 2 == 0)
			{
				I_h2 += 2 * f(x);
			}
			else
			{
				I_h2 += 4 * f(x);
			}
		}
		x += h;
		I_h2 += f(x);
		I_h2 *= h / 3;
		k++;
		cout << "k= " << k << setw(10) << "h= " << h << setw(10) << "I= " << I_h2 << endl;
	}
	return I_h2;
}
double SKF(double(*f)(double, double), double a, double b, double c, double d, double E)
{
	double ig_1, ig_2;
	ig_1 = 0;
	ig_2 = 0;
	int n = 2, m = 2;
	double hx = (b - a) / (2 * n),
		hy = (d - c) / (2 * m);
	int k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double x2 = a + (2 * i) * hx,
				y2 = c + (2 * j) * hy,
				x2_1 = x2 + hx,
				x2_2 = x2_1 + hx,
				y2_1 = y2 + hy,
				y2_2 = y2_1 + hy;
			ig_1 += f(x2, y2) + 4 * f(x2_1, y2) + f(x2_2, y2) + 4 * f(x2, y2_1) + 16 * f(x2_1, y2_1) + 4 * f(x2_2, y2_1) + f(x2, y2_2) + 4 * f(x2_1, y2_2) + f(x2_2, y2_2);
		}
	}
	ig_1 *= hx * hy / 9;
	k++;
	cout << "k= " << k << setw(10) << "hx= " << hx << setw(10) << "hy= " << hx << setw(10) << "I= " << ig_1 << endl;


	n *= 2;
	m *= 2;
	hx = (b - a) / (2 * n);
	hy = (d - c) / (2 * m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double x2 = a + (2 * i) * hx,
				y2 = c + (2 * j) * hy,
				x2_1 = x2 + hx,
				x2_2 = x2_1 + hx,
				y2_1 = y2 + hy,
				y2_2 = y2_1 + hy;
			ig_2 += f(x2, y2) + 4 * f(x2_1, y2) + f(x2_2, y2) + 4 * f(x2, y2_1) + 16 * f(x2_1, y2_1) + 4 * f(x2_2, y2_1) + f(x2, y2_2) + 4 * f(x2_1, y2_2) + f(x2_2, y2_2);
		}
	}
	ig_2 *= hx * hy / 9;
	k++;
	cout << "k= " << k << setw(10) << "hx= " << hx << setw(10) << "hy= " << hx << setw(10) << "I= " << ig_2 << endl;

	while (abs(ig_1 - ig_2) > 15 * E)
	{
		ig_1 = ig_2;
		ig_2 = 0;

		n *= 2;
		m *= 2;
		hx = (b - a) / (2 * n);
		hy = (d - c) / (2 * m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				double x2 = a + (2 * i) * hx,
					y2 = c + (2 * j) * hy,
					x2_1 = x2 + hx,
					x2_2 = x2_1 + hx,
					y2_1 = y2 + hy,
					y2_2 = y2_1 + hy;
				ig_2 += f(x2, y2) + 4 * f(x2_1, y2) + f(x2_2, y2) + 4 * f(x2, y2_1) + 16 * f(x2_1, y2_1) + 4 * f(x2_2, y2_1) + f(x2, y2_2) + 4 * f(x2_1, y2_2) + f(x2_2, y2_2);

			}
		}
			ig_2 *= hx * hy / 9;
			k++;
			cout << "k= " << k << setw(10) << "hx= " << hx << setw(10) << "hy= " << hx << setw(10) << "I= " << ig_2 << endl;	

	}
	return ig_1;
	
}
