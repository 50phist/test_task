#include<iostream>

void fnc(const double x[], double* dxdt, const int len)
{
	// Вектор-функция правых частей ОДУ
	for (int i = 0; i < len; i++)
	{
		switch (i)
		{
		case 0: // dwxdt
			dxdt[i] = x[i] + 2;
			break;
		case 1: // dwydt
			dxdt[i] = x[i] + 2;
			break;
		case 2: // dwzdt
			dxdt[i] = x[i] + 2;
			break;
		}
	}
}

int main()
{
	double x[]{1, 1, 1, 1};
	int len = sizeof(x) / sizeof(x[0]);
	double *dxdt = new double[len];

	fnc(x, dxdt, len);
	std::cout << "++++++++" << std::endl;
	for (int i = 0; i < len; i++)
	{
		std::cout << dxdt[i] << std::endl;
	}

}
