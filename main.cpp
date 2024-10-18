#include<iostream>
#include<vector>

std::vector<double> SU(std::vector<double> qref, std::vector<double> x)
{
	/* Вектор-функция правой части системы ОДУ работы СУ */

	double Jr = 5e-4, Hmax = 0.2, Mmax = 5e-3;

	std::vector<double> dudt(3);
	for (int i = 0; i < dudt.size(); i++)
	{
		switch (i)
		{
		case 0: // dHxdt
			dudt[i] = x[i] + 2;
			break;
		case 1: // dHydt
			dudt[i] = x[i] + 2;
			break;
		case 2: // dHzdt
			dudt[i] = x[i] + 2;
		}
	}
	return dudt;
}

std::vector<double> KA(std::vector<double>(*SU)(std::vector<double> qref, std::vector<double> x), 
					   std::vector<double> qref, std::vector<double> x)
{
	/* Вектор-функция правой части системы ОДУ динамики КА */
	
	double J[3][3]{ { 3.0, -0.1, 0.0 }, { -0.1, 2.3, 0.2 }, { 0.0, 0.2, 1.9 } };
	double J_inv[3][3]{ { 0.333822, 0.0146481, -0.0015419 }, { 0.0146481, 0.439442, -0.046257 }, { -0.0015419, -0.046257, 0.531185 } };
	double Jr = 5e-4;

	std::vector<double> dudt = SU(qref, x);
	std::vector <double> dxdt(x.size());
	for (int i = 0; i < dxdt.size(); i++)
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
		case 3: // dwxrdt
			dxdt[i] = x[i] + 2;
			break;
		case 4: // dwyrdt
			dxdt[i] = x[i] + 2;
			break;
		case 5: // dwzrdt
			dxdt[i] = x[i] + 2;
			break;
		case 6: // dq0dt
			dxdt[i] = x[i] + 2;
			break;
		case 7: // dq1dt
			dxdt[i] = x[i] + 2;
			break;
		case 8: // dq2dt
			dxdt[i] = x[i] + 2;
			break;
		case 9: // dq3dt
			dxdt[i] = x[i] + 2;
			break;
		}
	}
	return dxdt;
}

std::vector<double> GetSolver(std::vector<double> (*KA)(std::vector<double>(*SU)(std::vector<double> qref, std::vector<double> x), std::vector<double> qref, std::vector<double> x), 
							  std::vector<double> (*SU)(std::vector<double> qref, std::vector<double> x), 
							  std::vector<double> qref, std::vector<double> x, double ht = 1e-4)
{
	/* Численное интегрирование системы ОДУ (RungeKutt ode4) */
	
	std::vector<double> dxdt = KA(SU, qref, x);
	for (int i = 0; i < x.size(); i++)
	{
		x[i] += dxdt[i] * ht;
	}
	return x;
}

int main()
{
	// Целевая ориентация:      { q0, q1, q2, q2 }
	std::vector <double> qref = { 1,  0,  0,  0 };
	// Вектор состояния:     { wx, wy, wz, wxr, wyr, wzr, q0, q1, q2, q3 }
	std::vector <double> x = { 0,  0,  0,  0,   0,	 0,	  1,  0,  0,  0 };


	for (int i = 0; i < x.size(); i++)
	{
		x = GetSolver(KA, SU, qref, x);
	}

}
