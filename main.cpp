#include "vector_operations.h"
#include<iostream>
#include<vector>

std::vector<double> SU(std::vector<double> qref, std::vector<double> x)
{
	/* ������-������� ������ ����� ������� ��� ������ �� */

	double Jr = 5e-4, Hmax = 0.2, Mmax = 5e-3; 
	double ki = 100, kp = 10;

	std::vector <double> w, H, q, qerr;
	std::vector<double> dHdt;

	w = { x[0], x[1], x[2] };
	H = { x[3], x[4], x[5] };
	q = { x[6], x[7], x[8], x[9] };
	H *= Jr;
	
	// qerr = qref - q;
	qerr = quatProduct(qref, quatConjugate(q));
	qerr.erase(qerr.begin());
	dHdt = ki * qerr + kp * w;

	for (size_t i = 0; i < dHdt.size(); i++)
	{
		if ((abs(H[i]) <= Hmax) && (abs(dHdt[i]) <= Mmax))
		{
			dHdt[i] = 0;
		}
	}
	
	return dHdt;
}

std::vector<double> KA(std::vector<double>(*SU)(std::vector<double> qref, std::vector<double> x), 
					   std::vector<double> qref, std::vector<double> x)
{
	/* ������-������� ������ ����� ������� ��� �������� �� */

	std::vector < std::vector <double> > J(3, std::vector <double>(3)), J_inv(3, std::vector <double>(3));
	J = { { 3.0, -0.1, 0.0 }, { -0.1, 2.3, 0.2 }, { 0.0, 0.2, 1.9 } };
	J_inv = { { 0.333822, 0.0146481, -0.0015419 }, { 0.0146481, 0.439442, -0.046257 }, { -0.0015419, -0.046257, 0.531185 } };
	double Jr = 5e-4;

	std::vector <double> w, qw, H, q;
	std::vector<double> dxdt, dwdt, dwrdt, dqdt;

	w = { x[0], x[1], x[2] };
	qw = { 0, x[0], x[1], x[2] };
	H = { x[3], x[4], x[5] };
	q = { x[6], x[7], x[8], x[9] };
	H *= Jr;

	dxdt = SU(qref, x); // dHdt
	dwdt = -1 * J_inv * ( dxdt + crossProduct(w, J * w) + crossProduct(w, H) );
	dwrdt = dxdt / Jr;
	dqdt = 0.5 * quatProduct(q, qw);

	dxdt = dwdt;
	dxdt = add_vect(dxdt, dwrdt);
	dxdt = add_vect(dxdt, dqdt);

	return dxdt;
}

std::vector<double> GetSolver(std::vector<double> (*KA)(std::vector<double>(*SU)(std::vector<double> qref, std::vector<double> x), std::vector<double> qref, std::vector<double> x), 
							  std::vector<double> (*SU)(std::vector<double> qref, std::vector<double> x), 
							  std::vector<double> qref, std::vector<double> x, const double ht)
{
	/* ��������� �������������� ������� ��� (RungeKutt ode4) */

	std::vector<double> xi = x, dxdti, k1(x.size()), k2(x.size()), k3(x.size()), k4(x.size()), q;

	dxdti = KA(SU, qref, xi);
	for (size_t i = 0; i < dxdti.size(); i++)
	{
		k1[i] = dxdti[i] * ht;
		xi[i] = x[i] + k1[i] / 2;
	}
	dxdti = KA(SU, qref, xi);
	for (size_t i = 0; i < dxdti.size(); i++)
	{
		k2[i] = dxdti[i] * ht;
		xi[i] = x[i] + k2[i] / 2;
	}
	dxdti = KA(SU, qref, xi);
	for (size_t i = 0; i < dxdti.size(); i++)
	{
		k3[i] = dxdti[i] * ht;
		xi[i] = x[i] + k3[i];
	}
	dxdti = KA(SU, qref, xi);
	for (size_t i = 0; i < dxdti.size(); i++)
	{
		k4[i] = dxdti[i] * ht;
	}
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
	}

	q = { x[6], x[7], x[8], x[9] };
	q /= norm_vect(q);
	for (size_t i = 0; i < q.size(); i++)
	{
		x[6 + i] = q[i];
	}
	return x;
}

int main()
{
	// ������� ����������:      { q0, q1, q2, q2 }
	std::vector <double> qref = { 0.85355,  -0.14644,  0.35355,  0.35355 }; // {Rz = 45, Ry = 45, Rz = 0}
	// ������ ���������:     { wx, wy, wz, wxr, wyr, wzr, q0, q1, q2, q3 }
	std::vector <double> x = { 0,  0,  0,  0,   0,	 0,	  1,  0,  0,  0 }, dHdt;
	double t = 0, ht = 1e-3, hrec = 1e-1;
	int rec = (int)(hrec / ht);

	while (t < 20)
	{
		dHdt = SU(qref, x);
		if (rec == (int)(hrec / ht))
		{
			std::cout << "t =\t" << t << ";\tdHdt =\t" << dHdt[0] << "\t" << dHdt[1] << "\t" << dHdt[2] << ";\tq =\t" << x[6] << "\t" << x[7] << "\t" << x[8] << "\t" << x[9] << ";" << std::endl;
			std::cout << std::endl;
			rec = 0;
		}
		x = GetSolver(KA, SU, qref, x, ht);
		t += ht;
		rec++;
	}

}
