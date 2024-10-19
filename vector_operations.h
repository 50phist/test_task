#pragma once

#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include <vector>

// ���������� ��������� + ��� ��������
std::vector<double> operator+(const std::vector<double>& vect_1, const std::vector<double>& vect_2);
// ���������� ��������� - ��� ��������
std::vector<double> operator-(const std::vector<double>& vect_1, const std::vector<double>& vect_2);
// ���������� ��������� += ��� ��������
std::vector<double>& operator+=(std::vector<double>& vect_1, const std::vector<double>& vect_2);
// ���������� ��������� -= ��� ��������
std::vector<double>& operator-=(std::vector<double>& vect_1, const std::vector<double>& vect_2);
// ���������� ��������� * ��� ��������� ������� �� ������
std::vector<double> operator*(const double& scalar, const std::vector<double>& vect);
// ���������� ��������� * ��� ��������� ������� �� ������
std::vector<double> operator*(const std::vector<double>& vect, const double& scalar);
// ���������� ��������� *= ��� �������� (��������� �� ������)
std::vector<double>& operator*=(std::vector<double>& vect, const double& scalar);
// ���������� ��������� / ��� ������� ������� �� ������
std::vector<double> operator/(const std::vector<double>& vect, const double& scalar);
// ���������� ��������� /= ��� �������� (������� �� ������)
std::vector<double>& operator/=(std::vector<double>& vect, const double& scalar);
// ����� �������
double norm_vect(const std::vector<double>& vect);
// ���������� ������� b � ����� ������� a
std::vector<double> add_vect(const std::vector<double>& a, const std::vector<double>& b);
// ���������� ��������� + ��� ������
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// ���������� ��������� - ��� ������
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// ���������� ��������� * ��� ������
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// ���������� ��������� * ��� ������������ ������� � �������
std::vector<double> operator*(const std::vector<std::vector<double>>& matr, const std::vector<double>& vect);
// ���������� ��������� * ��� ��������� ������� �� ������
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr, const double& scalar);
// ���������� ��������� * ��� ��������� ������� �� �������
std::vector<std::vector<double>> operator*(const double& scalar, const std::vector<std::vector<double>>& matr);
// ��������������� �������
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matr);
// �������� ���������� ��� ������������ ������ � �������
std::vector<std::vector<double>> getSubmatrix(const std::vector<std::vector<double>>& matr, const size_t row, const size_t col);
// ���������� ������������ �������
double determinant(const std::vector<std::vector<double>>& matr);
// ���������� �������� �������
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matr);
// ��������� ������������ 3-� ������ ��������
std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b);
// ��������� ������������
double dotProduct(const std::vector<double>& a, const std::vector<double>& b);
// ������� ������������ �������� a * transpose(b)
std::vector<std::vector<double>> outerProduct(const std::vector<double>& a, const std::vector<double>& b);
// ���������� �����������
std::vector<double> quatConjugate(const std::vector<double>& a);
// ������������ ������������
std::vector<double> quatProduct(const std::vector<double>& a, const std::vector<double>& b);
// ���������� �������� ��������� -
std::vector<double> operator-(const std::vector<double>& vec);
// ���������� �������� ��������� -
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr);

#endif // !VECTOR_OPERATIONS_H
