#pragma once

#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include <vector>

// Перегрузка оператора + для векторов
std::vector<double> operator+(const std::vector<double>& vect_1, const std::vector<double>& vect_2);
// Перегрузка оператора - для векторов
std::vector<double> operator-(const std::vector<double>& vect_1, const std::vector<double>& vect_2);
// Перегрузка оператора += для векторов
std::vector<double>& operator+=(std::vector<double>& vect_1, const std::vector<double>& vect_2);
// Перегрузка оператора -= для векторов
std::vector<double>& operator-=(std::vector<double>& vect_1, const std::vector<double>& vect_2);
// Перегрузка оператора * для умножения вектора на скаляр
std::vector<double> operator*(const double& scalar, const std::vector<double>& vect);
// Перегрузка оператора * для умножения вектора на скаляр
std::vector<double> operator*(const std::vector<double>& vect, const double& scalar);
// Перегрузка оператора *= для векторов (умножение на скаляр)
std::vector<double>& operator*=(std::vector<double>& vect, const double& scalar);
// Перегрузка оператора / для деления вектора на скаляр
std::vector<double> operator/(const std::vector<double>& vect, const double& scalar);
// Перегрузка оператора /= для векторов (деление на скаляр)
std::vector<double>& operator/=(std::vector<double>& vect, const double& scalar);
// Норма вектора
double norm_vect(const std::vector<double>& vect);
// Добавление вектора b в конец вектора a
std::vector<double> add_vect(const std::vector<double>& a, const std::vector<double>& b);
// Перегрузка оператора + для матриц
std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// Перегрузка оператора - для матриц
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// Перегрузка оператора * для матриц
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2);
// Перегрузка оператора * для перемножения матрицы и вектора
std::vector<double> operator*(const std::vector<std::vector<double>>& matr, const std::vector<double>& vect);
// Перегрузка оператора * для умножения матрицы на скаляр
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr, const double& scalar);
// Перегрузка оператора * для умножения скаляра на матрицу
std::vector<std::vector<double>> operator*(const double& scalar, const std::vector<std::vector<double>>& matr);
// Транспонировать матрицу
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matr);
// Получить подматрицу без определенной строки и столбца
std::vector<std::vector<double>> getSubmatrix(const std::vector<std::vector<double>>& matr, const size_t row, const size_t col);
// Вычисление определителя матрицы
double determinant(const std::vector<std::vector<double>>& matr);
// Вычисление обратной матрицы
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matr);
// Векторное произведение 3-х мерных векторов
std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b);
// Скалярное произведение
double dotProduct(const std::vector<double>& a, const std::vector<double>& b);
// Внешнее произведение векторов a * transpose(b)
std::vector<std::vector<double>> outerProduct(const std::vector<double>& a, const std::vector<double>& b);
// Сопряжение кватерниона
std::vector<double> quatConjugate(const std::vector<double>& a);
// Произведение кватернионов
std::vector<double> quatProduct(const std::vector<double>& a, const std::vector<double>& b);
// Перегрузка унарного оператора -
std::vector<double> operator-(const std::vector<double>& vec);
// Перегрузка унарного оператора -
std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr);

#endif // !VECTOR_OPERATIONS_H
