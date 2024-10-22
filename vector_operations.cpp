#include "vector_operations.h"

#include <stdexcept>
#include <cmath>

std::vector<double> operator+(const std::vector<double>& vect_1, const std::vector<double>& vect_2) {
    if (vect_1.size() != vect_2.size()) {
        throw std::invalid_argument("–азмеры векторов не совпадают");
    }
    std::vector<double> result_vect;
    result_vect.reserve(vect_1.size());
    for (size_t i = 0; i < vect_1.size(); ++i)
        result_vect.push_back(vect_1[i] + vect_2[i]);
    return result_vect;
}

std::vector<double> operator-(const std::vector<double>& vect_1, const std::vector<double>& vect_2) {
    if (vect_1.size() != vect_2.size()) {
        throw std::invalid_argument("–азмеры векторов не совпадают");
    }
    std::vector<double> result_vect;
    result_vect.reserve(vect_1.size());
    for (size_t i = 0; i < vect_1.size(); ++i)
        result_vect.push_back(vect_1[i] - vect_2[i]);
    return result_vect;
}

std::vector<double>& operator+=(std::vector<double>& vect_1, const std::vector<double>& vect_2) {
    if (vect_1.size() != vect_2.size()) {
        throw std::invalid_argument("–азмеры векторов не совпадают");
    }
    for (size_t i = 0; i < vect_1.size(); ++i)
        vect_1[i] += vect_2[i];
    return vect_1;
}

std::vector<double>& operator-=(std::vector<double>& vect_1, const std::vector<double>& vect_2) {
    if (vect_1.size() != vect_2.size()) {
        throw std::invalid_argument("–азмеры векторов не совпадают");
    }
    for (size_t i = 0; i < vect_1.size(); ++i)
        vect_1[i] -= vect_2[i];
    return vect_1;
}

std::vector<double> operator*(const double& scalar, const std::vector<double>& vect) {
    std::vector<double> result_vect;
    result_vect.reserve(vect.size());
    for (const auto& element : vect)
        result_vect.push_back(element * scalar);
    return result_vect;
}

std::vector<double> operator*(const std::vector<double>& vect, const double& scalar) {
    return scalar * vect;
}

std::vector<double>& operator*=(std::vector<double>& vect, const double& scalar) {
    for (auto& element : vect)
        element *= scalar;
    return vect;
}

std::vector<double> operator/(const std::vector<double>& vect, const double& scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("ƒеление на ноль");
    }
    return (1.0 / scalar) * vect;
}

std::vector<double>& operator/=(std::vector<double>& vect, const double& scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("ƒеление на ноль");
    }
    for (auto& element : vect)
        element /= scalar;
    return vect;
}

double norm_vect(const std::vector<double>& vect) {
    double norm = 0.0;
    for (const auto& element : vect)
        norm += element * element;
    return std::sqrt(norm);
}

std::vector<double> add_vect(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result_vect = a;
    for (size_t i = 0; i < b.size(); i++)
        result_vect.push_back(b[i]);
    return result_vect;
}

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2) {
    if (matr_1.size() != matr_2.size() || matr_1.empty() || matr_1[0].size() != matr_2[0].size()) {
        throw std::invalid_argument("–азмеры матриц не совпадают");
    }
    std::vector<std::vector<double>> result_matr;
    result_matr.reserve(matr_1.size());
    for (size_t i = 0; i < matr_1.size(); ++i) {
        std::vector<double> row;
        row.reserve(matr_1[i].size());
        for (size_t j = 0; j < matr_1[i].size(); ++j)
            row.push_back(matr_1[i][j] + matr_2[i][j]);
        result_matr.push_back(row);
    }
    return result_matr;
}

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2) {
    if (matr_1.size() != matr_2.size() || matr_1.empty() || matr_1[0].size() != matr_2[0].size()) {
        throw std::invalid_argument("–азмеры матриц не совпадают");
    }
    std::vector<std::vector<double>> result_matr;
    result_matr.reserve(matr_1.size());
    for (size_t i = 0; i < matr_1.size(); ++i) {
        std::vector<double> row;
        row.reserve(matr_1[i].size());
        for (size_t j = 0; j < matr_1[i].size(); ++j)
            row.push_back(matr_1[i][j] - matr_2[i][j]);
        result_matr.push_back(row);
    }
    return result_matr;
}

std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr_1, const std::vector<std::vector<double>>& matr_2) {
    if (matr_1.empty() || matr_2.empty() || matr_1[0].size() != matr_2.size()) {
        throw std::invalid_argument("–азмеры матрицы не совпадают при умножении");
    }

    std::vector<std::vector<double>> result_matr(matr_1.size(), std::vector<double>(matr_2[0].size()));

    for (size_t i = 0; i < matr_1.size(); ++i) {
        for (size_t j = 0; j < matr_2[0].size(); ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < matr_2.size(); ++k) {
                sum += matr_1[i][k] * matr_2[k][j];
            }
            result_matr[i][j] = sum;
        }
    }
    return result_matr;
}

std::vector<double> operator*(const std::vector<std::vector<double>>& matr, const std::vector<double>& vect) {
    if (matr.empty() || matr[0].size() != vect.size()) {
        throw std::invalid_argument("–азмер матрицы не соответствует размеру вектора дл€ умножени€");
    }
    std::vector<double> result(matr.size());
    for (size_t i = 0; i < matr.size(); ++i) {
        result[i] = 0.0;
        for (size_t j = 0; j < matr[i].size(); ++j) {
            result[i] += matr[i][j] * vect[j];
        }
    }
    return result;
}

std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matr, const double& scalar) {
    std::vector<std::vector<double>> result(matr.size(), std::vector<double>(matr[0].size()));
    for (size_t i = 0; i < matr.size(); ++i) {
        for (size_t j = 0; j < matr[i].size(); ++j) {
            result[i][j] = matr[i][j] * scalar;
        }
    }
    return result;
}

std::vector<std::vector<double>> operator*(const double& scalar, const std::vector<std::vector<double>>& matr) {
    return matr * scalar;
}

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matr) {
    std::vector<std::vector<double>> result_matr(matr[0].size(), std::vector<double>(matr.size()));
    for (size_t i = 0; i < matr.size(); ++i)
        for (size_t j = 0; j < matr[i].size(); ++j)
            result_matr[j][i] = matr[i][j];
    return result_matr;
}

std::vector<std::vector<double>> getSubmatrix(
    const std::vector<std::vector<double>>& matr,
    const size_t row, const size_t col) {
    if (matr.size() <= 1 || matr[0].size() <= 1 || row >= matr.size() || col >= matr[0].size()) {
        throw std::invalid_argument("Ќедопустимый размер матрицы или индексы");
    }
    std::vector<std::vector<double>> result_matr(matr.size() - 1, std::vector<double>(matr[0].size() - 1));
    size_t r = 0;
    for (size_t i = 0; i < matr.size(); ++i) {
        if (i == row)
            continue;
        size_t c = 0;
        for (size_t j = 0; j < matr[0].size(); ++j) {
            if (j == col)
                continue;
            result_matr[r][c] = matr[i][j];
            ++c;
        }
        ++r;
    }
    return result_matr;
}

double determinant(const std::vector<std::vector<double>>& matr) {
    if (matr.empty() || matr.size() != matr[0].size()) {
        throw std::invalid_argument("ћатрица должна быть квадратной и непустой");
    }
    size_t N = matr.size();
    if (N == 1) {
        return matr[0][0];
    }
    else if (N == 2) {
        return matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0];
    }
    else {
        double det = 0.0;
        for (size_t j = 0; j < N; ++j) {
            double cofactor = (j % 2 == 0) ? 1.0 : -1.0;
            det += cofactor * matr[0][j] * determinant(getSubmatrix(matr, 0, j));
        }
        return det;
    }
}

std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matr) {
    if (matr.empty() || matr.size() != matr[0].size() || determinant(matr) == 0.0) {
        throw std::invalid_argument("Ќедопустимый размер матрицы или определитель равен нулю");
    }
    size_t N = matr.size();
    std::vector<std::vector<double>> inv_matr(N, std::vector<double>(N));
    double inv_det = 1.0 / determinant(matr);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double cofactor = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            inv_matr[j][i] = cofactor * determinant(getSubmatrix(matr, i, j)) * inv_det;
        }
    }
    return inv_matr;
}

std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != 3 || b.size() != 3) {
        throw std::invalid_argument("ƒл€ векторного произведени€ требуютс€ трехмерные векторы");
    }
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("¬екторы должны иметь одинаковый размер дл€ скал€рного произведени€");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

std::vector<std::vector<double>> outerProduct(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<std::vector<double>> result(a.size(), std::vector<double>(b.size()));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

std::vector<double> quatConjugate(const std::vector<double>& a) {
    std::vector<double> result_quat;
    if (a.size() != 4) {
        throw std::invalid_argument(" ватернион должен обладать размерностью - 4");
    }
    result_quat = - a;
    result_quat[0] = - a[0];
    return result_quat;
}

std::vector<double> quatProduct(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result_quat;
    if (a.size() != 4 || b.size() != 4) {
        throw std::invalid_argument("ƒл€ произведени€ кватернионов требуютс€ четырехмерные векторы");
    }
    return {
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]
    };
}

std::vector<double> operator-(const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = -vec[i];
    }
    return result;
}

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& matr) {
    std::vector<std::vector<double>> result(matr.size(), std::vector<double>(matr[0].size()));
    for (size_t i = 0; i < matr.size(); ++i) {
        for (size_t j = 0; j < matr[0].size(); ++j) {
            result[i][j] = -matr[i][j];
        }
    }
    return result;
}
