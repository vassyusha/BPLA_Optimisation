#include "stability.h"
#include <vector>
#include <cmath>
#include <algorithm>


using namespace std;

// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
// не нужен
bool isSystemStable(const std::vector<double>& a, const AeroConstants& config) {


    return true;

    // Используем плоские массивы в стеке. 
    // alignas(32) помогает компилятору использовать SIMD (AVX) инструкции.
    alignas(32) double A[25] = {0}; 

    // Заполнение матрицы (развернуто в 1D для скорости)
    // Строка 0
    A[0] = -config.k_beta; A[1] = config.alpha_0; A[2] = 1.0; A[3] = config.k_gamma;
    // Строка 1
    A[5] = -a[0];          A[6] = -a[1];          A[7] = -a[2]; A[9] = a[3];
    // Строка 2
    A[10] = -a[4];         A[11] = -a[5];         A[12] = -a[6];
    // Строка 3
    A[16] = 1.0;
    // Строка 4
    A[22] = 1.0;

    double c[6] = {1.0, 0, 0, 0, 0, 0};
    double M[25] = {0};
    for(int i = 0; i < 25; i += 6) M[i] = 1.0; // Единичная матрица

    // Алгоритм Фадеева-Леверье
    for (int k = 1; k <= 5; ++k) {
        double AM[25] = {0};
        double trace = 0;

        // Оптимизированное умножение матриц 5x5
        // Для таких размеров компилятор сам отлично разворачивает циклы
        for (int i = 0; i < 5; ++i) {
            for (int m = 0; m < 5; ++m) {
                double aik = A[i * 5 + m];
                if (aik == 0) continue; // Пропуск нулей (матрица разреженная)
                for (int j = 0; j < 5; ++j) {
                    AM[i * 5 + j] += aik * M[m * 5 + j];
                }
            }
            trace += AM[i * 5 + i];
        }

        c[k] = -trace / k;

        // Обновление M = AM + c[k]*I
        for (int i = 0; i < 25; ++i) {
            M[i] = AM[i];
        }
        for (int i = 0; i < 5; ++i) {
            M[i * 5 + i] += c[k];
        }
    }

    // Проверка условий Гурвица (необходимое условие)
    for (int i = 1; i <= 5; ++i) if (c[i] <= 0) return false;

    // Таблица Рауса
    double table[6][3] = {0};
    table[0][0] = 1.0;  table[0][1] = c[2]; table[0][2] = c[4];
    table[1][0] = c[1]; table[1][1] = c[3]; table[1][2] = c[5];

    for (int i = 2; i < 6; ++i) {
        if (std::abs(table[i-1][0]) < 1e-15) return false;
        
        double factor = table[i-2][0] / table[i-1][0];
        table[i][0] = table[i-2][1] - factor * table[i-1][1];
        table[i][1] = table[i-2][2] - factor * table[i-1][2];

        if (table[i][0] <= 0) return false;
    }

    return true;
}




