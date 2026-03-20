#include "sofama.h"
#include <cmath>
#include <random>

std::vector<double> runSofama(const SofamaParams& p, std::function<double(const std::vector<double>&)> targetFunc) {
    // 1-4. Инициализация популяции и индикаторов mu
    // 5. Выбор лучшей точки z*
    // 6-7. Выбор опорной и доп. точек (Pj)
    // 8-14. Анизотропная мутация (сильная/слабая) с распределением Коши
    // 15-16. Селекция (обновление mu)
    // 17. Вывод оптимальной точки
    return {0,1,2,3,4,5,6}; 
}