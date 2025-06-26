/*!
@file
@brief Краевая задача посадки на Луну
@details Автор: Степан Андреевич Борковой\n
         дата: 26.06.2025
*/
//--------------------------------------------------------------------------------------------
//#define FOR_DEN 999 //!< Если файл компилируется у Д.А. Тучина
//--------------------------------------------------------------------------------------------
#ifdef FOR_DEN
#include "estk.h"
#endif
//--------------------------------------------------------------------------------------------
/*!
@defgroup fun_output функции вывода
@details функции, предназначенные для вывода данных в файлы или в стандартный вывод
//--------------------------------------------------------------------------------------------
@defgroup utils разные функции
@details различные полезные функции
//--------------------------------------------------------------------------------------------
@defgroup ODE_solving решения ОДУ
@details Функции, решающие ОДУ методом Рунге-Кутта 4 порядка
//--------------------------------------------------------------------------------------------
@defgroup convert преобразования кеплеровых и декартовых векторов
@details Функции для перевода между кеплеровыми и декартовыми векторами состояния
*/
//--------------------------------------------------------------------------------------------

#include <math.h>

#include <stdio.h>

#include <stdlib.h>

#define pi 3.1415926
#define mu_moon 49028.01374e8 //!< Гравитационный параметр Луны из модели JGL075D1
#define TRUE 1
#define FALSE 0
//!Требуемое конечное положение аппарата(м, м, м/с, м/с)
const double x_end[4] = {
  0,
  1740140,
  0,
  0
};
const double theta_step = 0.000002; //!< Шаг производной по углу тангажа, рад
const double theta_dot_step = 0.0000005; //!< Шаг производной по изменению угла тангажа, рад/с
const double p_step = 5.0; //!< Шаг производной по тяге, Н
const double t_step = 4.0; //!< Шаг производной по времени торможения, с
//-------------------------------------------------------------------------------------------------------------------
/*!
@brief функция вывода матрицы
@param[in] A матрица размерности N x M
@param[in] N количество строк
@param[in] M количество столбцов
@ingroup fun_output
*/
void print_matrix(double * A, int N, int M) {
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      printf("%.4f ", A[4 * i + j]);
    }
    printf("\n");
  }
}

/*!
@brief функция вывода матрицы в файл
@param[in] A матрица размерности N x M
@param[in] N количество строк
@param[in] M количество столбцов
@param[in] file указатель на файл для вывода
@ingroup fun_output
*/
void fprint_matrix(double * A, int N, int M, FILE * file) {
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      fprintf(file, "%.4f ", A[4 * i + j]);
    }
    fprintf(file, "\n");
  }
}

/*!
@brief функция поиска максимума из 2 чисел
@param[in] a первое число для сравнения
@param[in] b второе число для сравнения
@return max набольшее из двух чисел
@ingroup utils
*/
double mymax(double a, double b) {
  double max;
  max = a;
  if (max < b) max = b;
  return max;
}

/*!
@brief функция поиска минимума из 2 чисел
@param[in] a первое число для сравнения
@param[in] b второе число для сравнения
@return min наименьшее из двух чисел
@ingroup utils
*/
double mymin(double a, double b) {
  double min;
  min = a;
  if (min > b) min = b;
  return min;
}

/*!
@brief функция, задающая уравнение движения
@param[in] x вектор состояния аппарата в декартовых координатах (x, y, v_x, v_y)
@param[in] p тяга аппарата
@param[in] u удельная тяга аппарата
@param[in] theta угол тангажа
@param[in] t время на момент расчёта
@param[out] res вектор для записи состояния аппарата
@param[in] m0 начальная масса аппарата
@ingroup ODE_solving
*/
void f(double * x, double p, double u, double theta, double t, double * res, double m0) {
  double r = pow(x[0], 2) + pow(x[1], 2);
  r = pow(r, 1.5);
  res[0] = x[2];
  res[1] = x[3];
  res[2] = -mu_moon * x[0] / r + u / (m0 * u / p - t) * cos(theta);
  res[3] = -mu_moon * x[1] / r + u / (m0 * u / p - t) * sin(theta);
}

/*!
@brief нахождение модуля разности между 2 векторами размерности 4
@param[in] x первый вектор
@param[in] y второй вектор
@return res результат работы функции
@ingroup utils
*/
double abs_diff(double x[4], double y[4]) {
  double res = 0.0;
  int i = 0;
  for (i = 0; i < 4; i++) {
    res += ((x[i] - y[i]) * (x[i] - y[i]));
    // printf("%f ", res);
  }
  // printf("\n");
  return sqrt(res);
}

/*!
@brief нахождение модуля вектора
@param[in] n размерность вектора
@param[in] x вектор
@return res результат работы функции
@ingroup utils
*/
double abs_vec(int n, double * x) {
  double res = 0.0;
  int i = 0;
  for (i = 0; i < n; i++) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}

//-----------------------------------------------------------------------------------------------------
#define MAXINVMATR 15 //!< Максимальная размерность матрицы для обращения
/*!
@brief  Обращение матрицы методом Гаусса
@param [in] A квадратная матрица A, записанная по строкам
@param [out] inv_A обратная матрица к матрице А
@param [in] dimA количество строк или столбцов квадратной матрицы
@return Функция возвращает
@return 0 -- успешное обращение матрицы
        1 -- размерность матрицы больше MAXINVMATR
        2 -- определитель матрицы меньше машинного эпсилона (1е-17)
@details Алгоритм из Б. П. Демидович, И. А. Марон. Основы вычислительной математики, 1960, 272с.
*/
int InvMatrixGauss(double * A, double * inv_A, int dimA) {
  double a[MAXINVMATR * MAXINVMATR], b[MAXINVMATR], x[MAXINVMATR], r3, s;
  int i, j, m, n, i0, n1, l, k, j0, m1, m2, m3;
  if (dimA > MAXINVMATR) return 1;
  m = dimA * dimA;
  n = dimA;
  for (i = 0; i < dimA; i += 1) {
    for (j = 0; j < m; j += 1)
      a[j] = A[j];
    for (j = 0; j < dimA; j += 1)
      b[j] = 0.;
    b[i] = 1.;
    n1 = n - 1;
    for (k = 0; k < n1; k += 1) {
      m1 = k * n;
      for (i0 = k + 1, l = k; i0 < n; i0 += 1) {
        if (fabs(a[i0 * n + k]) > fabs(a[l * n + k]))
          l = i0;
      }
      if (l != k) {
        m2 = l * n;
        for (j0 = k; j0 < n; j0 += 1) {
          r3 = a[m1 + j0];
          a[m1 + j0] = a[m2 + j0];
          a[m2 + j0] = r3;
        }
        r3 = b[k];
        b[k] = b[l];
        b[l] = r3;
      }
      m2 = m1 + k;
      /*for(int i = 0; i < 4; i++){
          for(int j = 0; j<4; j++){
              printf("%.10f ", a[4*i + j]);
          }
          printf("\n");
      }*/
      if (fabs(a[m2]) < 1.e-17) {
        return 2; //определитель матрицы меньше машэпсилона
      }
      for (i0 = k + 1; i0 < n; i0 += 1) {
        m3 = i0 * n;
        r3 = a[m3 + k] / a[m2];
        a[m3 + k] = 0.;
        for (j0 = k + 1; j0 < n; j0 += 1)
          a[m3 + j0] -= (r3 * a[m1 + j0]);
        b[i0] -= (r3 * b[k]);
      }
    }
    x[n1] = b[n1] / a[n1 * n + n1];
    for (i0 = n1 - 1; i0 >= 0; i0 -= 1) {
      s = 0.;
      for (j0 = i0 + 1; j0 < n; j0 += 1)
        s -= a[i0 * n + j0] * x[j0];
      x[i0] = (b[i0] + s) / a[i0 * n + i0];
    }
    for (j = 0; j < dimA; j += 1)
      inv_A[j * dimA + i] = x[j];
  }
  return 0;
}
//-----------------------------------------------------------------------------------------------------
/*!
@brief функция, расчитывающая один шаг методом Рунге-Кутта в дифференциальном уравнении
@param[out] x_res вектор для записи решения
@param[in] N размерность вектора состояния аппарата
@param[in] x вектор состояния аппарата в декартовых координатах (x, y, v_x, v_y) (м, м, м/с, м/с)
@param[in] t время на момент расчёта, с
@param[in] h шаг метода, с
@param[in] p тяга аппарата, Н
@param[in] u удельная тяга аппарата, м/с
@param[in] theta0 начальный угол тангажа, рад
@param[in] theta_dot скорость изменения угла тангажа, рад/с
@param[in] m0 начальная масса аппарата, кг
@param[in] memorize записывать или не записывать результат в файл
@param[in] file указатель на файл для вывода
@ingroup ODE_solving
*/
void runge_kutta_step(double * x_res, int N, double * x, double t, double h, double p, double u,
  double theta0, double theta_dot, double m0, int memorize, FILE * file) {
  int i;
  double * x_tmp = malloc(N * sizeof(double));
  double theta = theta0 + theta_dot * t;
  double k1[4], k2[4], k3[4], k4[4];
  f(x, p, u, theta, t, k1, m0);
  for (i = 0; i < N; i++) {
    x_tmp[i] = x[i] + k1[i] * h / 2;
  }
  f(x_tmp, p, u, theta0 + theta_dot * (t + h / 2), t + h / 2, k2, m0);
  for (i = 0; i < N; i++) {
    x_tmp[i] = x[i] + k2[i] * h / 2;
  }
  f(x_tmp, p, u, theta0 + theta_dot * (t + h / 2), t + h / 2, k3, m0);
  for (i = 0; i < N; i++) {
    x_tmp[i] = x[i] + k3[i] * h;
  }
  f(x_tmp, p, u, theta0 + theta_dot * (t + h), t + h, k4, m0);
  for (i = 0; i < N; i++) {
    x_res[i] = x[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
  }
  if (memorize) {
    fprintf(file, "%f %f\n", x_res[0], x_res[1]);
  }
  free(x_tmp);
}
//-----------------------------------------------------------------------------------------------------
/*!
@brief функция, расчитывающая полную траекторию аппарата методом Рунге-Кутта в дифференциальном уравнении
@param[in] x_prev вектор для записи решения
@param[in] N размерность вектора состояния аппарата
@param[in] x вектор состояния аппарата в декартовых координатах (x, y, v_x, v_y), (м, м, м/с, м/с)
@param[in] t_start время начала расчёта, с
@param[in] t_end время конца расчёта, с
@param[in] step шаг метода, с
@param[in] p тяга аппарата, Н
@param[in] u удельная тяга аппарата, м/с
@param[in] theta0 начальный угол тангажа, рад
@param[in] theta_dot скорость изменения угла тангажа, рад/с
@param[in] m0 начальная масса аппарата, кг
@param[in] memorize записывать или не записывать результат в файл
@param[in] file указатель на файл для вывода
@ingroup ODE_solving
*/
void runge_kutta_total(double * x_prev, int N, double * x, double t_start, double t_end, double step, double p, double u,
  double theta0, double theta_dot, double m0, int memorize, FILE * file) {
  double t = t_start;
  double x_next[4];
  int i;
  for (i = 0; i < N; i++) {
    x_prev[i] = x[i];
  }
  while ((t + step <= t_end)) { // с условием на невозможность перелёта по координате Х проблема -- производные должны иногда перелетать, чтобы корректно считались.
    runge_kutta_step(x_next, N, x_prev, t, step, p, u, theta0, theta_dot, m0, memorize, file);
    for (i = 0; i < N; i++) x_prev[i] = x_next[i];
    t += step;
  }

}

/*!
@brief нахождение времени изменения скорости по формуле Циолковского
@param[in] vhar абсолютное значение разности скоростей, м/с
@param[in] draft тяга аппарата, Н
@param[in] m масса аппарата, кг
@param[in] m_dot скорость изменения массы аппарата, кг/сек
@return время, необходимое для изменения скорости, сек
@ingroup utils
*/
double TsiolkovskyFormula(double vhar, double draft, double m, double m_dot) {
  double res = m / m_dot * (1. - exp(-vhar * m_dot / (draft)));
  return (res);
}

/*!
@brief функция для перемножения двух матриц
@param[in] A первая матрица размерности m1xn1
@param[in] m1 количество строк в матрице А
@param[in] n1 количество столбцов в матрице А
@param[in] B вторая матрица размерности m2xn2
@param[in] m2 количество строк в матрице В
@param[in] n2 количество столбцов в матрице B
@return res матрица размерности m1xn2, полученная в результате умножения A*B 
@ingroup utils
*/
double * matrix_multiply(double * A, int m1, int n1, double * B, int m2, int n2) {
  int i, j, k;
  double * res;
  if (n1 != m2) return NULL;
  //int dimA, dimB;
  //dimA = m1*n1;
  //dimB = m2*n2;
  res = malloc(m1 * n2 * sizeof(double));
  for (i = 0; i < m1; i++)
    for (j = 0; j < n2; j++) {
      res[i * n2 + j] = 0;
      for (k = 0; k < n1; k++)
        res[i * n2 + j] += A[i * n1 + k] * B[k * n2 + j];
    }
  return res;
}

/*!
@brief функция знака числа с плавающей точкой
@param[in] a число
@return -1 -- число меньше отрицательного машинного эпсилона (1е-17)
        1 -- число больше отрицательного машинного эпсилона (1е-17)
@ingroup utils
*/
int sign(double a) {
  if (a > -1e-17) return 1;
  else return -1;
}

/*!
@brief функция приведения произвольного числа к отрезку [-2*pi, 2*pi]
@param[in, out] M число
@ingroup utils
*/
double shrink_to_2pi(double M) {
  while ((M > 2 * pi) || (M < 0)) {
    M -= sign(M) * 2 * pi;
  }
  return M;
}
//
/*!
@brief функция, находящая скалярное произведение двумерных векторов
@param[in] x_1 координата х первого вектора
@param[in] y_1 координата y первого вектора
@param[in] x_2 координата х второго вектора
@param[in] y_2 координата y второго вектора
@return результат скалярного произведения
@ingroup utils
*/
double dot_product(double x_1, double y_1, double x_2, double y_2) {
  return x_1 * x_2 + y_1 * y_2;
}
/*!
@brief функция, находящая векторное произведение двумерных векторов
@param[in] x_1 координата х первого вектора
@param[in] y_1 координата y первого вектора
@param[in] x_2 координата х второго вектора
@param[in] y_2 координата y второго вектора
@return результат векторного произведения
@ingroup utils
*/
double cross_product(double x_1, double y_1, double x_2, double y_2) {
  return x_1 * y_2 - y_1 * x_2;
}
/*!
@brief функция, находящая Е из уравнения Кеплера
@param[in] e эксцентриситет орбиты
@param[in] E эксцентрическая аномалия, рад
@param[in] M средняя аномалия, рад
@return результат вычисления
@ingroup utils
*/
double calc_E(double e, double E, double M) {
  return E - e * sin(E) - M;
}
/*!
@brief функция, находящая производную по E уравнения Кеплера
@param[in] e эксцентриситет орбиты
@param[in] E эксцентрическая аномалия, рад
@return результат вычисления
@ingroup utils
*/
double calc_dE(double e, double E) {
  return 1 - e * cos(E);
}
/*!
@brief функция, находящая эксцентрическую аномалию из уравнения Кеплера методом Ньютона
@param[in] e эксцентриситет орбиты
@param[in] E0 начальное приближение, рад
@param[in] M средняя аномалия, рад
@param[in] eps требуемая точность
@return результат вычисления
@ingroup utils
*/
double Newton_Raphson(double e, double E0, double M, double eps) {
  double E = E0;
  while (fabs(calc_E(e, E, M) / calc_dE(e, E)) > eps) {
    E = E - calc_E(e, E, M) / calc_dE(e, E);
  }
  return E;
}

/*!
@brief перевод кеплеровых элементов орбиты в декартовый вектор состояния в двумерном случае
в нашем случае -- селеноцентрическая система координат, ось Ох направлена из центра Луны в направлении перицентра, ось Оу -- по направлению
движения аппарата в перицентре
алгоритм взят из
    M.Eng.ReneSchwarz(rene-schwarz.com):MemorandumSeries
    KeplerianOrbitElements -> CartesianStateVectors(Memorandum №1)
@param[in] N размерность требуемого вектора (4)
@param[in] K вектор кеплеровых элементов орбиты, состоящйи из:\n
    0: Большая полуось a [м]\n
    1: Эксцентриситет e [1]\n
    2: Аргумент перицентра (0 по построению СК) [рад]\n
    3: Долгота восходящего узла (не нужна в 2D) (LAN)[рад]\n
    4: Наклонение (не нужно в 2D) i [рад]\n
    5: Средняя аномалия M0 = M(t0) [rad] в момент времени t0 [JD]\n
@param[in] t0 начальный момент времени, с
@param[in] t момент времени t, если отличается от t0
@param[out] res возвращаемый четырехмерный вектор состояния в декартовой системе координат (м, м, м/с, м/с)
@ingroup convert
*/
void keplerian_2_cartesian_2D(int N, double * K, double t0, double t, double * res) {
  double delta_t, M, E, true_anomaly, r, tmp;
  delta_t = t - t0;
  M = K[5] + delta_t * sqrt(mu_moon / pow(K[0], 3));
  M = shrink_to_2pi(M);
  E = Newton_Raphson(K[1], M, M, 1e-5);
  true_anomaly = 2 * atan2(sqrt(1 + K[1]) * sin(E / 2), sqrt(1 - K[1]) * cos(E / 2));
  r = K[0] * (1 - K[1] * cos(E));
  res[1] = r * cos(true_anomaly);
  res[0] = -r * sin(true_anomaly);
  tmp = sqrt(mu_moon * K[0]) / r;
  res[3] = tmp * (-sin(E));
  res[2] = -tmp * (cos(E) * sqrt(1 - K[1] * K[1]));
  return;
}
/*!
@briefперевод элементов декартовой системы координат в кеплеровы элементы орбиты
в нашем случае -- селеноцентрическая система координат, ось Ох направлена из центра Луны в направлении перицентра, ось Оу -- по направлению
движения аппарата в перицентре
алгоритм взят из
    M.Eng.ReneSchwarz(rene-schwarz.com):MemorandumSeries
    KeplerianOrbitElements -> CartesianStateVectors(Memorandum №1)
@param[in] N размерность требуемого вектора (4)
@param[in] x вектор состояния аппарата (м, м, м/с, м/с)
@param[in] t0 начальный момент времени, с
@param[in] t момент времени t, если отличается от t0
@param[out] res возвращаемый вектор кеплеровых элементов положения аппарата на орбите
@ingroup convert
*/
void cartesian_2_keplerian_2D(int N, double * x, double t0, double t, double * res) {
  double /*delta_t,*/ r, v, r_dot_v, e_vec_x, e_vec_y, nu, E;
  //    delta_t = t - t0;
  res[2] = 0;
  res[3] = 0;
  res[4] = 0;
  r = sqrt(x[0] * x[0] + x[1] * x[1]);
  v = sqrt(x[2] * x[2] + x[3] * x[3]);
  res[0] = 1 / (2 / r - v * v / mu_moon);
  r_dot_v = dot_product(x[0], x[1], x[2], x[3]);
  e_vec_x = (v * v - mu_moon / r) * x[0] - r_dot_v * x[2];
  e_vec_y = (v * v - mu_moon / r) * x[1] - r_dot_v * x[3];
  res[1] = sqrt(e_vec_x * e_vec_x + e_vec_y * e_vec_y) / mu_moon;

  nu = atan2(x[1], x[0]); // Истинная аномалия
  E = atan2(sqrt(1 - res[1]) * sin(nu / 2), sqrt(1 + res[1]) * cos(nu / 2)) * 2; // Эксцентрическая аномалия

  res[5] = E - res[1] * sin(E); // Ур-е Кеплера
  return;
}
#ifdef FOR_DEN
int main_Stepan(void) {
    #else
    /*!
    @brief основная функция программы
    @details В функции представлена инициализация параметров при помощи перевода кеплеровых и декартовых элементов орбиты,\n
    нахождение начального приближения для решения уравнения методом Ньютона и его реализация. 
    */
    int main(void) {
      #endif

      double theta = 0; //начальный угол тангажа, радиан
      //	double theta_dot = 0;		//скорость изменения угла тангажа, радиан/с
      const int N = 4; //размерность неизвестного фазового вектора
      const double u = 3126.0; //удельная тяга, м/с
      //четырехмерный кинематический вектор начального положения и скорости орбитального движения (м, м/с)
      double x0[4], K[6];
      int p_iter;
      int i, j, iter, flag;
      double t_end;
      FILE * file;
      double t;
      double t_plus;
      double params[4];
      double x_res[4];
      double m0;
      //вектор конечного положения аппарата в Кеплеровых элементах орбиты
      for (i = 0; i < 6; i++) {
        K[i] = 0;
      }
      K[0] = 1796000;
      K[1] = 0.02283;
      t_end = 0;
      keplerian_2_cartesian_2D(6, K, 0, -t_end, x0);
      //время набора характеристической скорости по формуле Циолковского
      t = 0.;
      t_plus = 0.1; //0.2;
      flag = 0;
      theta = atan2(fabs(x0[3]), fabs(x0[2]));
      iter = 0;
      file = fopen("flight.txt", "w");
      if (file == NULL) {
        printf("Ошибка при открытии файла!\n");
      }
      fprint_matrix(params, 1, 4, file);
      m0 = 1250;
      t_end = 360;
      for (p_iter = 0; p_iter < 1; p_iter++) {
        double residual;
        printf("Counting parameters for T = %f\n", t_end);
        //printf("%f\n", t_end);
        //положение в момент времени t_end
        keplerian_2_cartesian_2D(6, K, 0, -t_end, x0);
        //printf("initial position vector: ");
        //print_matrix(x0, 1, 4);
        t = 0.;
        flag = 0;
        theta = atan2(fabs(x0[3]), fabs(x0[2]));
        params[0] = theta;
        params[1] = 0;
        params[2] = 3000;
        params[3] = TsiolkovskyFormula(sqrt(x0[2] * x0[2] + x0[3] * x0[3]), params[2], m0, params[2] / u);
        params[1] = -theta / params[3];
        //params -- параметры спутника для посадки, обозначены как вектор х в pdf
        printf("initial parameters vector: ");
        print_matrix(params, 1, 4);
        printf("initial position: ");
        print_matrix(x0, 1, 4);
        iter = 0;
        residual = 50000;
        fprintf(file, "0");
        runge_kutta_total(x_res, N, x0, t, t_end, t_plus, 0, u, params[0], params[1], m0, TRUE, file);
        while ((iter < 100) && (residual > 1)) {
          //пролетаем до момента времени t_end методом РК
          double diff[4];
          double A[4 * 4];
          double x_theta_plus[4], x_theta_minus[4], x_theta_dot_plus[4], x_theta_dot_minus[4], x_p_plus[4], x_p_minus[4], x_t_plus[4], x_t_minus[4];
          double corr[4];
          double * A_inv = malloc(16 * sizeof(double));
          double x0_plus[4];
          double x0_minus[4];\
          iter++;
          t = 0.;
          fprintf(file, "%d\n", iter);
          runge_kutta_total(x_res, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, TRUE, file);
          for (i = 0; i < N; i++) {
            diff[i] = 0.0;
          }
          for (i = 0; i < N; i++) {
            diff[i] = (x_res[i] - x_end[i]);
          }
          // Шаги для численного дифференцирования
          // Вычисление производных по theta0
          params[0] += theta_step;
          runge_kutta_total(x_theta_plus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[0] -= 2 * theta_step;
          runge_kutta_total(x_theta_minus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[0] += theta_step;
          for (i = 0; i < N; i++) {
            A[i * 4 + 0] = (x_theta_plus[i] - x_theta_minus[i]) / (2 * theta_step);
          }
          // Вычисление производных по theta_dot
          params[1] += theta_dot_step;
          runge_kutta_total(x_theta_dot_plus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[1] -= 2 * theta_dot_step;
          runge_kutta_total(x_theta_dot_minus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[1] += theta_dot_step;

          for (i = 0; i < N; i++) {
            A[i * 4 + 1] = (x_theta_dot_plus[i] - x_theta_dot_minus[i]) / (2 * theta_dot_step);
          }

          // Вычисление производных по p
          params[2] += p_step;
          runge_kutta_total(x_p_plus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[2] -= 2 * p_step;
          runge_kutta_total(x_p_minus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[2] += p_step;

          for (i = 0; i < N; i++) {
            A[i * 4 + 2] = (x_p_plus[i] - x_p_minus[i]) / (2 * p_step);
          }

          // Вычисление производных по t_k
          params[3] += t_step;
          runge_kutta_total(x_t_plus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[3] -= 2 * t_step;
          runge_kutta_total(x_t_minus, N, x0, t, params[3], t_plus, params[2], u, params[0], params[1], m0, FALSE, file);
          params[3] += t_step;

          for (i = 0; i < N; i++) {
            A[i * 4 + 3] = (x_t_plus[i] - x_t_minus[i]) / (2 * t_step);
          }
          for (i = 0; i < 16; i++) {
            A_inv[i] = 0;
          }
          flag = InvMatrixGauss(A, A_inv, 4);
          if (flag != 0) {
            printf("Matrix inversion error!\n");
            break;
          }
          for (i = 0; i < 4; i++) {
            corr[i] = 0;
            for (j = 0; j < N; j++) {
              corr[i] -= A_inv[i * 4 + j] * diff[j];
            }
          }
          for (i = 0; i < 4; i++) {
            params[i] += corr[i];
          }

          free(A_inv);
          //  printf("%d\n", flag);
          //  print_matrix(A, 4, 4);
          //	printf("initial position: ");
          //	print_matrix(x0, 1, 4);
          //	printf("time %f\n", t);
          printf("final position with given parameters: ");
          print_matrix(x_res, 1, 4);
          //	printf("difference in vectors: ");
          //	print_matrix(diff, 1, 4);
          //	printf("absolute difference: %f\n", abs_diff(x_res, x_end));
          //	printf("Flag = %d \n", flag);
          //	printf("parameters vector after correction: ");
          //	fprint_matrix(params, 1, 4, file);
          //	printf("new initial posiiton: ");
          //	print_matrix(x0, 1, 4);
          residual = sqrt((diff[0] * diff[0] + diff[1] * diff[1]) + 10 * (diff[2] * diff[2] + diff[3] * diff[3]));
          printf("residual functional = %f\n", residual);
          //	printf("%f %f\n", sqrt(diff[0]*diff[0]+diff[1]*diff[1]), sqrt(diff[2]*diff[2]+diff[3]*diff[3]));

          //printf("%d iteration end\n\n", iter);
        }
        t_end -= 5;
        printf("end of an algorithm after %d iterations\n", iter);
        printf("final position vector: ");
        print_matrix(x_res, 1, 4);
        printf("final parameters vector: ");
        printf("%f %f %f %f\n\n", params[0] * 180 / pi, params[1] * 180 / pi, params[2], params[3]);
      }

      fclose(file);
      return 0;
    }