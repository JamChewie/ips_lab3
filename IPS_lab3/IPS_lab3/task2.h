#ifndef TASK2_H
#define TASK2_H
#include <iostream>

namespace TASK2 
{
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
	// количество строк в исходной квадратной матрице
	const int MATRIX_SIZE = 1500;

	/// Функция InitMatrix() заполняет переданную в качестве 
	/// параметра квадратную матрицу случайными значениями
	/// matrix - исходная матрица СЛАУ
	void InitMatrix(double** matrix)
	{
		for (int i = 0; i < MATRIX_SIZE; ++i)
		{
			matrix[i] = new double[MATRIX_SIZE + 1];
		}

		for (int i = 0; i < MATRIX_SIZE; ++i)
		{
			for (int j = 0; j <= MATRIX_SIZE; ++j)
			{
				matrix[i][j] = rand() % 2500 + 1;
			}
		}
	}

	/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
	/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
	/// последний столбей матрицы - значения правых частей уравнений
	/// rows - количество строк в исходной матрице
	/// result - массив ответов СЛАУ
	double SerialGaussMethod(double **matrix, const int rows, double* result)
	{
		int k;
		double koef;
		// прямой ход метода Гаусса
		auto start = clock() / 1000.0;	
		for (k = 0; k < rows; ++k)
		{
			//
			for (int i = k + 1; i < rows; ++i)
			{
				koef = -matrix[i][k] / matrix[k][k];

				for (int j = k; j <= rows; ++j)
				{
					matrix[i][j] += koef * matrix[k][j];
				}
			}
		}
		auto end = clock() / 1000.0;
		auto duration = end - start;
		// обратный ход метода Гаусса
		result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
		for (k = rows - 2; k >= 0; --k)
		{
			result[k] = matrix[k][rows];
			
			for (int j = k + 1; j < rows; ++j)
			{
				result[k] -= matrix[k][j] * result[j];
			}

			result[k] /= matrix[k][k];
		}
		return duration;
	}


	void run_task2()
	{
		srand((unsigned)time(0));

		double **test_matrix = new double*[MATRIX_SIZE];

		// массив решений СЛАУ
		double *result = new double[MATRIX_SIZE];

		// инициализация матрицы
		InitMatrix(test_matrix);

		auto duration = SerialGaussMethod(test_matrix, MATRIX_SIZE, result);

		for (auto i = 0; i < MATRIX_SIZE; ++i)
		{
			delete[]test_matrix[i];
		}

		std::cout << "Solution" << std::endl;
		for (auto i = 0; i < MATRIX_SIZE; ++i)
		{
			std::cout << "x(" << i << ") = " << result[i] << std::endl;
		}
		delete[] result;
		std::cout << duration << " Seconds" << std::endl;
	}
}
#endif // !TASK2_H
