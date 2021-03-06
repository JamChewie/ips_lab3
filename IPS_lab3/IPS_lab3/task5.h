#ifndef TASK5_H
#define TASK5_H
#include <iostream>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
namespace TASK5
{

	// ���������� ����� � �������� ���������� �������
	const int MATRIX_SIZE = 1500;

	/// ������� InitMatrix() ��������� ���������� � �������� 
	/// ��������� ���������� ������� ���������� ����������
	/// matrix - �������� ������� ����
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

	/// ������� SerialGaussMethod() ������ ���� ������� ������ 
	/// matrix - �������� ������� �������������� ���������, �������� � ����,
	/// ��������� ������� ������� - �������� ������ ������ ���������
	/// rows - ���������� ����� � �������� �������
	/// result - ������ ������� ����
	double SerialGaussMethod(double **matrix, const int rows, double* result)
	{
		int k;
		double koef;
		// ������ ��� ������ ������
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
		// �������� ��� ������ ������
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

	double ParallelGaussMethod(double **matrix, const int rows, double* result)
	{
		int k;		
		// ������ ��� ������ ������
		auto start = clock() / 1000.0;
		for (k = 0; k < rows; ++k)
		{
			//
			cilk_for(int i = k + 1; i < rows; ++i)
			{
				double koef = -matrix[i][k] / matrix[k][k];

				for (int j = k; j <= rows; ++j)
				{
					matrix[i][j] += koef * matrix[k][j];
				}
			}
		}
		auto end = clock() / 1000.0;
		auto duration = end - start;
		// �������� ��� ������ ������
		result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];
		for (k = rows - 2; k >= 0; --k)
		{
			cilk::reducer_opadd<double> result_k(matrix[k][rows]);

			cilk_for(int j = k + 1; j < rows; ++j)
			{
				//result[k] -= matrix[k][j] * result[j];
				result_k -= matrix[k][j] * result[j];
			}

			result[k] = result_k.get_value() / matrix[k][k];
		}
		return duration;
	}

	// ������ � �������� ��������
	void run_task5_test_matrix()
	{
		srand((unsigned)time(0));

		int i;

		// ���-�� ����� � �������, ���������� � �������� �������
		const int test_matrix_lines = 4;

		double **test_matrix = new double*[test_matrix_lines];

		// ���� �� �������
		for (i = 0; i < test_matrix_lines; ++i)
		{
			// (test_matrix_lines + 1)- ���������� �������� � �������� �������,
			// ��������� ������� ������� ������� ��� ������ ����� ���������, �������� � ����
			test_matrix[i] = new double[test_matrix_lines + 1];
		}

		// ������ ������� ����
		double *result = new double[test_matrix_lines];

		// ������������� �������� �������
		test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
		test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
		test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
		test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;
		auto durationSerial = SerialGaussMethod(test_matrix, test_matrix_lines, result);
		auto durationParallel = ParallelGaussMethod(test_matrix, test_matrix_lines, result);

		for (i = 0; i < test_matrix_lines; ++i)
		{
			delete[]test_matrix[i];
		}

		delete[] result;
		std::cout << "Duration serial = " << durationSerial << std::endl
			<< "Duration parallel = " << durationParallel << std::endl;
	}
	 
	// ������ � �������� ������� MATRIX_SIZE
	void run_task5_MATRIX_SIZE()
	{
		srand((unsigned)time(0));

		double **test_matrix = new double*[MATRIX_SIZE];

		// ������ ������� ����
		double *result = new double[MATRIX_SIZE];

		// ������������� �������
		InitMatrix(test_matrix);

		auto durationSerial = SerialGaussMethod(test_matrix, MATRIX_SIZE, result);
		auto durationParallel = ParallelGaussMethod(test_matrix, MATRIX_SIZE, result);
		
		for (auto i = 0; i < MATRIX_SIZE; ++i)
		{
			delete[]test_matrix[i];
		}			
		delete[] result;
		std::cout << "Duration serial = " << durationSerial << " seconds" << std::endl
			  << "Duration parallel = " << durationParallel << " seconds" << std::endl;
	}
	
}
#endif // !TASK5_H
