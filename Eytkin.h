#pragma once
#include <iostream>

#include <string.h>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace std;
class Eytkin
{
private:
double start, stop, step;
int qknots;
double* x, * y,*nex,*ney;
double** EytkinMatr;

public:

	Eytkin()
	{
		qknots = 11;
		start = 0;
		stop = 0.5;
		step = (stop - start) / 10;

		x = new double[qknots];
		y = new double[qknots];

		EytkinMatr = new double* [5];
		for (int i = 0; i < 5; i++)
		{
			EytkinMatr[i] = new double[9];
		}

		

		ZeroingMatr();
		ZeroingMas(x, qknots);
		ZeroingMas(y, qknots);
		FunctionCalculation();
	}
	void ZeroingMatr()
	{
		for (int i = 0; i < 5; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				EytkinMatr[i][j] = 0;
			}
		}
	}
	void ZeroingMas(double* mas, int n)
	{
		for (int i = 0; i < n; i++)
		{
			mas[i] = 0;
		}
	}


	void FunctionCalculation()
	{
		x[0] = start;
		y[0] = sqrt(asin(x[0]));
		for (int i = 1; i < qknots; i++)
		{
			x[i] = x[i - 1] + step;
			y[i] = sqrt(asin(x[i]));
		}

		cout << "X" << "\t" << "Y" << endl;
		cout << setprecision(4);
		for (int i = 0; i < qknots; i++)
		{
			cout << x[i] << "\t" << y[i] << endl;
		}
	}


	double FillingEytkinTable(double X)
	{

		cout <<"X = "<<X<<endl;
	
		nex = new double[5];
		nex[0] = x[0];
		nex[1] = x[2];
		nex[2] = x[5];
		nex[3] = x[7];
		nex[4] = x[10];
		ney = new double[5];
		ney[0] = y[0];
		ney[1] = y[2];
		ney[2] = y[5];
		ney[3] = y[7];
		ney[4] = y[10];

		for (int i = 0; i < 5; i++)
		{
			for (int j = 1; j < 6; j++)
			{
				if (i == (j - 1))
				{
					EytkinMatr[i][j] = round((X - nex[i]) * 10000) / 10000;
					
				}
				else
				{
					EytkinMatr[i][j] = round((nex[i] - nex[j - 1]) * 100000) / 100000;
					
				}
			}
		}


		for (int i = 0; i < 5; i++)
		{
			double tmp = 1;
			for (int j = 1; j < 6; j++)
			{
				
				tmp *= round(EytkinMatr[i][j]*10000)/10000;

			}

			EytkinMatr[i][6] = tmp;
			EytkinMatr[i][7] = ney[i];
			
			if(EytkinMatr[i][6]==0)
			{
				EytkinMatr[i][7] = 0;
			}
			else if(EytkinMatr[i][6] != 0)
			{
				EytkinMatr[i][8] = round((EytkinMatr[i][7] / EytkinMatr[i][6])*10000)/10000;
			}
		
			
			EytkinMatr[i][0] = i;

		}

		return X;
	}

	
	

	double DoubleLagrangPolynom()
	{
		double mult = 1;
		double sum = 0;
		double polynom=0;
		for (int i = 0; i < 5; i++)
		{
			for (int j = 1; j < 6; j++)
			{
				if (i == (j - 1))
				{
					mult *= EytkinMatr[i][j];
				}
			}
		}
		
		for (int i = 0; i < 5; i++)
		{

			sum += EytkinMatr[i][8];

		}
	
		cout << "Sum=" << sum<<endl;
		cout << "Proiz" << mult << endl;
		polynom = mult * sum;
	

		return polynom;
	}

	void CallLagrang()
	{
		double** FinaleTable;
		FinaleTable = new double* [11];
		for (int i = 0; i < 11; i++)
		{
			FinaleTable[i] = new double[5];
		}


		FinaleTable[1][0] = FillingEytkinTable(x[2] + 1.5 * step);
		FinaleTable[1][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[1][1] = sqrt(asin((x[2] + 1.5 * step)));
		FinaleTable[1][3] = abs(FinaleTable[1][1] - FinaleTable[1][2]);
		FinaleTable[1][4] = (FinaleTable[1][3] / FinaleTable[1][1]) * 100;


		FinaleTable[2][0] = FillingEytkinTable(x[2] + 1.5 * step);
		FinaleTable[2][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[2][1] = sqrt(asin((x[2] + 1.5 * step)));
		FinaleTable[2][3] = abs(FinaleTable[2][1] - FinaleTable[2][2]);
		FinaleTable[2][4] = (FinaleTable[2][3] / FinaleTable[2][1]) * 100;


		FinaleTable[3][0]= FillingEytkinTable(x[6]);
		FinaleTable[3][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[3][1] = sqrt(asin((x[6])));
		FinaleTable[3][3] = abs(FinaleTable[3][1] - FinaleTable[3][2]);
		FinaleTable[3][4] = (FinaleTable[3][3] / FinaleTable[3][1]) * 100;


		FinaleTable[4][0] = FillingEytkinTable(x[10]-0.5*step);
		FinaleTable[4][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[4][1] = sqrt(asin((x[10]-0.5*step)));
		FinaleTable[4][3] = abs(FinaleTable[4][1] - FinaleTable[4][2]);
		FinaleTable[4][4] = (FinaleTable[4][3] / FinaleTable[4][1]) * 100;

		FinaleTable[5][0] = FillingEytkinTable(x[8] - 2.5 * step);
		FinaleTable[5][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[5][1] = sqrt(asin((x[8] - 2.5 * step)));
		FinaleTable[5][3] = abs(FinaleTable[5][1] - FinaleTable[5][2]);
		FinaleTable[5][4] = (FinaleTable[5][3] / FinaleTable[5][1]) * 100;

		FinaleTable[6][0] = FillingEytkinTable(x[3] + 0.2*step);
		FinaleTable[6][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[6][1] = sqrt(asin((x[3] + 0.2*step)));
		FinaleTable[6][3] = abs(FinaleTable[6][1] - FinaleTable[6][2]);
		FinaleTable[6][4] = (FinaleTable[6][3] / FinaleTable[6][1]) * 100;


		FinaleTable[7][0] = FillingEytkinTable(x[5] - 3.5*step);
		FinaleTable[7][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[7][1] = sqrt(asin((x[5] - 3.5 * step)));
		FinaleTable[7][3] = abs(FinaleTable[7][1] - FinaleTable[7][2]);
		FinaleTable[7][4] = (FinaleTable[7][3] / FinaleTable[7][1]) * 100;

		FinaleTable[8][0] = FillingEytkinTable(x[7] +0.3 * step);
		FinaleTable[8][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[8][1] = sqrt(asin((x[7] + 0.3 * step)));
		FinaleTable[8][3] = abs(FinaleTable[8][1] - FinaleTable[8][2]);
		FinaleTable[8][4] = (FinaleTable[8][3] / FinaleTable[8][1]) * 100;


		FinaleTable[9][0] = FillingEytkinTable(x[1] +0.9 * step);
		FinaleTable[9][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[9][1] = sqrt(asin((x[1] + 0.9 * step)));
		FinaleTable[9][3] = abs(FinaleTable[9][1] - FinaleTable[9][2]);
		FinaleTable[9][4] = (FinaleTable[9][3] / FinaleTable[9][1]) * 100;

		FinaleTable[10][0] = FillingEytkinTable(x[10] - 0.9 * step);
		FinaleTable[10][2] = DoubleLagrangPolynom();
		Printf(EytkinMatr, 9,5);
		FinaleTable[10][1] = sqrt(asin((x[10] - 0.9 * step)));
		FinaleTable[10][3] = abs(FinaleTable[10][1] - FinaleTable[10][2]);
		FinaleTable[10][4] = (FinaleTable[10][3] / FinaleTable[10][1]) * 100;

		cout << endl;
		cout << "X\t" << "   " << "Истинное\t"<< "Лагранж\t" << " " << "Абсолютная"<<" Относительная" <<endl;

		Printf(FinaleTable, 5,10);

		

	}
	void Printf(double** mas, int m,int n)
	{
		
		for (int i = 0; i <n; i++)
		{

			
			cout << endl;
			cout << setprecision(10);
			for (int j = 0; j <m; j++)
			{
				
				
					
				
			printf("%1.7f  ", mas[i][j]);
			
			}
		}

		cout << endl;
		cout << "\n";
	}


};

