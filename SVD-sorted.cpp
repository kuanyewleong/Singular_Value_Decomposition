// SVD-sorted.cpp
//

#include "stdafx.h"
#include "sort.h"


int main(int argc, char* argv[]) {
	typedef double Real_t;
	
	// a fast approach to read from file and push_back into a vector
	char c; // to flush the commas
	double v,w,x,y,z;
	vector<double> data;

	ifstream file("test_output1.csv");

	// check whether if the file is available
	if (!file.is_open())
	{
		cout << "File not found!\n";
		return 1;
	}

	while (file >> v >> c >> w >> c >> x >> c >> y >> c >> z)
	{
		data.push_back(v);
		data.push_back(w);
		data.push_back(x);
		data.push_back(y);
		data.push_back(z);

	}

	// set size of the matrix, M x N (M dimensions, N samples)
	int MRow = 5, NCol = data.size() / MRow;
	int m = NCol, n = MRow; // swap the matrix size mxn, needed later	
	vector<double> all_sum, mean, Y;

	// compute means of N samples of each column M
	for (int i = 0; i < MRow; i++)
	{
		double sum = 0.0;
		int p = i;
		for (int j = 0; j < NCol; j++)
		{
			sum += data[p];
			p += MRow;
		}
		all_sum.push_back(sum);
	}

	for (int i = 0; i < MRow; i++)
	{
		mean.push_back(all_sum[i] / NCol);
	}

	// minus the mean from the corresponding data point	
	for (int i = 0; i < NCol; i++)
	{
		int p = i*MRow;
		for (int j = 0; j < MRow; j++)
		{
			data[p] = data[p] - mean[j];
			p++;
		}

	}
	
	// finally data samples are zero-centred, Y
	for (int i = 0; i < MRow*NCol; i++)
	{
		Y.push_back(data[i] / sqrt(NCol - 1));
	}

	// Transpose the Y vector into M matrix
	Real_t* M = new Real_t[m*n];
	for (int i = 0; i < n; i++)
	{
		int counter1 = i * m;
		int counter2 = 1 * i;
		for (int j = 0; j < m; j++)
		{
			M[counter1] = Y[counter2];
			counter1++;
			counter2 += 5;
		}
	}		

	// Get the row and column right
	int k = (m<n ? m : n);

	Real_t* u = new Real_t[m*k];
	Real_t* S = new Real_t[k];
	Real_t* PC = new Real_t[k*n];

	{ // Compute SVD
		int INFO = 0;
		char JOBU = 'S';
		char JOBVT = 'S';
		int wssize = 3 * (m<n ? m : n) + (m>n ? m : n);
		int wssize1 = 5 * (m<n ? m : n);
		wssize = (wssize>wssize1 ? wssize : wssize1);
		Real_t* wsbuf = new Real_t[wssize];		
		svd(&JOBU, &JOBVT, &m, &n, &M[0], &m, &S[0], &u[0], &m, &PC[0], &k, wsbuf, &wssize, &INFO);
		delete[] wsbuf;
	}

	{ // Check Error
		Real_t max_err = 0;
		for (size_t i0 = 0; i0<m; i0++)
			for (size_t i1 = 0; i1<n; i1++) {
				Real_t E = M[i1*m + i0];
				for (size_t i2 = 0; i2<k; i2++) E -= u[i2*m + i0] * S[i2] * PC[i1*k + i2];
				if (max_err<fabs(E)) max_err = fabs(E);
			}
		//cout << max_err << '\n';
	}

	
	// Push S into a vector
	vector<double> SV;
	for (int i = 0; i < n; i++)
		SV.push_back(S[i]);

	// Sort the singular values (SV) from large to small
	vector<size_t> sorted_indexes;
	vector<double> b;
	sort(SV, b, sorted_indexes);

	// SV after sorting
	cout << "\S = \n";
	for (int j = 0; j<SV.size(); j++)
	{		
		cout << sorted_indexes[j] << " " << SV[sorted_indexes[j]] << "\n";
	}
	
	// U is not needed here
		
	// Group the Principal Components (PC)
	vector<vector<double> > PC_groups(k);
			
	for (int i = 0; i < k; i++)
	{	
		int j1 = k * i;
		for (int j0 =0; j0 < k ; j0++)
		{				
			PC_groups[i].push_back(PC[j1]);
			j1++;
		}
	}
		
	// Sort PC's columns according to the sorted_indexes order
	// PC is transposed to PC' as well
	vector<double> sorted_PC;
	//cout << "\n\nPC' sorted = \n";
	for (int i = 0; i < k; i++)
	{		
		for (int j = 0; j < k; j++)
		{
			sorted_PC.push_back(PC_groups[j][sorted_indexes[i]]);			
		}
	}	// end of PC sorting

	/*
	// Show PC' and data
	int PC_counter = 0;
	for (int i = 0; i < n; i++)
	{		
		for (int j = 0; j < n; j++)
		{
			cout << sorted_PC[PC_counter] << " ";
			PC_counter++;
		}
		cout << "\n";
	}
	*/
	
	//cout << "Data = \n";
	// Transpose the data for projection calculation
	int data_counter = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			//cout << data[data_counter] << " ";
			data_counter++;
		}
		//cout << "\n";
	}	// end of data transpose

	// Calculate signal for projection
	// Signal = PC' * data (matrix multiplication, dot product of PC' . data)
	vector<double> signal;

	for (int h = 0; h < n; h++)
	{
		int data_iter_1 = 0;
		for (int i = 0; i <m; i++)
		{
			double product = 0.0;
			int PC_iter_1 = h * n;
			for (int j = 0; j < n; j++)
			{				
				product += sorted_PC[PC_iter_1] * data[data_iter_1];
				data_iter_1++;
				PC_iter_1++;
			}
			signal.push_back(product);
		}
	}	//end of projection calculation
	
	/*
	cout << "size of signal: " << signal.size() << "\n";
	
	for (int i = 0; i < signal.size(); i++)
		cout << signal[i] << " ";
	*/

	// ----- Dimensions Reduction ----- //
	// compute the total energy based on the sum of squares of the retained singular values
	
	// sum of squares of the fisrt 3 singular values as ordered in the sorted_indexes
	double SS_SV_1to3 = 0.0;
	for (int i = 0; i < 3; i++)
	{
		SS_SV_1to3 += pow(SV[sorted_indexes[i]], 2);
	}
	// sum of squares of all the singular values
	double SS_SV_all = 0.0;
	for (int i = 0; i < (int)SV.size(); i++)
	{
		SS_SV_all += pow(SV[i], 2);
	}
	// Percentage of energy in the first 3 sorted singular values
	double E_percentage = SS_SV_1to3 / SS_SV_all;
	// (we would like to retain enough singular values to make up at least or more than 90 % of the energy)
	cout.precision(4);
	if (E_percentage >= 0.9)
		cout << "Energy percentage is high: " << fixed << E_percentage << endl;
	else
	{
		cout << "Energy percentage is low: " << fixed << E_percentage << endl;
		return 1; // not possible to just use 3 dimensions for this data, hence termination is called
	}
	
	
	// Get points of the first 3 dimensions (to be used as X, Y, Z in 3D coordinate space)
	// and write the 3D points into a CSV file to be used later
	ofstream PC_file; PC_file.open("PC_3D.csv");	

	for (int i = 0; i < m; i++)
	{
		int p = i * 1;
		for (int j = 0; j < 3; j++)
		{
			cout << signal[p] << " ";
			PC_file << signal[p] << ",";
			p += m;
		}
		cout << "\n";
		PC_file << "\n";
	}
	PC_file.close();
	
	// clear the following from the memory	
	delete[] u;
	delete[] S;
	delete[] PC;
	delete[] M;

	return 0;
}
