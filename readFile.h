#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <iterator>

using namespace std;

class Read {
public:
	void readCSV(istream &input, vector< vector<string> > &output)
	{
		string csvLine;
		// read every line from the stream
		while (getline(input, csvLine))
		{
			istringstream csvStream(csvLine);
			vector<string> csvColumn;
			string csvElement;
			// read every element from the line that is separated by commas
			// and put it into the vector 
			while (getline(csvStream, csvElement, ','))
			{
				csvColumn.push_back(csvElement);
			}

			output.push_back(csvColumn);
		}
	}
};