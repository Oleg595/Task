#pragma once
#include <vector>
#include "Linear.h"

class Task {
public:
	Task(double L, std::vector<double> const& lengths, std::vector<size_t> const& number_of_products);
	std::pair<double, double> solve();

	~Task();

private:
	void all_cuts(double L, std::vector<double> const& lengths, std::vector<int>& combinations, size_t index = 0);
	Limitations create_limitations() const;

	static int number_of_cuts(double length, double cur_length);
	static double mmin(std::vector<double> lengths, size_t index);
	static int matrix_rang(Matrix& mat);

	std::vector<std::vector<int>> all_combs;
	std::vector<size_t> number_of_products;
};
