#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "Linear.h"
#include "Simplex.h"
#include "Task.h"

bool next_combination(std::vector<size_t>& vec, size_t n) {
    size_t k = vec.size();
    if (n == k) {
        return false;
    }
    int i;
    for (i = k - 1; i >= 0 && vec[i] == n - 1; --i, --n) {}
    if (i < 0) {
        return false;
    }
    ++vec[i];
    for (size_t j = i + 1; j < k; ++j) {
        vec[j] = vec[j - 1] + 1;
    }
    return true;
}

std::vector<size_t> create_vector(size_t size) {
    std::vector<size_t> vec(size);
    for (size_t i = 0; i < size; ++i) {
        vec[i] = i;
    }
    return vec;
}

std::vector<std::vector<size_t>> combinations_of_cuts(std::vector<double> const& lengths, std::vector<size_t> const& products_numbers) {
    std::vector<std::vector<size_t>> all_combinations;
    return all_combinations;
}

int number_of_cuts(double length, double cur_length) {
    return (int)(length / cur_length);
}

void print(std::vector<int> const& combinations) {
    for (size_t i = 0; i < combinations.size() - 1; ++i) {
        if (combinations[i] != 0) {
            return;
        }
    }
    for (auto elem : combinations) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

std::vector<std::vector<int>> all_combs;

double mmin(std::vector<double> const& lengths, size_t index) {
    double min = lengths[index];

    for (size_t i = index + 1; i < lengths.size(); ++i) {
        if (min > lengths[i]) {
            min = lengths[i];
        }
    }
    return min;
}

void all_cuts(double L, std::vector<double> const& lengths, std::vector<int>& combinations, size_t index = 0) {
    if (index >= lengths.size()) {
        return;
    }
    double length = L;
    int noc = number_of_cuts(L, lengths[index]);
    combinations[index] = noc;

    while (combinations[index] > 0) {
        length = L - combinations[index] * lengths[index];

        if (length < mmin(lengths, index)) {
            /*for (auto elem : combinations) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;*/
            all_combs.push_back(combinations);
        }

        all_cuts(length, lengths, combinations, index + 1);
        --combinations[index];
    }
    length = L;

    all_cuts(length, lengths, combinations, index + 1);
}

Limitations create_limitations(std::vector<size_t> const& number_of_products) {
    size_t size = all_combs.size();
    size_t comb_size = all_combs[0].size();
    Limitations limitations;

    for (size_t i = 0; i < comb_size; ++i) {
        std::vector<double> limitation(size);
        for (size_t j = 1; j < size; ++j) {
            limitation[j - 1] = all_combs[j][i];
        }
        limitation[size - 1] = number_of_products[i];
        limitations.add_limitations({ limitation, LT::LT_GT });
    }

    for (size_t i = 0; i < limitations.limititation_size(); ++i) {
        /*for (auto elem : limitations.limitations[i].first) {
            std::cout << elem << " ";
        }*/
        std::cout << limitations.limitations[i].first[limitations.limitations[i].first.size() - 1];
        std::cout << std::endl;
    }

    return limitations;
}

double sum(std::vector<double> const& lengths, std::vector<size_t> const& number_of_products) {
    size_t size = lengths.size();
    if (size != number_of_products.size()) {
        return -1;
    }
    double sum = 0;
    for (size_t i = 0; i < size; ++i) {
        sum += lengths[i] * number_of_products[i];
    }
    return sum;
}

int matrix_rank(Matrix& mat) {
    size_t n = mat.get_n(), m = mat.get_m();
    double EPS = 1E-5;
    int rank = std::max(n, m);
    std::vector<char> line_used(n);
    for (size_t i = 0; i < m; ++i) {
        size_t j;
        for (j = 0; j < n; ++j)
            if (!line_used[j] && abs(mat[j][i]) > EPS)
                break;
        if (j == n)
            --rank;
        else {
            line_used[j] = true;
            for (size_t p = i + 1; p < m; ++p)
                mat[j][p] /= mat[j][i];
            for (size_t k = 0; k < n; ++k)
                if (k != j && abs(mat[k][i]) > EPS)
                    for (size_t p = i + 1; p < m; ++p)
                        mat[k][p] -= mat[j][p] * mat[k][i];
        }
    }
    return rank;
}

int main() {
    std::vector<size_t> vector = { 1, 2 ,3, 4, 5, 6 };
    std::vector<double> lengths{ 0.60, 0.68, 0.83, 1.61, 1.67, 1.79, 2.80, 3.25, 3.25, 3.70, 3.95 };
    //std::vector<size_t> numbers_of_products{ 249, 60, 97, 76, 72, 18, 43, 5424, 450, 515, 28 };
    std::vector<size_t> numbers_of_products{ 2 * 249, 2 * 60, 2 * 97, 2 * 76, 2 * 72, 2 * 18, 2 * 43, 2 * 5424, 2 * 450, 2 * 515, 2 * 28 };
    double constexpr L = 11.7;
    std::cout << "All: " << sum(lengths, numbers_of_products) << std::endl;

    // std::vector<double> lenghts{ 250, 200, 150 };
    // std::vector<size_t> numbers_of_products{ 249, 60, 97, 76, 72, 18, 43, 5424, 450, 515, 28 };
    // double constexpr L = 800;

    /*std::vector<int> combinations(lenghts.size());

    std::vector<int> first_comb(lenghts.size());
    all_combs.push_back(first_comb);
    all_cuts(L, lenghts, combinations);
    /*for (size_t i = 1; i < all_combs.size(); ++i) {
        for (auto elem : all_combs[i]) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }*/

    /*Limitations limitations = create_limitations(numbers_of_products);
    std::vector<double> obj_func(limitations.limitations[0].first.size() - 1);
    for (auto& elem : obj_func) {
        elem = 1;
    }
    std::vector<bool> var_signs(limitations.limitations[0].first.size() - 1);
    for (size_t i = 0; i < var_signs.size(); ++i) {
        var_signs[i] = true;
    }

    std::cout << "Counter))): " << all_combs.size() - 1 << std::endl;

    Linear linear(obj_func, limitations, var_signs);
    auto mat = linear.get_matrix();
    auto obj = linear.get_obj_func();
    auto b = linear.get_b();
    //std::cout << "Rank: " << matrix_rank(mat) << std::endl;
    std::cout << "Matrix params: " << mat.get_n() << " | " << mat.get_m() << std::endl;
    std::cout << "b: " << b.size() << std::endl;
    std::cout << "Func: " << obj.size() << std::endl;
    Simplex simplex(mat, linear.get_b(), linear.get_obj_func(), TT::TT_MIN);
    std::vector<double> optimal = simplex.answer_func();
    std::vector<double> opt = linear.back_to_original_vars(optimal);
    // std::vector<double> optimal = linear.solve_task();
    double sum = 0;
    double sum_int = 0;
    for (auto elem : opt) {
        //std::cout << "Elem: " << elem << std::endl;
        sum += elem;
        sum_int += std::ceil(elem);
    }*/
    // std::cout << "Answer: " << sum << " | " << sum_int << std::endl;

    Task task(L, lengths, numbers_of_products);
    auto answer = task.solve();

    std::cout << "Answer: " << answer.first << " | " << answer.second << std::endl;

    return 0;
}