#include "Task.h"
#include "Simplex.h"

Task::Task(double L, std::vector<double> const& lengths, std::vector<size_t> const& number_of_products) : number_of_products{ number_of_products } {
    all_combs.push_back(std::vector<int>(lengths.size()));
    std::vector<int> combinations(lengths.size());
    all_cuts(L, lengths, combinations);
}

Task::~Task() {}

void Task::all_cuts(double L, std::vector<double> const& lengths, std::vector<int>& combinations, size_t index) {
    if (index >= lengths.size()) {
        return;
    }
    double length = L;
    int noc = number_of_cuts(L, lengths[index]);
    combinations[index] = noc;

    while (combinations[index] > 0) {
        length = L - combinations[index] * lengths[index];

        if (length < mmin(lengths, index)) {
            all_combs.push_back(combinations);
        }

        all_cuts(length, lengths, combinations, index + 1);
        --combinations[index];
    }
    length = L;

    all_cuts(length, lengths, combinations, index + 1);
}

int Task::number_of_cuts(double length, double cur_length) {
    return (int)(length / cur_length);
}

double Task::mmin(std::vector<double> lengths, size_t index) {
    double min = lengths[index];

    for (size_t i = index + 1; i < lengths.size(); ++i) {
        if (min > lengths[i]) {
            min = lengths[i];
        }
    }
    return min;
}

int Task::matrix_rang(Matrix& mat) {
    size_t n = mat.get_n(), m = mat.get_m();
    double EPS = 1E-5;
    int rang = std::max(n, m);
    std::vector<char> line_used(n);
    for (size_t i = 0; i < m; ++i) {
        size_t j;
        for (j = 0; j < n; ++j)
            if (!line_used[j] && abs(mat[j][i]) > EPS)
                break;
        if (j == n)
            --rang;
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
    return rang;
}

Limitations Task::create_limitations() const {
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

    return limitations;
}

std::pair<double, double> Task::solve() {
    Limitations limitations = create_limitations();
    std::vector<double> obj_func(limitations.limitations[0].first.size() - 1);
    for (auto& elem : obj_func) {
        elem = 1;
    }
    std::vector<bool> var_signs(limitations.limitations[0].first.size() - 1);
    for (size_t i = 0; i < var_signs.size(); ++i) {
        var_signs[i] = true;
    }

    Linear linear(obj_func, limitations, var_signs);
    auto mat = linear.get_matrix();
    auto obj = linear.get_obj_func();
    auto b = linear.get_b();
    // std::cout << "Rank: " << matrix_rank(mat) << std::endl;
    // std::cout << "Matrix params: " << mat.get_n() << " | " << mat.get_m() << std::endl;
    // std::cout << "b: " << b.size() << std::endl;
    // std::cout << "Func: " << obj.size() << std::endl;
    Simplex simplex(mat, linear.get_b(), linear.get_obj_func(), TT::TT_MIN);
    std::vector<double> optimal = simplex.answer_func();
    std::vector<double> opt = linear.back_to_original_vars(optimal);
    double sum = 0;
    double sum_int = 0;
    for (auto elem : opt) {
        //std::cout << "Elem: " << elem << std::endl;
        sum += elem;
        sum_int += std::ceil(elem);
    }

    return { sum, sum_int };
}
