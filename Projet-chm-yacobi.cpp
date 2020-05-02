#include "Projet-chm-yacobi.h"

template<typename T>
class Matrix {
public:
    Matrix(std::size_t size) : elements(size, std::vector<T>(size)) {}

    Matrix(std::vector<std::vector<T>> &&v) : elements(v) {}

    Matrix<T> getWithoutRowAndCol(std::size_t ni, std::size_t nj) {
        if (elements.size() == 1) {
            return Matrix<T>(0, 0);
        }
        std::vector<std::vector<T>> result;
        for (std::size_t i = 0; i < elements.size(); ++i) {
            if (ni == i) continue;
            result.push_back({});
            for (std::size_t j = 0; j < elements[i].size(); ++j) {
                if (nj == j) continue;
                result.back().push_back(elements[i][j]);
            }
        }
        return std::move(result);
    }

    std::vector<std::vector<T>> const &getValue() const {
        return elements;
    }

    void addElement(std::size_t i, std::size_t j, T element) {
        if (i >= elements.size() || j >= elements[i].size()) return;
        elements[i][j] = element;
    }

    void print() const {
        std::cout.precision(PRECISION_PRINT);
        std::cout << "[" << std::endl;
        bool is_first_row = true;
        for (auto const &v : elements) {
            if (!is_first_row) std::cout << "," << std::endl;
            std::cout << "    [";
            bool is_first = true;
            for (auto const &el : v) {
                if (!is_first) std::cout << ", ";
                std::cout << el;
                is_first = false;
            }
            std::cout << "]";
            is_first_row = false;
        }
        std::cout << std::endl << "]" << std::endl;
    }

    Matrix<T> getB() const {
        std::vector<std::vector<T>> result(elements.size(), std::vector<T>(elements.size()));
        for (std::size_t i = 0; i < elements.size(); ++i) {
            for (std::size_t j = 0; j < elements[0].size(); ++j) {
                result[i][j] = (i == j ? 0 : -elements[i][j] / elements[i][i]);
            }
        }
        return std::move(result);
    }

    std::vector<T> getVector(std::vector<T> const &b) const {
        std::vector<T> result(elements.size());
        for (std::size_t i = 0; i < elements.size(); ++i) {
            result[i] = b[i] / elements[i][i];
        }
        return result;
    }

    std::vector<T> prod(std::vector<T> const &v) {
        std::vector<T> result(elements.size());
        for (std::size_t i = 0; i < elements.size(); ++i) {
            for (std::size_t j = 0; j < elements.size(); ++j) {
                result[i] += v[j] * elements[i][j];
            }
        }
        return result;
    }

private:
    std::vector<std::vector<T>> elements;
};

bool is_zeros(const std::vector<type>& vec)
{
    bool res = true;
    for (size_t i = 0; i < vec.size(); ++i)
    {
        res &= std::abs(vec[i]) < ZERO_EPS;
    }
    std::cout << "Check: " << res << std::endl;
    return res;
}

type dist_vectors(const std::vector<type>& first, const std::vector<type>& second)
{
    type res = 0;
    for (size_t i = 0; i < first.size(); ++i)
    {
        type first_i = i >= first.size() ? 0 : first[i];
        type second_i = i >= second.size() ? 0 : second[i];
        type delta = first_i - second_i;
        res = std::max(res, std::abs(delta)); // += delta * delta;
    }
    return res;
}

template<typename T>
void sum(std::vector<T> const &first, std::vector<T> const &second, std::vector<T> &dst) {
    for (std::size_t i = 0; i < dst.size(); ++i) {
        dst[i] = first[i] + second[i];
    }
}

template<typename T = type>
Matrix<T> generateDiagMatrix(std::size_t size) {
    std::vector<std::vector<T>> result(size, std::vector<T>(size));
    static std::mt19937 gen(SEED); // rd()
    static std::uniform_real_distribution<> dis(0., 30);
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            result[i][j] = (i == j ? 0 : dis(gen));
        }
    }
    for (std::size_t i = 0; i < size; ++i) {
        double rowSum = 0;
        for (auto const &elem : result[i]) { rowSum += elem; }
        result[i][i] = rowSum + dis(gen) + 1.;
    }
    return std::move(result);
}

template<typename T = type>
Matrix<T> generateRandomMatrix(std::size_t size) {
    std::vector<std::vector<T>> result(size, std::vector<T>(size));
    //std::random_device rd;
    static std::mt19937 gen(SEED); // rd()
    static std::uniform_real_distribution<> dis(0., 1.);
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            result[i][j] = dis(gen);
        }
    }
    return std::move(result);
}

template<typename T = type>
Matrix<T> generateGilbertMatrix(std::size_t size) {
    std::vector<std::vector<T>> result(size, std::vector<T>(size));
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            result[i][j] = 1.0 / (1.0 + i + j);
        }
    }
    return std::move(result);
}

template<typename T = type>
std::vector<T> generateRandomVector(std::size_t size) {
    std::vector<T> result(size);
    static std::mt19937 gen(SEED); // rd()
    static std::uniform_real_distribution<> dis(0., 200.);
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = dis(gen);
    }
    return result;
}

template <typename T>
Matrix<T> readMatrix(size_t size)
{
    std::vector<std::vector<T>> result(size, std::vector<T>(size));
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            std::cin >> result[i][j];
        }
    }
    return std::move(result);
}

template<typename T>
void printVector(char const *preffix, std::vector<T> const &vec, std::size_t ind = 0) {
    std::cout.precision(PRECISION_PRINT);
    std::cout << preffix;
    if (ind > 0) std::cout << ind;
    std::cout << " = [";
    bool is_first = true;
    for (auto const &el : vec) {
        if (!is_first) std::cout << ", ";
        std::cout << el;
        is_first = false;
    }
    std::cout << "]" << std::endl;
}

template<typename T>
void yacobi_iter(Matrix<T> const &a, std::vector<T> const &b, std::size_t cntIter = 10000) {
    Matrix<type> B = a.getB();
    std::vector<type> c = a.getVector(b);
    std::cout << "B matrix = " << std::endl;
    B.print();
    printVector("c", c);
    std::vector<type> xs = (is_zeros(c) ? std::vector<type>(b.size(), 0.5) : c);
    for (std::size_t iter = 1; iter <= cntIter; ++iter) {
        printVector("xs_", xs, iter);
        std::vector<type> prev = xs;
        sum(B.prod(xs), c, xs);
        if (std::abs(dist_vectors(prev, xs)) <= EPS)
        {
            std::cout << "Finished by EPS" << std::endl;
            break;
        }
    }
}

int main() {
    std::size_t size = 10;
    std::vector<double> b = {2.656991643595489894e-01, 5.803798583735642058e-02, 5.624955118666850051e-01 ,1.931421196895460879e-01, 5.775747262761687928e-02, 4.643785469426497947e-01, 3.272472661949272776e-01, 7.331189716378443411e-01, 8.798678821028841357e-01, 9.433030545066206640e-01};
    std::cout << "Good (diag) matrix:" << std::endl;
    {
        Matrix<type> a = readMatrix<type>(size);
        a.print();
        printVector("b", b);
        yacobi_iter(a, b);
    }
    std::cout << "---------------------------- end ----------------------------" << std::endl;
    std::cout << "Normal (random) matrix:" << std::endl;
    {
        Matrix<type> a = readMatrix<type>(size);
        a.print();
        printVector("b", b);
        yacobi_iter(a, b, 200);
    }
    std::cout << "---------------------------- end ----------------------------" << std::endl;
    std::cout << "Bad (gilbert) matrix:" << std::endl;
    {
        Matrix<type> a = readMatrix<type>(size);
        a.print();
        printVector("b", b);
        std::cout << std::endl << "Execute algo" << std::endl;
        yacobi_iter(a, b, 200);
    }
    return 0;
}