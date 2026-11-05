// See hicx/dgeev_selected.hpp for why and what. The eigensolver is LAPACK's
// dgeev from the dependency prefix; this file only picks the requested
// columns from its result.

#include "hicx/dgeev_selected.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
// Fortran character arguments carry a hidden trailing length.
void dgeev_(const char*, const char*, const int*, double*, const int*, double*, double*, double*,
            const int*, double*, const int*, double*, const int*, int*, std::size_t, std::size_t);
}

namespace hicx {
namespace {

EigenResult to_result(int n, const std::vector<int>& which, const std::vector<double>& wr,
                      const std::vector<double>& wi,
                      const std::function<const double*(int)>& column,
                      std::vector<double>* values_imaginary) {
    EigenResult result;
    result.values.assign(which.size(), 0.0);
    result.vectors.assign(which.size(), {});
    if (values_imaginary != nullptr) {
        values_imaginary->assign(which.size(), 0.0);
    }
    for (std::size_t entry = 0; entry < which.size(); ++entry) {
        const int index = which[entry] - 1;
        if (index < 0 || index >= n) {
            continue;
        }
        // A conjugate pair occupies columns (j, j+1) as (real, imaginary);
        // scipy builds both complex vectors from them and hicPCA keeps the real
        // part, which is column j for both members.
        int source = index;
        if (wi[static_cast<std::size_t>(index)] < 0.0 && index > 0) {
            source = index - 1;
        }
        const double* values = column(source);
        result.values[entry] = wr[static_cast<std::size_t>(index)];
        result.vectors[entry].assign(values, values + n);
        if (values_imaginary != nullptr) {
            (*values_imaginary)[entry] = wi[static_cast<std::size_t>(index)];
        }
    }
    return result;
}

EigenResult selected_impl(int n, double* a, const std::vector<int>& which,
                          std::vector<double>* values_imaginary, const std::function<void()>& release) {
    if (n == 0) {
        release();
        EigenResult empty;
        empty.values.assign(which.size(), 0.0);
        empty.vectors.assign(which.size(), {});
        if (values_imaginary != nullptr) {
            values_imaginary->assign(which.size(), 0.0);
        }
        return empty;
    }
    const char cV = 'V';
    const char cN = 'N';
    std::vector<double> wr(static_cast<std::size_t>(n));
    std::vector<double> wi(static_cast<std::size_t>(n));
    int info = 0;
    int lwork = -1;
    double query = 0.0;
    double dummy = 0.0;
    dgeev_(&cN, &cV, &n, a, &n, wr.data(), wi.data(), &dummy, &n, &dummy, &n, &query, &lwork,
           &info, 1, 1);
    if (info != 0) {
        throw std::runtime_error("dgeev workspace query failed");
    }
    lwork = static_cast<int>(query);
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    std::vector<double> vr(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
    dgeev_(&cN, &cV, &n, a, &n, wr.data(), wi.data(), &dummy, &n, vr.data(), &n, work.data(),
           &lwork, &info, 1, 1);
    if (info != 0) {
        throw std::runtime_error("dgeev failed with info " + std::to_string(info));
    }
    work.clear();
    work.shrink_to_fit();
    release();
    return to_result(
        n, which, wr, wi,
        [&](int column) {
            return vr.data() + static_cast<std::size_t>(column) * static_cast<std::size_t>(n);
        },
        values_imaginary);
}

}  // namespace

EigenResult dgeev_selected_eigenvectors(DenseSymmetric& matrix, const std::vector<int>& which) {
    const int n = static_cast<int>(matrix.size());
    return selected_impl(n, matrix.data(), which, nullptr,
                         [&matrix] { matrix.release(); });
}

EigenResult dgeev_selected_general(int n, double* matrix, const std::vector<int>& which,
                                   std::vector<double>* values_imaginary) {
    return selected_impl(n, matrix, which, values_imaginary, [] {});
}

}  // namespace hicx
