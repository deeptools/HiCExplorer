#ifndef HICFILECPP_DETAIL_NORMS_HPP
#define HICFILECPP_DETAIL_NORMS_HPP

// The normalization vectors and expected values Juicer tools compute, ported
// from juicebox.tools.utils.norm.NormalizationCalculations, ZeroScale,
// final2.FinalScale and juicebox.tools.utils.original.ExpectedValueCalculation.
//
// Every vector computation is a template over the precision of its vectors:
// double reproduces Juicer tools 1.22.01, which writes version 8 files, and
// float reproduces Juicer tools 2.20.00, which writes version 9 files and
// keeps its normalization vectors in float arrays. Scalars and accumulations
// stay double in both, as in Java.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace hicfilecpp::detail {

// One stored record of an intra-chromosomal matrix, in bins, x <= y.
struct NormRecord {
    int32_t x = 0;
    int32_t y = 0;
    float counts = 0.0f;
};

// commons-math 2 StatUtils.percentile (Juicer tools 1.22.01) and commons-math3
// DescriptiveStatistics.getPercentile with its default LEGACY estimation
// (2.20.00). p is a percentage in (0, 100].
double percentileMath2(std::vector<double> values, double p);
double percentileMath3(std::vector<double> values, double p);

class ExpectedValueCalculation {
public:
    // lengths are indexed by chromosome index; index 0, the whole-genome
    // pseudo-chromosome, takes no part.
    ExpectedValueCalculation(const std::vector<int64_t>& lengths, int32_t gridSize, std::string type);

    void addDistance(int32_t chrIndex, int32_t bin1, int32_t bin2, double weight);
    bool hasData() const { return !chromosomeCounts_.empty(); }
    void computeDensity();

    int32_t gridSize() const { return gridSize_; }
    const std::string& type() const { return type_; }
    const std::vector<double>& densityAvg() const { return densityAvg_; }
    const std::map<int32_t, double>& chrScaleFactors() const { return chrScaleFactors_; }

private:
    std::vector<int64_t> lengths_;
    int32_t gridSize_;
    std::string type_;
    int64_t numberOfBins_ = 0;
    std::vector<double> actualDistances_;
    std::map<int32_t, double> chromosomeCounts_;
    std::vector<double> densityAvg_;
    std::map<int32_t, double> chrScaleFactors_;
};

// NormalizationCalculations.computeVC
template <class T>
std::vector<T> computeVC(const std::vector<NormRecord>& records, int64_t size);

// NormalizationCalculations.getSumFactor
template <class T>
double sumFactor(const std::vector<NormRecord>& records, const std::vector<T>& norm);

// NormalizationCalculations.computeKR; math3 selects the percentile of 2.20.00.
template <class T>
std::vector<T> computeKR(const std::vector<NormRecord>& records, int64_t size, bool math3);

// NormalizationCalculations.computeMMBA: FinalScale followed by
// ZeroScale.normalizeVectorByScaleFactor.
template <class T>
std::vector<T> computeScale(const std::vector<NormRecord>& records, int64_t size);

// ExpectedValueCalculation.addDistancesFromIterator
template <class T>
void addDistancesFromRecords(ExpectedValueCalculation& ev, int32_t chrIndex, const std::vector<NormRecord>& records,
                             const std::vector<T>& vector);

extern template std::vector<double> computeVC(const std::vector<NormRecord>&, int64_t);
extern template std::vector<float> computeVC(const std::vector<NormRecord>&, int64_t);
extern template double sumFactor(const std::vector<NormRecord>&, const std::vector<double>&);
extern template double sumFactor(const std::vector<NormRecord>&, const std::vector<float>&);
extern template std::vector<double> computeKR(const std::vector<NormRecord>&, int64_t, bool);
extern template std::vector<float> computeKR(const std::vector<NormRecord>&, int64_t, bool);
extern template std::vector<double> computeScale(const std::vector<NormRecord>&, int64_t);
extern template std::vector<float> computeScale(const std::vector<NormRecord>&, int64_t);
extern template void addDistancesFromRecords(ExpectedValueCalculation&, int32_t, const std::vector<NormRecord>&,
                                             const std::vector<double>&);
extern template void addDistancesFromRecords(ExpectedValueCalculation&, int32_t, const std::vector<NormRecord>&,
                                             const std::vector<float>&);

}  // namespace hicfilecpp::detail

#endif  // HICFILECPP_DETAIL_NORMS_HPP
