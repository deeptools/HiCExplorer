// The HDF5 reading and writing and the thread pool come from coolercpp.

#ifndef HIC2COOL_IO_HPP
#define HIC2COOL_IO_HPP

#include <coolercpp/fastio.hpp>
#include <coolercpp/parallel.hpp>

namespace hic2cool::detail {

namespace h5 = coolercpp::fastio;
using coolercpp::ThreadPool;

}  // namespace hic2cool::detail

#endif  // HIC2COOL_IO_HPP
