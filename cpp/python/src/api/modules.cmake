# The hicexplorer modules whose public names hicexplorer._hicx implements,
# one per line, sorted (cpp/python/API_BINDINGS.md).
#
# An entry is the module path below hicexplorer with '.' replaced by '_'
# (hicexplorer.lib.cnb is lib_cnb). Its bindings are src/api/<entry>.cpp,
# which defines
#
#     void hicx::python::bind_<entry>(pybind11::module_& module);
#
# and the extension offers them as hicexplorer._hicx.<entry>.
set(HICX_API_MODULES
    chicAggregateStatistic
    chicDifferentialTest
    chicExportData
    chicQualityControl
    chicSignificantInteractions
    chicViewpoint
    chicViewpointBackgroundModel
    lib_viewpoint
    readBed
    utilities
)
