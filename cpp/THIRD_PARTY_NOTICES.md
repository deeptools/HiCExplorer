# Third-party code in HiCExplorer 4 (cpp/)

This file lists the code of other projects that is part of the C++ tree,
where it came from, and the licence it is under. The licence texts are
reproduced in full; the directories holding the code keep their own copy.

## coolercpp

- Location: `core/io/coolercpp/` (headers in `include/coolercpp/`, sources in `src/`), compiled into libhicx4core.
- Origin: coolercpp, git@github.com:MHH-Bioinformatics-Hematology/coolercpp.git (local clone `~/src/coolercpp`), commit `2078ff1fc038a999fd30013968c3c742a251e63e` ("Merge fastio: threaded chunk io in the library"), version 0.4.0. This is the commit HiCExplorer pinned (`HICX_COOLERCPP_COMMIT`); the repository HEAD at the time of vendoring (`a90fdeb`) is 12 commits newer and was not used.
- Changes: `include/coolercpp/version.hpp.in` became `version.hpp` with the project version written out. Tests, documentation, examples, the equivalence harness and the CMake packaging were not copied.
- Licence: BSD 3-Clause, Copyright (c) 2026, Joachim Wolff. The licence text was added to coolercpp after the pinned commit and is taken from its HEAD.
- coolercpp contains code ported from cooler 0.10.2 (https://github.com/open2c/cooler), BSD 3-Clause, Copyright (c) 2015-2023, Cooler developers. The ported files are `api.hpp`/`api.cpp`, `rangequery.hpp`/`rangequery.cpp`, `create.hpp`/`create.cpp`, `region.hpp`/`region.cpp` and `util.hpp`/`util.cpp`. Behaviour reproduced from h5py, numpy, pandas and scipy (all BSD 3-Clause) is implemented from their documented semantics; no code of theirs is included.

```
BSD 3-Clause License

Copyright (c) 2015-2023, Cooler developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

The coolercpp licence itself (`core/io/coolercpp/LICENSE`):

```
BSD 3-Clause License

Copyright (c) 2026, Joachim Wolff
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## hicfilecpp (repository hiccpp)

- Location: `core/io/hicfilecpp/` (headers in `include/hicfilecpp/`, sources in `src/`), compiled into libhicx4core.
- Origin: hiccpp, git@github.com:MHH-Bioinformatics-Hematology/hiccpp.git (local clone `~/src/hiccpp`), commit `15a2aedfb24e51d6c2d44b1833ed7cb78115215a` ("Write caller-supplied normalization vectors, report a missing normalized section, 0.4.0"), library version 0.4.0. This is the commit HiCExplorer pinned (`HICX_HICFILECPP_COMMIT`); the repository HEAD at the time of vendoring (`9af9af6`) is 9 commits newer, renames the library to `hiccpp` and was not used.
- Changes: `include/hicfilecpp/version.hpp.in` became `version.hpp` with the project version written out. Tests, documentation, examples, the equivalence harness and the CMake packaging were not copied.
- Licence: MIT, Copyright (c) 2026 Joachim Wolff. The licence text was added to hiccpp after the pinned commit and is taken from its HEAD.
- hiccpp contains code ported from Juicer tools and Juicebox 1.22.01 and 2.20.00 (https://github.com/aidenlab/Juicebox, MIT, Copyright (c) 2016 Aiden Lab): the `.hic` writer (`src/writer.cpp`, from `Preprocessor`, `MatrixPP`, `MatrixZoomDataPP` and `NormalizationVectorUpdater`) and the normalization vectors and expected values (`src/norms.hpp`/`src/norms.cpp`, from `NormalizationCalculations`, `ZeroScale`, `FinalScale` and `ExpectedValueCalculation`). The reader (`src/reader.cpp`) reproduces the query semantics of straw/hicstraw 1.3.1 (https://github.com/aidenlab/straw, MIT, Copyright (c) 2016 Aiden Lab). hic2cool (https://github.com/4dn-dcic/hic2cool, MIT, Harvard Medical School) served as a second reading of the version 8 layout. The licence of Juicebox and straw is reproduced below the hiccpp licence.

```
MIT License

Copyright (c) 2026 Joachim Wolff

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

The licence of Juicebox and of straw (identical text in both repositories):

```
The MIT License (MIT)

Copyright (c) 2016 Aiden Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## hic2cool (repository hicxchange)

- Location: `core/io/hic2cool/` (header in `include/hic2cool/`, sources in `src/`), compiled into libhicx4core. The `.hic` to cool and cool to `.hic` conversion is reachable only through `hicConvertFormat`; no converter executable is built.
- Origin: hicxchange, git@github.com:MHH-Bioinformatics-Hematology/hicxchange.git (local clone `~/src/hicxchange`), commit `48d0042042f3a1d823847c7096d9afba37fed15d` ("cool2hic: accept a complete cooler group whose format attribute is bytes"), library version 2.0.0. This is the commit HiCExplorer pinned (`HICX_HIC2COOL_COMMIT`); the repository HEAD at the time of vendoring (`91502d4`) is a rebased line with the same subject, renames the library to `hicxchange` and was not used.
- Changes: `include/hic2cool/version.hpp.in` became `version.hpp` with the project version written out. The command line tools (`apps/`), the Python module and package, the tests and test data, the examples and the CMake packaging were not copied.
- Licence: MIT, Copyright (c) 2017 Harvard Medical School and Copyright (c) 2026 Joachim Wolff. hicxchange is a C++ fork of hic2cool (https://github.com/4dn-dcic/hic2cool), the licence text below is the one of that project extended by the 2026 notice.

```
The MIT License

Copyright (c) 2017 Harvard Medical School
Copyright (c) 2026 Joachim Wolff (C++ rewrite, .hic version 9 support and cool2hic)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

## Arrowhead (Juicer tools)

- Location: `tools/arrowhead_impl.hpp`/`tools/arrowhead_impl.cpp`, used by `hicArrowhead`.
- Origin: the Arrowhead contact domain caller of Juicebox (https://github.com/aidenlab/Juicebox, `juicebox/tools/utils/juicer/arrowhead/`), ported from its Java source; the deviations are listed at the top of `arrowhead_impl.hpp`. Method: Rao et al., Cell 159 (2014), doi:10.1016/j.cell.2014.11.021.
- Licence: MIT, Copyright (c) 2016 Aiden Lab; the text is reproduced in the hicfilecpp section above.

## Stripenn

- Location: `tools/stripes_impl.hpp`/`tools/stripes_impl.cpp`, `tools/stripes_numpy.hpp`/`tools/stripes_numpy.cpp`, `tools/stripes_npy.hpp` and `tools/hicDetectStripes.cpp`.
- Origin: Stripenn 1.1.65.22 (https://github.com/VahediLab/stripenn), whose stripe detection stages are ported in the order of their floating point operations. Method: Yoon et al., Stripenn detects architectural stripes from chromatin conformation data using computer vision, Nature Communications (2022). The scikit-image and OpenCV routines Stripenn calls are implemented from their documented behaviour; no code of those two projects is included.
- Licence: MIT. The licence file of the repository carries the copyright line reproduced below unchanged.

```
Copyright (c) 2020 The Python Packaging Authority

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## pandas

- Location: `core/src/text_table.cpp` (the float parser `precise_xstrtod` of `read_csv`, its NA spellings and column type inference) and `core/src/qc_html.cpp` (the HTML of `Styler.to_html` used by `hicPrepareQCreport`).
- Origin: pandas 2.2.3 (https://github.com/pandas-dev/pandas), the version of the HiCExplorer 3.7.6 reference environment. These parts are ported so that tables are read and written with the same bits as in HiCExplorer 3.7.6.
- Licence: BSD 3-Clause.

```
BSD 3-Clause License

Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team
All rights reserved.

Copyright (c) 2011-2023, Open source contributors.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Unidecode

- Location: `core/src/unidecode_table.cpp`, generated by `scripts/generate_unidecode_table.py`; used by `hicx::unidecode` (`core/src/utilities.cpp`), which replaces non-ASCII characters in names as HiCExplorer 3.7.6 did.
- Origin: the transliteration tables of Unidecode 1.3.8 (https://github.com/avian2/unidecode). Original character transliteration tables: Copyright 2001, Sean M. Burke. Python code and later additions: Copyright 2024, Tomaž Šolc.
- Licence: GNU General Public License version 2 or (at the option of the licensee) any later version. HiCExplorer distributes the tables under version 3 of that licence, the licence of HiCExplorer (`../LICENSE`).

## HiCMatrix

- Location: `core/include/hicx/hic_matrix.hpp` and the matrix layer around it (`obsexp_ops.hpp`, `h5_file.hpp` and their sources).
- Origin: HiCMatrix 17.2 (https://github.com/deeptools/HiCMatrix), the matrix library of HiCExplorer 3, by Lucille Lopez-Delisle, Joachim Wolff, Leily Rabbani, Vivek Bhardwaj and Fidel Ramirez.
- Licence: GNU General Public License version 3, the licence of HiCExplorer (`../LICENSE`).

## Spectra

- Location: `third_party/spectra/include/Spectra/` (54 headers), the eigensolver of `hicPCA --eigenSolver lanczos`.
- Origin: https://github.com/yixuan/spectra, release v1.2.0 (commit `6841bcb`), copied unchanged from the conda-forge package spectra-cpp 1.2.0, which has no build for linux-aarch64.
- Licence: Mozilla Public License 2.0, Copyright (C) 2016-2025 Yixuan Qiu. The licence text is in `third_party/spectra/LICENSE`; every header carries the MPL-2.0 notice.
- SHA256 of the sorted per-file SHA256 list of the headers: `4fdc79096af4b182e72111dcc06b9671709d936c861dfa1e7c33c94095cd9c98`.
