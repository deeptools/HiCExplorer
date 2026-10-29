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
- hiccpp re-implements behaviour of straw/hicstraw 1.3.1 (https://github.com/aidenlab/straw, MIT, Aiden Lab; `src/reader.cpp`), of Juicer tools and Juicebox 1.22.01 and 2.20.00 (https://github.com/aidenlab/Juicebox, MIT, Aiden Lab; `src/writer.cpp`, `src/norms.cpp`) and used hic2cool (https://github.com/4dn-dcic/hic2cool, MIT, Harvard Medical School) as a second reading of the version 8 layout. All are MIT licensed like hiccpp; no code is copied from them.

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

