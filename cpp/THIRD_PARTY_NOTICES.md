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
