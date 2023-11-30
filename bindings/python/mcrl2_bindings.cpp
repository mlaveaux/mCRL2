// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <pybind11/pybind11.h>

#include "mcrl2/process/parse.h"
#include "mcrl2/process/process_specification.h"
#include "mcrl2/lps/linearise.h"

namespace py = pybind11;

using namespace mcrl2::process;
using namespace mcrl2::lps;

struct LinearProcess {
  LinearProcess(const std::string& mcrl2_spec)
    : linear_spec(linearise(mcrl2_spec))
  {}

  stochastic_specification linear_spec;
};

PYBIND11_MODULE(mcrl2py, m) {
  py::class_<LinearProcess>(m, "LinearProcess")
      .def(py::init<const std::string &>());
}