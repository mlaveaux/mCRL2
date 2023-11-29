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
#include "mcrl2/utilities/logger.h"

namespace py = pybind11;

using namespace py::literals;
using namespace mcrl2::process;
using namespace mcrl2::lps;
using namespace mcrl2::log;

/// \brief The interface for a linear process specification.
struct LinearProcess {
  LinearProcess(const std::string& mcrl2_spec, t_lin_options options)
    : linear_spec(linearise(mcrl2_spec, options))
  {}

  stochastic_specification linear_spec;
};

struct LabelledTransitionSystem {
  LabelledTransitionSystem(const LinearProcess& lps) 
  {
    const auto& stochastic_lpsspec = lps.linear_spec;
  }
};

/// \brief File output class.
///
/// Ensure that prints are going to sys.stdout and sys.stderr
class python_output: public output_policy
{
  public:
    python_output()
    {
      //mcrl2::log::logger().register_output_policy
    }

    virtual ~python_output()
    {}

    virtual void output(const log_level_t level, const time_t timestamp, const std::string& msg, const bool print_time_information) override
    {
      py::print("{}", formatter::format(level, timestamp, msg, print_time_information));
    }
};

PYBIND11_MODULE(mcrl2py, m) {

  py::enum_<t_lin_method>(m, "Linearisation")
    .value("stack", lmStack)
    .value("regular", lmRegular)
    .value("regular2", lmRegular2)
    .export_values();

  py::class_<t_lin_options>(m, "LinearisationOptions")
    .def(py::init<>())
    .def_readwrite("method", &t_lin_options::lin_method)
    .def_readwrite("no_intermediate_cluster", &t_lin_options::no_intermediate_cluster)
    .def_readwrite("final_cluster",  &t_lin_options::final_cluster)
    .def_readwrite("newstate",  &t_lin_options::newstate)
    .def_readwrite("binary", &t_lin_options::binary)
    .def_readwrite("statenames", &t_lin_options::statenames)
    .def_readwrite("norewrite", &t_lin_options::norewrite);

  py::class_<LinearProcess>(m, "LinearProcess")
      .def(py::init<const std::string &, t_lin_options>(), "spec"_a, "options"_a = t_lin_options());
}