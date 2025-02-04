// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief Supperted preorders for LTSes
 * \author Jan Friso Groote, Bas Ploeger, Muck van Weerdenburg, Jeroen Keiren
 */

#ifndef MCRL2_LTS_LTS_PREORDER_H
#define MCRL2_LTS_LTS_PREORDER_H

#include <string>
#include "mcrl2/utilities/exception.h"

namespace mcrl2
{
namespace lts
{

/** \brief LTS preorder relations.
 * \details This enumerated type defines preorder relations on LTSs.
 * They can be used to decide whether one LTS is behaviourally
 * contained in another LTS. */
enum class lts_preorder
{
  lts_pre_none,   /**< Unknown or no preorder */
  lts_pre_sim,    /**< Strong simulation preorder */
  lts_pre_ready_sim,    /**< Strong ready simulation preorder */  
  lts_pre_trace,  /**< Strong trace preorder */
  lts_pre_weak_trace,   /**< Weak trace preorder */
  lts_pre_trace_anti_chain,  /**< Trace preorder based on anti chains */
  lts_pre_weak_trace_anti_chain, /**< Weak trace preorder based on anti chains */
  lts_pre_failures_refinement,    /**< Failures refinement based on anti chains */
  lts_pre_weak_failures_refinement, /**< Weak failures refinement based on anti chains */
  lts_pre_failures_divergence_refinement, /**< Failures divergence refinement based on anti chains, which is automatically weak */
  lts_preorder_min = lts_pre_none,
  lts_preorder_max = lts_pre_failures_divergence_refinement
};

/** \brief Determines the preorder from a string.
 * \details The following strings may be used:
 * \li "sim" for strong simulation preorder;
 * \li "trace" for strong trace preorder;
 * \li "weak-trace" for weak trace preorder.
 *
 * \param[in] s The string specifying the preorder.
 * \return The preorder type specified by \a s.
 * If \a s is none of the above values then \a lts_pre_none is returned. */
inline
lts_preorder parse_preorder(std::string const& s)
{
  if (s == "unknown")
  {
    return lts_preorder::lts_pre_none;
  }
  else if (s == "sim")
  {
    return lts_preorder::lts_pre_sim;
  }
  else if (s == "ready-sim")
  {
    return lts_preorder::lts_pre_ready_sim;
  }  
  else if (s == "trace")
  {
    return lts_preorder::lts_pre_trace;
  }
  else if (s == "weak-trace")
  {
    return lts_preorder::lts_pre_weak_trace;
  }
  else if (s == "trace-ac")
  {
    return lts_preorder::lts_pre_trace_anti_chain;
  }
  else if (s == "weak-trace-ac")
  {
    return lts_preorder::lts_pre_weak_trace_anti_chain;
  }
  else if (s == "failures")
  {
    return lts_preorder::lts_pre_failures_refinement;
  }
  else if (s == "weak-failures")
  {
    return lts_preorder::lts_pre_weak_failures_refinement;
  }
  else if (s == "failures-divergence")
  {
    return lts_preorder::lts_pre_failures_divergence_refinement;
  }
  else
  {
    throw mcrl2::runtime_error("unknown preorder " + s);
  }
}

// \overload
inline
std::istream& operator>>(std::istream& is, lts_preorder& eq)
{
  try
  {
    std::string s;
    is >> s;
    eq = parse_preorder(s);
  }
  catch (mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/** \brief Gives the short name of a preorder.
 * \param[in] pre The preorder type.
 * \return A short string representing the preorder specified by \a
 * pre. The returned value is one of the strings listed for
 * \ref parse_preorder. */
inline
std::string print_preorder(const lts_preorder pre)
{
  switch(pre)
  {
    case lts_preorder::lts_pre_none:
      return "unknown";
    case lts_preorder::lts_pre_sim:
      return "sim";
    case lts_preorder::lts_pre_ready_sim:
      return "ready-sim";      
    case lts_preorder::lts_pre_trace:
      return "trace";
    case lts_preorder::lts_pre_weak_trace:
      return "weak-trace";
    case lts_preorder::lts_pre_trace_anti_chain:
      return "trace-ac";
    case lts_preorder::lts_pre_weak_trace_anti_chain:
      return "weak-trace-ac";
    case lts_preorder::lts_pre_failures_refinement:
      return "failures";
    case lts_preorder::lts_pre_weak_failures_refinement:
      return "weak-failures";
    case lts_preorder::lts_pre_failures_divergence_refinement:
      return "failures-divergence";
    default:
      throw mcrl2::runtime_error("unknown preorder");
  }
}

inline
std::ostream& operator<<(std::ostream& os, const lts_preorder pre)
{
  os << print_preorder(pre);
  return os;
}

/** \brief Gives a description of a preorder.
 * \param[in] pre The preorder type.
 * \return A string describing the preorder specified by \a
 * eq. */
inline
std::string description(const lts_preorder pre)
{
  switch(pre)
  {
    case lts_preorder::lts_pre_none:
      return "default void preorder";
    case lts_preorder::lts_pre_sim:
      return "strong simulation preorder";
    case lts_preorder::lts_pre_ready_sim:
      return "strong ready simulation preorder";      
    case lts_preorder::lts_pre_trace:
      return "strong trace preorder";
    case lts_preorder::lts_pre_weak_trace:
      return "weak trace preorder";
    case lts_preorder::lts_pre_trace_anti_chain:
      return "trace preorder based on an anti chain algorithm";
    case lts_preorder::lts_pre_weak_trace_anti_chain:
      return "weak trace preorder based on an anti chain algorithm";
    case lts_preorder::lts_pre_failures_refinement:
      return "failures refinement";
    case lts_preorder::lts_pre_weak_failures_refinement:
      return "weak failures refinement";
    case lts_preorder::lts_pre_failures_divergence_refinement:
      return "failures divergence refinement (automatically weak)";
    default:
      throw mcrl2::runtime_error("unknown preorder");
  }
}

} // namespace lts
} // namespace mcrl2

#endif // MCRL2_LTS_LTS_PREORDER_H



