// Author(s): Hector Joao Rivera Veruzco
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief Type recording the probabilistic equivalence reductions supported by the LTS library.
 * \author Hector Joao Rivera Verduzco
 */

#ifndef MCRL2_LTS_LTS_PROBABILISTIC_EQUIVALENCE_H
#define MCRL2_LTS_LTS_PROBABILISTIC_EQUIVALENCE_H

#include <string>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lts
{

/** \brief LTS equivalence relations.
 * \details This enumerated type defines equivalence relations on
 * LTSs. They can be used to reduce an LTS or decide whether two LTSs
 * are equivalent. */
enum lts_probabilistic_equivalence
{
  lts_probabilistic_eq_none,             /**< Unknown or no equivalence */
  lts_probabilistic_bisim_bem,          /**< Probabilistic bisimulation equivalence using the O(mn (log n + log m)) algorithm [Bier] */
  lts_probabilistic_bisim_grv           /**< Probabilistic bisimulation equivalence using the O(m(log n)) algorithm by Groote, Rivera Verduzco and de Vink */
};

/** \brief Determines the equivalence from a string.
 * \details The following strings may be used:
 * \li "none" for identity equivalence;
 * \li "pbisim" for Probabilistic bisimulation equivalence using the O(mn (log n + log m)) algorithm [Bier];
 *
 * \param[in] s The string specifying the equivalence.
 * \return The equivalence type specified by \a s.
 * If \a s is none of the above values then \a lts_probabilistic_eq_none is returned. */
inline
lts_probabilistic_equivalence parse_probabilistic_equivalence(const std::string& s)
{
  if (s == "none")
  {
    return lts_probabilistic_eq_none;
  }
  else if (s == "pbisim-bem")
  {
    return lts_probabilistic_bisim_bem;
  }
  else if (s == "pbisim")
  {
	  return lts_probabilistic_bisim_grv;
  }
  else
  {
    throw mcrl2::runtime_error("Unknown equivalence " + s + ".");
  }
}

// \overload
inline
std::istream& operator>>(std::istream& is, lts_probabilistic_equivalence& eq)
{
  try
  {
    std::string s;
    is >> s;
    eq = parse_probabilistic_equivalence(s);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/** \brief Gives the short name of an equivalence.
 * \param[in] eq The equivalence type.
 * \return A short string representing the equivalence specified by \a
 * eq. The returned value is one of the strings listed for
 * \ref parse_equivalence. */
inline std::string print_probabilistic_equivalence(const lts_probabilistic_equivalence& eq)
{
  switch (eq)
  {
    case lts_probabilistic_eq_none:
      return "none";
    case lts_probabilistic_bisim_bem:
      return "pbisim-bem";
    case lts_probabilistic_bisim_grv:
      return "pbisim";
    default:
      throw mcrl2::runtime_error("Unknown equivalence.");
  }
}

// \overload
inline
std::ostream& operator<<(std::ostream& os, const lts_probabilistic_equivalence& eq)
{
  os << print_probabilistic_equivalence(eq);
  return os;
}

/** \brief Gives a description of an equivalence.
 * \param[in] eq The equivalence type.
 * \return A string describing the equivalence specified by \a
 * eq. */
inline std::string description(const lts_probabilistic_equivalence& eq)
{
  switch(eq)
  {
    case lts_probabilistic_eq_none:
      return "identity equivalence";
    case lts_probabilistic_bisim_bem:
      return "probabilistic bisimulation equivalence using the O(mn (log n + log m)) algorithm by Baier, Engelen and Majster-Cederbaum, 2000";
	case lts_probabilistic_bisim_grv:
		return "probabilistic bisimulation equivalence using the O(m(log n)) algorithm by Groote, Rivera-Verduzco and de Vink, 2017";
    default:
      throw mcrl2::runtime_error("Unknown equivalence.");
  }
}


/** \brief LTS preorder relations.
 * \details This enumerated type defines equivalence relations on
 * LTSs. They can be used to reduce an LTS or decide whether two LTSs
 * are equivalent. */
enum lts_probabilistic_preorder
{
  lts_probabilistic_pre_none,             /**< Unknown or no preorder */
};

/** \brief Determines the probabilistic preorder from a string.
 * \details The following strings may be used:
 * \li "none" for identity preorder;
 *
 * \param[in] s The string specifying the equivalence.
 * \return The equivalence type specified by \a s.
 * If \a s is none of the above values then \a lts_probabilistic_pre_none is returned. */
inline
lts_probabilistic_preorder parse_probabilistic_preorder(std::string const& s)
{
  if (s == "none")
  {
    return lts_probabilistic_pre_none;
  }
  else
  {
    throw mcrl2::runtime_error("Unknown preorder " + s + ".");
  }
}

// \overload
inline
std::istream& operator>>(std::istream& is, lts_probabilistic_preorder& eq)
{
  try
  {
    std::string s;
    is >> s;
    eq = parse_probabilistic_preorder(s);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

/** \brief Gives the short name of an preorder.
 * \param[in] pre The preorder type.
 * \return A short string representing the preorder specified by \a
 * pre. The returned value is one of the strings listed for
 * \ref parse_preorder. */
inline std::string print_probabilistic_preorder(const lts_probabilistic_preorder& pre)
{
  switch (pre)
  {
    case lts_probabilistic_pre_none:
      return "none";
    default:
      throw mcrl2::runtime_error("Unknown preorder.");
  }
}

// \overload
inline
std::ostream& operator<<(std::ostream& os, const lts_probabilistic_preorder& pre)
{
  os << print_probabilistic_preorder(pre);
  return os;
}

/** \brief Gives a description of an preorder.
 * \param[in] pre The preorder type.
 * \return A string describing the preorder specified by \a
 * pre. */
inline std::string description(const lts_probabilistic_preorder& pre)
{
  switch (pre)
  {
    case lts_probabilistic_pre_none:
      return "identity preorder";
    default:
      throw mcrl2::runtime_error("Unknown preorder.");
  }
}

} // namespace mcrl2::lts


#endif // MCRL2_LTS_LTS_PROBABILISTIC_EQUIVALENCE_H



