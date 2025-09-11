// Author(s): Jan Friso Groote, Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pbesrewr.cpp

#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/pbes/detail/instantiate_global_variables.h"
#include "mcrl2/pbes/detail/ppg_rewriter.h"
#include "mcrl2/pbes/detail/ppg_traverser.h"
#include "mcrl2/pbes/normalize.h"
#include "mcrl2/pbes/pbes_input_tool.h"
#include "mcrl2/pbes/pbes_output_tool.h"
#include "mcrl2/pbes/pbes_rewriter_tool.h"
#include "mcrl2/pbes/pbes_rewriter_type.h"
#include "mcrl2/pbes/rewrite.h"
#include "mcrl2/pbes/rewriters/bqnf_rewriter.h"
#include "mcrl2/pbes/rewriters/enumerate_quantifiers_rewriter.h"
#include "mcrl2/pbes/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/pbes/rewriters/pfnf_rewriter.h"
#include "mcrl2/pbes/rewriters/quantifiers_inside_rewriter.h"
#include "mcrl2/pbes/rewriters/simplify_quantifiers_rewriter.h"
#include "mcrl2/pbes/srf_pbes.h"
#include "mcrl2/utilities/input_output_tool.h"

using namespace mcrl2;
using namespace mcrl2::log;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::utilities;
using pbes_system::tools::pbes_input_tool;
using pbes_system::tools::pbes_output_tool;
using pbes_system::tools::pbes_rewriter_tool;
using data::tools::rewriter_tool;
using utilities::tools::input_output_tool;

class pbes_rewriter : public pbes_input_tool<pbes_output_tool<pbes_rewriter_tool<rewriter_tool<input_output_tool> > > >
{
  protected:
    using super = pbes_input_tool<pbes_output_tool<pbes_rewriter_tool<rewriter_tool<input_output_tool>>>>;

    /// \brief Returns the types of rewriters that are available for this tool.
    std::set<pbes_system::pbes_rewriter_type> available_rewriters() const override
    {
      std::set<pbes_system::pbes_rewriter_type> result = super::available_rewriters();
      result.insert(pbes_system::pbes_rewriter_type::quantifier_one_point);
      return result;
    }

  public:
    pbes_rewriter()
      : super(
        "pbesrewr",
        "Jan Friso Groote and Wieger Wesselink",
        "rewrite and simplify a PBES",
        "Rewrite the PBES in INFILE, remove quantified variables and write the resulting PBES to OUTFILE. "
        "If INFILE is not present, stdin is used. If OUTFILE is not present, stdout is used."
      )
    {}

    bool run() override
    {
      using namespace pbes_system;
      using namespace utilities;

      mCRL2log(verbose) << "pbesrewr parameters:" << std::endl;
      mCRL2log(verbose) << "  input file:         " << m_input_filename << std::endl;
      mCRL2log(verbose) << "  output file:        " << m_output_filename << std::endl;
      mCRL2log(verbose) << "  pbes rewriter:      " << m_pbes_rewriter_type << std::endl;

      // load the pbes
      pbes p;
      load_pbes(p, input_filename(), m_pbes_input_format);

      // data rewriter
      data::rewriter datar(p.data(), rewrite_strategy());

      // pbes rewriter
      switch (rewriter_type())
      {
        case pbes_rewriter_type::simplify:
        {
          simplify_quantifiers_data_rewriter<data::rewriter> pbesr(datar);
          //simplify_data_rewriter<data::rewriter> pbesr(datar);
          pbes_rewrite(p, pbesr);
          break;
        }
        case pbes_rewriter_type::quantifier_all:
        {
          bool enumerate_infinite_sorts = true;
          enumerate_quantifiers_rewriter pbesr(datar, p.data(), enumerate_infinite_sorts);
          pbes_rewrite(p, pbesr);
          break;
        }
        case pbes_rewriter_type::quantifier_finite:
        {
          bool enumerate_infinite_sorts = false;
          enumerate_quantifiers_rewriter pbesr(datar, p.data(), enumerate_infinite_sorts);
          pbes_rewrite(p, pbesr);
          break;
        }
        case pbes_rewriter_type::quantifier_inside:
        {
          quantifiers_inside_rewriter pbesr;
          pbes_rewrite(p, pbesr);
          break;
        }
        case pbes_rewriter_type::quantifier_one_point:
        {
          // apply the one point rule rewriter
          one_point_rule_rewriter pbesr;
          bool innermost = false;
          replace_pbes_expressions(p, pbesr, innermost); // use replace, since the one point rule rewriter does the recursion itself

          // post processing: apply the simplifying rewriter
          simplify_data_rewriter<data::rewriter> simp(datar);
          pbes_rewrite(p, simp);
          break;
        }
        case pbes_rewriter_type::pfnf:
        {
          pfnf_rewriter pbesr;
          pbes_system::normalize(p);
          pbes_rewrite(p, pbesr);
          break;
        }
        case pbes_rewriter_type::ppg:
        {
          //bool bqnf = detail::is_bqnf(p);
          //std::clog << "bqnf_traverser says: p is " << (bqnf ? "" : "NOT ") << "in BQNF." << std::endl;
          bool ppg = pbes_system::detail::is_ppg(p);
          if (ppg)
          {
            mCRL2log(log::verbose) << "PBES is already a PPG." << std::endl;
          }
          else
          {
            mCRL2log(log::verbose) << "Rewriting..." << std::endl;
            pbes q = pbes_system::detail::to_ppg(p);
            mCRL2log(log::verbose) << "Rewriting done." << std::endl;
            ppg = pbes_system::detail::is_ppg(q);
            if (!ppg)
            {
              throw(std::runtime_error("The result PBES if not a PPG!"));
            }
            p = q;
          }
          break;
        }
        case pbes_rewriter_type::srf:
        {
          pbes_system::detail::instantiate_global_variables(p);
          pbes_system::algorithms::normalize(p);
          auto result = pbes2srf(p,true);
          save_pbes(result.to_pbes(), output_filename(), m_pbes_output_format);  
          return true;
          break;
        }
        case pbes_rewriter_type::pre_srf:
        {
          pbes_system::detail::instantiate_global_variables(p);
          pbes_system::algorithms::normalize(p);
          auto result = pbes2pre_srf(p);
          save_pbes(result.to_pbes(), output_filename(), m_pbes_output_format);  
          return true;
          break;    
        }
        case pbes_rewriter_type::bqnf_quantifier:
        {
          bqnf_rewriter pbesr;
          pbes_rewrite(p, pbesr);
          break;
        }
      }
      
      save_pbes(p, output_filename(), m_pbes_output_format);
      return true;
    }

};

int main(int argc, char* argv[])
{
  return pbes_rewriter().execute(argc, argv);
}
