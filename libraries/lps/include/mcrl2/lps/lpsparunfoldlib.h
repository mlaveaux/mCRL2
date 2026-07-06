// Author(s): Frank Stappers, Thomas Neele
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/lpsparunfoldlib.h
/// \brief Unfold process parameters in mCRL2 process specifications.

#ifndef MCRL2_LPS_LPSPARUNFOLDLIB_H
//Fileinfo
#define MCRL2_LPS_LPSPARUNFOLDLIB_H

#include "mcrl2/data/consistency.h"
#include "mcrl2/data/representative_generator.h"
#include "mcrl2/data/unfold_pattern_matching.h"
#include "mcrl2/lps/detail/lps_algorithm.h"
#include "mcrl2/lps/replace_capture_avoiding.h"

namespace mcrl2::lps
{

/// \brief Element in the cache that keeps track of the information for a single
///        unfolded sort, say s.
struct unfold_cache_element
{
  mcrl2::data::function_symbol_vector affected_constructors; // constructors for s
  mcrl2::data::basic_sort fresh_basic_sort; // new sort introduced to represent s
  mcrl2::data::function_symbol_vector new_constructors; // constructors for fresh_basic_sort
  mcrl2::core::identifier_string case_function_name; // name of the case function introduced to deal with sort s
  // case functions for sort s. Note this may be more due to, e.g., adding distribution rules.
  // invariant: case_functions.front() is the case function
  // C: fresh_basic_sort # s # ... # s -> s
  std::map<mcrl2::data::sort_expression, mcrl2::data::function_symbol> case_functions;
  mcrl2::data::function_symbol determine_function;  // Det function for s
  mcrl2::data::function_symbol_vector projection_functions; // pi functions for s
};

namespace detail
{
  /// \brief Fresh variable generator for the arguments of a function symbol.
  ///
  /// Intended for use to generate variables in the left hand side of a data equation.
  /// Therefore, reuse of fresh variables is allowed in multiple equations.
  /// This generator reuses such variable names.
  class data_equation_argument_generator
  {
    public:
      data_equation_argument_generator(
        data::set_identifier_generator& identifier_generator)
        : m_identifier_generator(identifier_generator)
      {}

      /// Generate argument variables for f.
      data::variable_vector arguments(const data::function_symbol& f)
      {
        data::variable_vector result;
        std::map<data::sort_expression, std::size_t> used_arguments;
        if (is_function_sort(f.sort()))
        {
          for(const data::sort_expression& arg_sort: atermpp::down_cast<data::function_sort>(f.sort()).domain())
          {
            // Check if we already used the sort for this function symbol.
            if(used_arguments.find(arg_sort) == used_arguments.end())
            {
              used_arguments[arg_sort] = 0;
            }

            // Check if we need to generate a new variable.
            if(m_variables[arg_sort].size() == used_arguments[arg_sort])
            {
              m_variables[arg_sort].emplace_back(m_identifier_generator("d"), arg_sort);
            }

            // Push new variable
            result.push_back(m_variables[arg_sort].at(used_arguments[arg_sort]));
            ++used_arguments[arg_sort];
          }
        }
        return result;
      }

    protected:
      data::set_identifier_generator& m_identifier_generator;
      std::map<data::sort_expression, data::variable_vector> m_variables;
  };

  class unfold_data_manager
  {
    private:
    /// \brief cache for previously unfolded sorts.
    /// facilitates reuse of previously introduced sorts and function symbols.
    std::map< mcrl2::data::sort_expression , unfold_cache_element >& m_cache;

    data::data_specification& m_dataspec;

    /// \brief set of identifiers to use during fresh variable generation
    mcrl2::data::set_identifier_generator m_identifier_generator;

    /// \brief generator for arguments in left hand side of data equations
    detail::data_equation_argument_generator m_data_equation_argument_generator;

    /// \brief a generator for default data expressions of a given sort;
    mcrl2::data::representative_generator m_representative_generator;

    /// \brief Boolean indicating whether rewrite rules may be added that could make
    ///        the data specification inconsistent.
    bool m_possibly_inconsistent;

    public:
    unfold_data_manager(std::map< mcrl2::data::sort_expression , unfold_cache_element >& cache, data::data_specification& dataspec, bool possibly_inconsistent)
    : m_cache(cache)
    , m_dataspec(dataspec)
    , m_identifier_generator()
    , m_data_equation_argument_generator(m_identifier_generator)
    , m_representative_generator(dataspec)
    , m_possibly_inconsistent(possibly_inconsistent)
    {}

    const data::data_specification& dataspec() const { return m_dataspec; }

    void add_used_identifier(const core::identifier_string& id)
    {
      m_identifier_generator.add_identifier(id);
    }

    void add_used_identifiers(const std::set<core::identifier_string>& ids)
    {
      m_identifier_generator.add_identifiers(ids);
    }

    data::set_identifier_generator& id_gen()
    {
      return m_identifier_generator;
    }

    unfold_cache_element& get_cache_element(const data::sort_expression& sort);

    bool is_cached(const data::sort_expression& sort) const
    {
      return m_cache.find(sort) != m_cache.end();
    }

    bool is_constructor(const data::function_symbol& f) const
    {
      return utilities::detail::contains(m_dataspec.constructors(), f);
    }

    const std::vector<data::function_symbol>& get_constructors(const data::sort_expression& sort)
    {
      return get_cache_element(sort).affected_constructors;
    }

    data::data_expression create_cases(const data::data_expression& target, const data::data_expression_vector& rhss)
    {
      unfold_cache_element& cache_elem = get_cache_element(target.sort());
      data::application first_arg(cache_elem.determine_function, target);

      data::data_expression_vector args;
      args.push_back(first_arg);
      for (const data::data_expression& rhs: rhss)
      {
        args.push_back(rhs);
      }

      // Create the case function on demand:
      // the output sort may not have been creation of case functions.
      const data::function_symbol case_func = create_case_function(target.sort(), rhss.front().sort());
      return data::application(case_func, args);
    }

    const data::function_symbol_vector& get_projection_funcs(const data::function_symbol& f)
    {
      return get_cache_element(f.sort().target_sort()).projection_functions;
    }

    /** \brief  Generates a fresh basic sort given a sort expression.
      * \param  sort This sort's name will be used to derive a name for the new sort
      * \return A fresh basic sort.
    **/
    mcrl2::data::basic_sort generate_fresh_basic_sort(const data::sort_expression& sort);

    /** \brief  Generates a fresh name for a constructor or mapping.
      * \param  str a string value. The value is used to generate a fresh
      *         name for a constructor or mapping.
      * \post   A fresh name for a constructor or mapping is generated.
      * \return A fresh name for a constructor or mapping.
    **/
    mcrl2::core::identifier_string generate_fresh_function_symbol_name(const std::string& str);

    /** \brief Generates variable of type sort based on a given string str.
      * \param str a string value. The value is used to generate a fresh
      *         variable name.
      * \param sort The sort of the variable to generate.
      * \post   A fresh variable is generated, which has an unique name.
      * \return A fresh variable.
    **/
    mcrl2::data::variable generate_fresh_variable(std::string str, const data::sort_expression& sort);


    /** \brief  Creates the determine function.
      * \param sort The sort on which this function operates
      * \return A function that maps constructors to the fresh basic sort
    **/
    void create_determine_function(const data::sort_expression& sort);

    /** \brief  Creates projection functions for the unfolded process parameter.
      * \param sort The sort on which these functions operate
      * \return A function that returns the projection functions for the
      *         constructors of the unfolded process parameter.
    **/
    void create_projection_functions(const data::sort_expression& sort);

    /** \brief  Determines the constructors that are affected with the unfold
      *         process parameter.
      * \param sort The sort for which to find constructors
      * \post The constructors that are affected with the unfold
      *         process parameter are stored in m_affected_constructors
    **/
    void determine_affected_constructors(const data::sort_expression& sort);

    /** \brief  Creates a set of constructors for the fresh basic sort
      * \param sort The sort for which to create analogous constructors
      * \return The constructors that are created for the fresh basic sort
    **/
    void create_new_constructors(const data::sort_expression& sort);

    /** \brief  Creates the case function with number of arguments determined by
     *          the number of affected constructors, the sort of the arguments and
     *          result are determined by sort..
     *  \param  det_sort The sort whose constructor is determined in the first argument
     *  \param  output_sort The sort of the arguments and return sort of the case function
      * \return A function that returns the corresponding constructor given the
      *         case selector and constructors.
    **/
    data::function_symbol create_case_function(const data::sort_expression& det_sort, const data::sort_expression& output_sort);

    /** \brief Create distribution rules for distribution_functions over case_functions
    **/
    void create_distribution_law_over_case(const data::sort_expression& sort,
      const data::function_symbol& function_for_distribution,
      data::function_symbol case_function);

    /** \brief Create the data equations for case functions */
    void generate_case_function_equations(
      const data::sort_expression& sort,
      const data::function_symbol& case_function);

    /** \brief Create the data equations for the determine function */
    void generate_determine_function_equations(const data::sort_expression& sort);

    /** \brief Create the data equations for the projection functions */
    void generate_projection_function_equations(const data::sort_expression& sort);

    static bool char_filter(char c)
    {
      // Put unwanted characters here
      return c==' ' || c==':' || c==',' || c=='|'
             || c=='>' || c=='[' || c==']' || c=='@'
             || c=='.' || c=='{' || c=='}' || c=='#'
             || c=='%' || c=='&' || c=='*' || c=='!'
             || c=='(' || c==')'
             ;
    }

    std::string filter_illegal_characters(std::string in) const
    {
      in.resize(std::remove_if(in.begin(), in.end(), &char_filter) - in.begin());
      return in;
    }
  };

  /// @brief Class for unfolding expressions f(a1,...,an) based on the pattern-matching
  /// rewrite rules that define f.
  class pattern_match_unfolder
  {
    private:
    unfold_data_manager& m_datamgr;
    data::representative_generator m_repgen;
    std::map<data::function_symbol, bool> m_is_pattern_matching;
    std::map<data::function_symbol, data::data_equation> m_new_eqns;

    /// @brief Finds all rewriting equations for f
    data::data_equation_vector find_equations(const data::function_symbol& f)
    {
      data::data_equation_vector result;
      for (const data::data_equation& eqn: m_datamgr.dataspec().equations())
      {
        if (data::detail::get_top_fs(eqn.lhs()) == f)
        {
          result.push_back(eqn);
        }
      }
      return result;
    }

    /// @brief Checks whether f is defined by pattern matching
    bool is_pattern_matching(const data::function_symbol& f)
    {
      auto find_result = m_is_pattern_matching.find(f);
      if (find_result != m_is_pattern_matching.end())
      {
        return find_result->second;
      }
      data::data_equation_vector eqns = find_equations(f);
      bool result = std::all_of(eqns.begin(), eqns.end(),
        [&](auto& eqn) { return data::is_pattern_matching_rule(m_datamgr, eqn); });
      m_is_pattern_matching.insert_or_assign(f, result);
      return result;
    }

    data::data_expression unfolded_expr(const data::function_symbol& f, const data::data_expression_vector& args)
    {
      data::data_equation new_eqn;
      auto find_result = m_new_eqns.find(f);
      if (find_result != m_new_eqns.end())
      {
        new_eqn = find_result->second;
      }
      else
      {
        new_eqn = data::unfold_pattern_matching(f, find_equations(f), m_datamgr, m_repgen, m_datamgr.id_gen());
        m_new_eqns.emplace(f, new_eqn);
      }

      std::map<data::variable, data::data_expression> sigma;
      auto it1 = atermpp::down_cast<data::application>(new_eqn.lhs()).begin();
      auto it2 = args.begin();
      while (it2 != args.end())
      {
        sigma[atermpp::down_cast<data::variable>(*it1++)] = *it2++;
      }
      data::map_substitution sigma1(sigma);


      return data::replace_all_variables(new_eqn.rhs(), sigma1);
    }

    public:
    pattern_match_unfolder(unfold_data_manager& datamgr)
    : m_datamgr(datamgr)
    , m_repgen(datamgr.dataspec())
    {}

    bool is_constructor(const data::function_symbol& f)
    {
      return m_datamgr.is_constructor(f);
    }

    /// @brief Checks whether expr is of the shape Det(h(a1,...,an)) or pi(h(a1,...,an)),
    /// where h is defined by patter matching, and, if so, unfolds h(a1,...,an).
    bool is_det_or_pi(const data::application& expr) const
    {
      using utilities::detail::contains;

      const data::function_symbol f = data::detail::get_top_fs(expr);
      // If f is not unary, then it is certainly unequal to Det or pi
      if (f == data::function_symbol() || expr.size() != 1)
      {
        return false;
      }
      const data::sort_expression& arg_sort = *atermpp::down_cast<data::function_sort>(f.sort()).domain().begin();
      // The argument of f must be of an unfolded sort
      if (!m_datamgr.is_cached(arg_sort))
      {
        return false;
      }
      const unfold_cache_element& cache_elem = m_datamgr.get_cache_element(arg_sort);
      if (f != cache_elem.determine_function && !contains(cache_elem.projection_functions, f))
      {
        return false;
      }

      return true;
    }

    bool can_unfold(const data::data_expression& x)
    {
      if (!data::is_application(x))
      {
        return false;
      }

      const data::function_symbol f = data::detail::get_top_fs(atermpp::down_cast<data::application>(x));
      if (f == data::function_symbol() || !is_pattern_matching(f))
      {
        return false;
      }
      auto udm = m_datamgr.dataspec().mappings();
      const data::data_specification& dataspec = m_datamgr.dataspec();
      if (std::find_if(udm.begin(), udm.end(),
          [&](const auto& f2){ return f.name() == f2.name() && dataspec.equal_sorts(f.sort(), f2.sort()); }) == udm.end())
      {
        // f is not a mapping, but likely a constructor
        return false;
      }
      if (!matches_only_known_sorts(f))
      {
        return false;
      }
      return true;
    }

    std::vector<std::size_t> pattern_matching_args(const data::function_symbol& f)
    {
      std::vector<std::size_t> result;
      if (!data::is_function_sort(f.sort()))
      {
        return result;
      }
      data::data_equation_vector eqns = find_equations(f);

      std::size_t f_num_args = atermpp::down_cast<data::function_sort>(f.sort()).domain().size();
      for (std::size_t arg = 0; arg < f_num_args; arg++)
      {
        for (const data::data_equation& eq: eqns)
        {
          const data::application& lhs_appl = atermpp::down_cast<data::application>(eq.lhs());
          if (!data::is_variable(lhs_appl[arg]))
          {
            result.push_back(arg);
            break;
          }
        }
      }
      return result;
    }

    /// @brief Determines whether f pattern matches on argument arg
    /// @pre this->can_unfold(f)
    bool matches_only_known_sorts(const data::function_symbol& f)
    {
      if (!data::is_function_sort(f.sort()))
      {
        return true;
      }

      const data::sort_expression_list& domain = atermpp::down_cast<data::function_sort>(f.sort()).domain();
      const data::sort_expression_vector domain_vec{domain.begin(), domain.end()};
      for (const std::size_t i: pattern_matching_args(f))
      {
        if (!m_datamgr.is_cached(domain_vec[i]))
        {
          return false;
        }
      }
      return true;
    }

    /// @brief Unfolds expr if it is of the shape h(a1,...,an) and h is defined by pattern matching
    /// @pre this->can_unfold(x)
    template <class T>
    void operator()(T& result, const data::application& x)
    {
      const data::data_expression_vector args(x.begin(), x.end());
      result = unfolded_expr(data::detail::get_top_fs(x), args);
    }
  };


  template <template <class> class Builder>
  struct replace_pattern_match_builder: public Builder<replace_pattern_match_builder<Builder> >
  {
    using super = Builder<replace_pattern_match_builder<Builder>>;
    using super::enter;
    using super::leave;
    using super::apply;
    using super::update;

    pattern_match_unfolder& m_unfolder;
    bool m_currently_recursing = false;
    std::size_t m_current_depth = 0;

    /// \brief Maximum number of times an expression is unfolded before recursion
    ///        stops; unfolding this often is enough to rewrite Det() and pi().
    static constexpr std::size_t max_unfold_depth = 3;

    replace_pattern_match_builder(pattern_match_unfolder& unfolder)
      : m_unfolder(unfolder)
    {}

    bool is_applied_to_constructor(const data::application& x)
    {
      using utilities::detail::contains;
      const std::vector<std::size_t> pattern_matching_args = m_unfolder.pattern_matching_args(data::detail::get_top_fs(x));
      for(const std::size_t i: pattern_matching_args)
      {
        if (m_unfolder.is_constructor(data::detail::get_top_fs(x[i])))
        {
          return true;
        }
      }
      return false;
    }

    template <class T>
    void apply(T& result, const data::application& x)
    {
      if (m_currently_recursing)
      {
        if (m_current_depth >= max_unfold_depth || m_unfolder.is_constructor(data::detail::get_top_fs(x)))
        {
          // Stop recursing after unfolding three times or when meeting a constructor
          // In the latter case, we have done enough to rewrite Det() and pi()
          result = x;
          return;
        }
        else if (data::is_if_application(x))
        {
          // In if-statements, do not traverse the condition
          data::data_expression branch1;
          super::apply(branch1, x[1]);
          data::data_expression branch2;
          super::apply(branch2, x[2]);

          data::make_application(result,
            data::if_(x.sort()),
            x[0],
            branch1,
            branch2
          );
        }
        else if (m_unfolder.can_unfold(x) && is_applied_to_constructor(x))
        {
          data::data_expression intermediate_result;
          m_unfolder(intermediate_result, x);
          // Recursively apply unfolding
          m_current_depth++;
          super::apply(result, intermediate_result);
          m_current_depth--;
        }
        else
        {
          super::apply(result, x);
        }
      }
      else
      {
        // Determine whether we see Det(f(..)) or pi(f(..))
        // If so, unfold and start recursing
        if (m_unfolder.is_det_or_pi(x) && m_unfolder.can_unfold(x[0]))
        {
          data::data_expression intermediate_result1;
          data::data_expression intermediate_result2;
          m_unfolder(intermediate_result1, atermpp::down_cast<data::application>(x[0]));

          m_currently_recursing = true;
          m_current_depth++;
          super::apply(intermediate_result2, intermediate_result1);
          m_current_depth--;

          mCRL2log(log::debug) << "Unfolded " << x[0] << " into " << intermediate_result2 << std::endl;
          data::make_application(result, x.head(), intermediate_result2);
          m_currently_recursing = false;
        }
        else
        {
          super::apply(result, x);
        }
      }
    }
  };

  inline
  data::data_expression unfold_pattern_matching(const data::data_expression& x,
                               pattern_match_unfolder& unfolder
                              )
  {
    data::data_expression result;
    replace_pattern_match_builder<data::data_expression_builder>(unfolder).apply(result, x);
    return result;
  }
} // namespace detail

class lpsparunfold: public detail::lps_algorithm<lps::stochastic_specification>
{
  using super = typename detail::lps_algorithm<lps::stochastic_specification>;

public:

  // old parameter, case function, determinizing parameter, replacement expressions
  using case_func_replacement = std::tuple<data::variable,
      std::map<data::sort_expression, data::function_symbol>,
      data::variable,
      data::data_expression_vector>;

  /** \brief  Constructor for lpsparunfold algorithm.
   * \param[in] spec which is a valid mCRL2 process specification.
   * \param[in,out] cache Cache to store information for reuse.
   * \param[in] alt_case_placement If true, case functions are placed at a higher level.
   * \param[in] possibly_inconsistent If true, case functions over Booleans are replaced by a disjunction of
   *conjunctions. For this to be correct, the unfolded sort needs to satisfy some restrictions.
   * \post   The content of mCRL2 process specification analysed for useful information and class variables are set.
   **/
  lpsparunfold(lps::stochastic_specification& spec,
      std::map<data::sort_expression, unfold_cache_element>& cache,
      bool alt_case_placement = false,
      bool possibly_inconsistent = false,
      bool unfold_pattern_matching = true);

  /** \brief  Applies lpsparunfold algorithm on a process parameter of an mCRL2 process specification .
   *  \pre algorithm has not been called before.
   *  \param[in] parameter_at_index An integer value that represents the index value of an process parameter.
   *  \post   The process parameter at index parameter_at_index is unfolded in the mCRL2 process specification.
   **/
  void algorithm(std::size_t parameter_at_index);

private:
  /// \brief set to true when the algorithm has been run once; as the algorithm should
  /// run only once...
  bool m_run_before = false;

  /// @brief Bookkeeper for recogniser and projection functions.
  detail::unfold_data_manager m_datamgr;
  detail::pattern_match_unfolder m_pattern_unfolder;

  /// \brief The process parameter that needs to be unfold.
  mcrl2::data::variable m_unfold_parameter;

  /// \brief The process parameters that are inserted.
  mcrl2::data::variable_vector m_injected_parameters;

  /// \brief Boolean to indicate if alternative placement of case functions should be used.
  bool m_alt_case_placement;

  /// \brief Indicates whether functions defined by pattern matching that occur in the scope of
  /// a Det or pi in a state update should be unfolded.
  bool m_unfold_pattern_matching;

  // data::data_expression apply_case_function(const data::data_expression& expr, const case_func_replacement&
  // case_funcs);
  case_func_replacement parameter_case_function();

  /** \brief  Get the process parameter at given index
   * \param  index The index of the parameter which must be obtained.
   * \return the process parameter at given index.
   **/
  mcrl2::data::variable process_parameter_at(std::size_t index);

  /** \brief  substitute function for replacing process parameters with unfolded process parameters functions.
   * \return substitute function for replacing process parameters with unfolded process parameters functions.
   **/
  std::map<mcrl2::data::variable, mcrl2::data::data_expression> parameter_substitution();

  data::data_expression apply_function(const data::function_symbol& f, const data::data_expression& de) const;

  /** \brief unfolds a data expression into a vector of process parameters
   * \param  de the data expression
   * \return The following vector: < det(de), pi_0(de), ... ,pi_n(de) >
   **/
  mcrl2::data::data_expression_vector unfold_constructor(const mcrl2::data::data_expression& de);

  /** \brief substitute unfold process parameter in the linear process
   * \param  parameter_at_index the parameter index
   * \post the process parameter at given index is unfolded in the the linear process
   **/
  void update_linear_process(std::size_t parameter_at_index);

  /** \brief substitute unfold process parameter in the initialization of the linear process
   * \param  parameter_at_index the parameter index
   * \post the initialization for the linear process is updated by unfolding the parameter at given index is unfolded
   **/
  void update_linear_process_initialization(std::size_t parameter_at_index);

  // Applies 'process unfolding' to a sequence of summands.
  void unfold_summands(mcrl2::lps::stochastic_action_summand_vector& summands);
};


template <template <class> class Builder, template <template <class> class, class, class> class Binder>
struct parunfold_replacement: public
Binder<Builder, parunfold_replacement<Builder, Binder>, parunfold_replacement<Builder, Binder>>
{
  using super = Binder<Builder, parunfold_replacement<Builder, Binder>, parunfold_replacement<Builder, Binder>>;
  using super::enter;
  using super::leave;
  using super::apply;
  using super::update;

  data::detail::capture_avoiding_substitution_updater<parunfold_replacement<Builder, Binder>> sigma1;
  lpsparunfold::case_func_replacement case_funcs;
  data::data_expression current_replacement;

  parunfold_replacement(const lpsparunfold::case_func_replacement& case_funcs,
                        data::set_identifier_generator& id_generator)
  : super(sigma1)
  , sigma1(*this, id_generator)
  {
    this->case_funcs = case_funcs;
  }

  template <class T>
  void apply(T& result, const data::application& x)
  {
    if (current_replacement != data::data_expression() || data::is_and(x) || data::is_or(x) || data::is_not(x) || data::is_imp(x) || data::is_if_application(x))
    {
      // if no placement of case functions is underway, or we are still traversing the regular boolean operators, we continue as usual
      super::apply(result, x);
    }
    else
    {
      // place the case functions here
      apply_case_function(result, x);
    }
  }

  void apply_case_function(data::data_expression& result, const data::application& expr)
  {
    auto& [par, case_f, det_f, replacements] = case_funcs;

    if (data::find_free_variables(expr).count(par) == 0)
    {
      // variable to be replaced does not occur here
      // make sure to still apply the substitutions necessary for the capture avoiding tricks
      // NB: stack overflow happens if type of second argument is 'data::data_expression'.
      super::apply(result, expr);
    }
    else
    {
      data::data_expression_vector args;
      args.push_back(det_f);

      for (const data::data_expression& r : replacements)
      {
        current_replacement = r;
        data::data_expression arg;
        super::apply(arg, expr);
        args.push_back(arg);
      }
      current_replacement = data::data_expression();

      if (case_f.find(expr.sort()) == case_f.end())
      {
        throw mcrl2::runtime_error("Case function with target sort " +
                                   data::pp(expr.sort()) + " not declared.");
      }
      result = data::application(case_f[expr.sort()], args);
    }
  }

  // Substitution application
  data::data_expression operator()(const data::variable& x)
  {
    if (current_replacement == data::data_expression())
    {
      return x;
    }
    if (std::get<0>(case_funcs) != x)
    {
      return x;
    }
    return current_replacement;
  }
};

template <template <class> class Builder, template <template <class> class, class, class> class Binder>
parunfold_replacement<Builder, Binder>
apply_parunfold_replacement_builder(const lpsparunfold::case_func_replacement& case_funcs,
                                    data::set_identifier_generator& id_generator)
{
  return parunfold_replacement<Builder, Binder>(case_funcs, id_generator);
}

template <typename T>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void insert_case_functions(T& x,
    const lpsparunfold::case_func_replacement& cfv,
    data::set_identifier_generator& id_generator)
{
  apply_parunfold_replacement_builder<lps::data_expression_builder, lps::detail::add_capture_avoiding_replacement>(cfv, id_generator).update(x);
}

template <typename T>
  requires(!std::is_base_of_v<atermpp::aterm, T>)
void insert_case_functions(T& x,
    const lpsparunfold::case_func_replacement& cfv
  )
{
  data::set_identifier_generator id_generator;
  id_generator.add_identifiers(lps::find_identifiers(x));
  auto& [par, case_f, det_f, replacements] = cfv;

  // assumes all case functions have the same name and that par.sort() is in the map of case functions
  id_generator.add_identifier(case_f.find(par.sort())->second.name());
  id_generator.add_identifier(det_f.name());
  for (const data::data_expression& r: replacements)
  {
    id_generator.add_identifiers(data::find_identifiers(r));
  }

  insert_case_functions(x, cfv, id_generator);
}

} // namespace mcrl2::lps


#endif
