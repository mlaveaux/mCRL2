// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file linearise.cpp
/// \brief This file contains code to transform an mCRL2 process
///        into a linear process.

/* This file contains the implementation of an mCRL2 lineariser.

   It is based on the implementation of the mCRL lineariser, on which work
   started on 12 juli 1997.  This lineariser was based on the CWI technical
   report "The Syntax and Semantics of Timed mCRL", by J.F. Groote.

   Everybody is free to use this software, provided it is not changed.

   In case problems are encountered when using this software, please report
   them to J.F. Groote, TU/e, Eindhoven, J.F.Groote@tue.nl

   This software comes as it is. I.e. the author assumes no responsibility for
   the use of this software.
*/

// The following define allows collecting and printing statistics about the operations performed in the linearizer.
// Currently, we collect and print information about:
// the operations performed on processes,
// the number of summands and the total number of actions in the multications in these summands.
// Note that all logs are written to std::cout, so if this flag is enabled, tool output can not be streamed through std::cout
//#define MCRL2_LOG_LPS_LINEARISE_STATISTICS 1

//mCRL2 data
#include "mcrl2/atermpp/aterm.h"
#include "mcrl2/data/substitutions/maintain_variables_in_rhs.h"
#include "mcrl2/data/fourier_motzkin.h"
#include "mcrl2/data/enumerator.h"
#include "mcrl2/data/real_utilities.h"

// linear process libraries.
#include "mcrl2/lps/constelm.h"
#include "mcrl2/lps/linearise.h"
#include "mcrl2/lps/linearise_allow_block.h"
#include "mcrl2/lps/linearise_utility.h"
#include "mcrl2/lps/linearise_rename.h"
#include "mcrl2/lps/linearise_communication.h"
#include "mcrl2/lps/replace_capture_avoiding_with_an_identifier_generator.h"
#include "mcrl2/lps/sumelm.h"

#include "mcrl2/lps/detail/ultimate_delay.h"

// Process libraries.
#include "mcrl2/process/alphabet_reduce.h"
#include "mcrl2/process/balance_nesting_depth.h"
#include "mcrl2/process/process_expression.h"


// For Aterm library extension functions
using namespace atermpp;
using namespace mcrl2;
using namespace mcrl2::core;
using namespace mcrl2::core::detail;
using namespace mcrl2::data;
using namespace mcrl2::data::detail;
using namespace mcrl2::lps;
using namespace mcrl2::process;

/*  Preamble */

using processstatustype = enum
{
  unknown,
  mCRL,
  mCRLdone,
  mCRLbusy,
  mCRLlin,
  pCRL,
  multiAction,
  GNF,
  GNFalpha,
  GNFbusy,
  error
};

/**************** Definitions of object class  ***********************/

using objecttype = enum
{
  none,
  _map,
  func,
  act,
  proc,
  variable_,
  sorttype,
  multiact
};

class objectdatatype
{
  public:
    identifier_string objectname;
    process::action_label_list multi_action_names;
    bool constructor=false;
    process_expression representedprocess;
    process_identifier process_representing_action; /* for actions target sort is used to
                                                       indicate the process representing this action. */
    process_expression processbody;
    std::set <variable> free_variables;
    bool free_variables_defined=false;
    variable_list parameters;
    variable_list old_parameters;
    processstatustype processstatus=unknown;
    objecttype object=none;
    bool canterminate=false;
    bool containstime=false;

    objectdatatype()=default;
    objectdatatype(const objectdatatype& o)=default;
    objectdatatype& operator=(const objectdatatype& o)=default;
    ~objectdatatype()=default;
};


class specification_basic_type
{
  public:
    process::action_label_list acts;     /* storage place for actions */
    std::set < variable > global_variables; /* storage place for free variables occurring
                                   in processes ranging over data */
    variable_list initdatavars; /* storage place for free variables in
                                   init clause */
    data_specification data;    /* contains the data specification for the current process.  */

    specification_basic_type(const specification_basic_type& )=delete;
    specification_basic_type& operator=(const specification_basic_type&)=delete;

  private:
    class stackoperations;
    class stacklisttype;
    class enumtype;
    class enumeratedtype;

    std::vector < process_equation > procs;
    /* storage place for processes,
       uses alt, seq, par, lmer, cond,sum,
       com, bound, at, name, delta,
       tau, hide, rename, encap */
    mcrl2::data::rewriter rewr; /* The rewriter used while linearising */
    action terminationAction;   /* A list of length one with the action that denotes termination */
    process_identifier terminatedProcId; /* A process identifier of which the body consists of the termination
                                            action */
    process_identifier tau_process;
    process_identifier delta_process;
    std::vector < process_identifier > seq_varnames; /* Contains names of processes which represent a sequence
                                                            of process variables */
    std::vector < std::vector < process_instance_assignment > > representedprocesses; /* contains the sequences of process
                                                         instances that are represented by the variables in seq_varnames */
    t_lin_options options;
    bool timeIsBeingUsed = false;
    bool stochastic_operator_is_being_used = false;
    bool fresh_equation_added = false;
    std::map < aterm, objectdatatype > objectdata; // It is important to guarantee that the objects will not
                                                   // be moved to another place when the object data structure grows. This
                                                   // is because objects in this datatype  are passed around by reference.

    set_identifier_generator fresh_identifier_generator;
    std::vector < enumeratedtype > enumeratedtypes;
    stackoperations* stack_operations_list;

  public:
    specification_basic_type(const process::action_label_list& as,
                             const std::vector< process_equation >& ps,
                             const variable_list& idvs,
                             const data_specification& ds,
                             const std::set < data::variable >& glob_vars,
                             const t_lin_options& opt,
                             const process_specification& procspec):
      acts(),
      global_variables(glob_vars),
      data(ds),
      rewr(data,opt.rewrite_strategy),
      options(opt)
    {
      // find_identifiers does not find the identifiers in the enclosed data specification.
      fresh_identifier_generator.add_identifiers(process::find_identifiers(procspec));
      // So, the identifiers in the data type must be added explicitly.
      fresh_identifier_generator.add_identifiers(data::find_identifiers(ds.equations()));
      fresh_identifier_generator.add_identifiers(data::find_identifiers(ds.user_defined_aliases()));
      fresh_identifier_generator.add_identifiers(data::find_identifiers(ds.sorts()));
      fresh_identifier_generator.add_identifiers(data::find_identifiers(ds.constructors()));
      fresh_identifier_generator.add_identifiers(data::find_identifiers(ds.mappings()));

      stack_operations_list=nullptr;
      acts=as;
      storeact(acts);
      procs=ps;
      storeprocs(procs);
      initdatavars=idvs;
      // The terminationAction and the terminatedProcId must be defined after initialisation of
      // data as otherwise fresh name does not work properly.
      terminationAction=action(action_label(fresh_identifier_generator("Terminate"),sort_expression_list()),data_expression_list());
      terminatedProcId=process_identifier(fresh_identifier_generator("Terminated**"), variable_list());
// /* Changed delta() to DeltaAtZero on 24/12/2006. Moved back in spring 2007, as this introduces unwanted time constraints. */
      insert_process_declaration(
        terminatedProcId,
        variable_list(),
        seq(terminationAction,delta()),
        pCRL,false,false);
      delta_process=newprocess(variable_list(),delta(),pCRL,false,false);
      tau_process=newprocess(variable_list(),tau(),pCRL,true,false);
    }

    ~specification_basic_type()
    {
      while (stack_operations_list!=nullptr)
      {
        stackoperations* temp=stack_operations_list->next;
        delete stack_operations_list;
        stack_operations_list=temp;
      }

    }

  private:

    static
    data_expression real_times_optimized(const data_expression& r1, const data_expression& r2)
    {
      if (r1==sort_real::real_zero() || r2==sort_real::real_zero())
      {
        return sort_real::real_zero();
      }
      if (r1==sort_real::real_one())
      {
        return r2;
      }
      if (r2==sort_real::real_one())
      {
        return r1;
      }
      return sort_real::times(r1,r2);
    }

    static
    process_expression delta_at_zero()
    {
      return at(delta(), sort_real::real_(0));
    }

    bool isDeltaAtZero(const process_expression& t)
    {
      if (!is_at(t))
      {
        return false;
      }
      if (!is_delta(at(t).operand()))
      {
        return false;
      }
      return RewriteTerm(at(t).time_stamp())==sort_real::real_one();
    }

    /***************** temporary helper function to compare substitutions ******************/

    template <class Expression, class Substitution>
    Expression replace_variables_capture_avoiding_alt(const Expression& e, Substitution& sigma)
    {
      return process::replace_variables_capture_avoiding_with_an_identifier_generator(e, sigma, fresh_identifier_generator);
    }


    /*****************  retrieve basic objects  ******************/

    void detail_check_objectdata(const process_identifier& o) const
    {
      if (objectdata.count(o)==0)
      {
        throw mcrl2::runtime_error("Fail to recognize " + process::pp(process_identifier(o)) + 
                                   ". Most likely due to unguarded recursion in a process equation. ");
      }
    }

    objectdatatype& objectIndex(const process_identifier& o)
    {
      detail_check_objectdata(o);
      assert(objectdata.find(o)->second.parameters==o.variables());
      return objectdata.find(o)->second;
    }

    const objectdatatype& objectIndex(const process_identifier& o) const
    {
      detail_check_objectdata(o);
      assert(objectdata.find(o)->second.parameters==o.variables());
      return objectdata.find(o)->second;
    }

    void addString(const identifier_string& str)
    {
      fresh_identifier_generator.add_identifier(str);
    }

    static
    process_expression action_list_to_process(const action_list& ma)
    {
      if (ma.size()==0)
      {
        return tau();
      }
      if (ma.size()==1)
      {
        return ma.front();
      }
      return process::sync(ma.front(),action_list_to_process(ma.tail()));
    }

    /// Convert the process expression to an action list.
    /// Note: The resulting actions are not sorted.
    static
    action_list to_action_list(const process_expression& p)
    {
      if (is_tau(p))
      {
        return action_list();
      }

      if (is_action(p))
      {
        return action_list({ action(p) });
      }

      if (is_sync(p))
      {
        return to_action_list(process::sync(p).left())+to_action_list(process::sync(p).right());
      }
      assert(0);
      return action_list();
    }

    /// Convert the process expression to a sorted action list.
    static
    action_list to_sorted_action_list(const process_expression& p)
    {
      return sort_actions(to_action_list(p));
    }

    static
    action_label_list getnames(const process_expression& multiAction)
    {
      if (is_action(multiAction))
      {
        return action_label_list({ action(multiAction).label() });
      }
      assert(is_sync(multiAction));
      return getnames(process::sync(multiAction).left())+getnames(process::sync(multiAction).right());
    }

    // Returns a list of variables with the same sort as the expressions in the list.
    // If the expression is a variable not occurring in the occurs_set that variable
    // is used.
    variable_list make_parameters_rec(const data_expression_list& l,
                                      std::set < variable>& occurs_set)
    {
      variable_list result;
      for (const data_expression& e: l)
      {
        /* if the current argument of the multi-action is a variable that does
         not occur in result, use this variable. This is advantageous, when joining
         processes to one linear process where variable names are joined. If this
         is not being done (as happened before 4/1/2008) very long lists of parameters
         can occur when linearising using regular2 */
        if (is_variable(e) && std::find(occurs_set.begin(),occurs_set.end(),e)==occurs_set.end())
        {
          const variable& v = atermpp::down_cast<variable>(e);
          result.push_front(v);
          occurs_set.insert(v);
        }
        else
        {
          result.push_front(variable(get_fresh_variable("a",e.sort())));
        }
      }
      return reverse(result);
    }

    variable_list getparameters_rec(const process_expression& multiAction,
                                    std::set < variable>& occurs_set)
    {
      if (is_action(multiAction))
      {
        return make_parameters_rec(action(multiAction).arguments(),occurs_set);
      }
      assert(is_sync(multiAction));
      return getparameters_rec(process::sync(multiAction).left(),occurs_set)+
             getparameters_rec(process::sync(multiAction).right(),occurs_set);
    }

    variable_list getparameters(const process_expression& multiAction)
    {
      std::set < variable > occurs_set;
      return getparameters_rec(multiAction,occurs_set);
    }

    static
    data_expression_list getarguments(const action_list& multiAction)
    {
      data_expression_list result;
      for (const action& a: multiAction)
      {
        result=reverse(a.arguments()) + result;
      }
      return reverse(result);
    }

    static
    action_list makemultiaction(const process::action_label_list& actionIds, const data_expression_list& args)
    {
      action_list result;
      data_expression_list::const_iterator e_walker=args.begin();
      for (const process::action_label& l: actionIds)
      {
        std::size_t arity=l.sorts().size();
        data_expression_list temp_args;
        for (std::size_t i=0 ; i< arity; ++i,++e_walker)
        {
          assert(e_walker!=args.end());
          temp_args.push_front(*e_walker);
        }
        temp_args=reverse(temp_args);
        result.push_front(action(l,temp_args));
      }
      assert(e_walker==args.end());
      return reverse(result);
    }

    objectdatatype& addMultiAction(const process_expression& multiAction, bool& isnew)
    {
      const process::action_label_list actionnames=getnames(multiAction);

      isnew=(objectdata.count(actionnames)==0);

      if (isnew)
      {
        objectdatatype object;

        // tempvar is needed as objectdata can change during a call
        // of getparameters.
        const variable_list templist=getparameters(multiAction);
        object.parameters=templist;
        object.object=multiact;
        // must separate assignment below as
        // objectdata may change as a side effect of make
        // multiaction.
        const action_list tempvar=makemultiaction(actionnames, variable_list_to_data_expression_list(object.parameters));
        object.processbody=action_list_to_process(tempvar);
        object.free_variables=std::set<variable>(object.parameters.begin(), object.parameters.end());
        object.free_variables_defined=true;

        objectdata[actionnames]=object;
      }
      return objectdata.find(actionnames)->second;
    }

    const std::set<variable>& get_free_variables(objectdatatype& object)
    {
      if (!object.free_variables_defined)
      {
        object.free_variables=process::find_free_variables(object.processbody);
        object.free_variables_defined=true;
      }
      return object.free_variables;
    }

    void insertvariable(const variable& var, const bool mustbenew)
    {
      addString(var.name());

      if (objectdata.count(var.name())>0  && mustbenew)
      {
        throw mcrl2::runtime_error("Variable " + data::pp(var) + " already exists. ");
      }

      objectdatatype object;
      object.objectname=var.name();
      object.object=variable_;
      objectdata[var.name()]=object;
    }

    void insertvariables(const variable_list& vars, const bool mustbenew)
    {
      for (const variable& v: vars)
      {
        insertvariable(v,mustbenew);
      }
    }

    template <class SUBSTITUTION>
    std::set<data::variable> sigma_variables(const SUBSTITUTION& sigma)
    {
      std::set<data::variable> result;
      for (typename SUBSTITUTION::const_iterator i = sigma.begin(); i != sigma.end(); ++i)
      {
        std::set<data::variable> V = data::find_free_variables(i->second);
        V.erase(i->first);
        result.insert(V.begin(), V.end());
      }
      return result;
    }


    /************ upperpowerof2 *********************************************/

    static
    std::size_t upperpowerof2(std::size_t i)
    /* function yields n for the smallest value n such that
       2^n>=i. This constitutes the number of bits necessary
       to represent a number smaller than i. i is assumed to
       be at least 1. */
    {
      std::size_t n=0;
      std::size_t powerof2=1;
      for (; powerof2< i ; n++)
      {
        powerof2=2*powerof2;
      }
      return n;
    }

    data_expression RewriteTerm(const data_expression& t)
    {
      if (!options.norewrite)
      {
        if (fresh_equation_added)
        {
          rewr=rewriter(data,options.rewrite_strategy);
          fresh_equation_added=false;
        }
        return rewr(t);
      }
      return t;
    }

    data_expression_list RewriteTermList(const data_expression_list& t)
    {
      data_expression_vector v;
      for(const data_expression& d: t)
      {
        v.push_back(RewriteTerm(d));
      }
      return data_expression_list(v.begin(),v.end());
    }

    assignment_list rewrite_assignments(const assignment_list& t)
    {
      assignment_vector v;
      for(const assignment& a: t)
      {
        v.emplace_back(a.lhs(), RewriteTerm(a.rhs()));
      }
      return assignment_list(v.begin(),v.end());
    }

    action RewriteAction(const action& t)
    {
      return action(t.label(),RewriteTermList(t.arguments()));
    }

    process_instance_assignment RewriteProcess(const process_instance_assignment& t)
    {
      return process_instance_assignment(t.identifier(),rewrite_assignments(t.assignments()));
    }

    process_expression RewriteMultAct(const process_expression& t)
    {
      if (is_tau(t))
      {
        return t;
      }

      if (is_action(t))
      {
        return RewriteAction(action(t));
      }

      assert(is_sync(t));  // A multi action is a sequence of actions with a sync operator in between.
      return process::sync(RewriteMultAct(process::sync(t).left()),RewriteMultAct(process::sync(t).right()));
    }


    process_expression pCRLrewrite(const process_expression& t)
    {
      if (options.norewrite)
      {
        return t;
      }

      if (is_if_then(t))
      {
        const data_expression new_cond=RewriteTerm(down_cast<if_then>(t).condition());
        const process_expression new_then_case=pCRLrewrite(down_cast<if_then>(t).then_case());
        if (new_cond==sort_bool::true_())
        {
          return new_then_case;
        }
        return if_then(new_cond,new_then_case);
      }

      if (is_seq(t))
      {
        /* only one summand is needed */
        return seq(
                 pCRLrewrite(seq(t).left()),
                 pCRLrewrite(seq(t).right()));
      }

      if (is_at(t))
      {
        const data_expression& atTime=RewriteTerm(down_cast<at>(t).time_stamp());
        const process_expression t1=pCRLrewrite(down_cast<at>(t).operand());
        return at(t1,atTime);
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(t);
        return stochastic_operator(sto.variables(),RewriteTerm(sto.distribution()),pCRLrewrite(sto.operand()));
      }

      if (is_delta(t) || is_tau(t))
      {
        return t;
      }

      if (is_action(t))
      {
        return RewriteAction(atermpp::down_cast<process::action>(t));
      }

      if (is_process_instance_assignment(t))
      {
        return RewriteProcess(process_instance_assignment(t));
      }

      if (is_sync(t))
      {
        return RewriteMultAct(t);
      }

      assert(0); // Expected a term in pCRL format, using only basic process operators.
      return t;
    }

    /************ storeact ****************************************************/

    objectdatatype& insertAction(const action_label& actionId)
    {
      if (objectdata.count(actionId)>0)
      {
        throw mcrl2::runtime_error("Action " + process::pp(actionId) + " is added twice. This is an internal error in the lineariser. Please report. ");
      }

      const identifier_string& str=actionId.name();
      addString(str);

      objectdatatype object;
      object.objectname=str;
      object.object=act;
      object.process_representing_action=process_identifier();

      objectdata[actionId]=object;
      return objectdata.find(actionId)->second;
    }

    void storeact(const process::action_label_list& acts)
    {
      for (const process::action_label& l: acts)
      {
        insertAction(l);
      }
    }

    /************ storeprocs *************************************************/

    objectdatatype& insert_process_declaration(
      const process_identifier& procId,
      const variable_list& parameters,
      const process_expression& body,
      processstatustype s,
      const bool canterminate,
      const bool containstime)
    {
      assert(procId.variables().size()==parameters.size());
      const std::string str=procId.name();
      addString(str);

      if (objectdata.count(procId)>0)
      {
        throw mcrl2::runtime_error("Process " + process::pp(procId) + " is added twice. This is an internal error in the lineariser. Please report. ");
      }

      objectdatatype object;
      object.objectname=procId.name();
      object.object=proc;
      object.processbody=body;
      object.free_variables_defined=false;
      object.canterminate=canterminate;
      object.containstime=containstime;
      object.processstatus=s;
      object.parameters=parameters;
      insertvariables(parameters,false);
      objectdata[procId]=object;
      return objectdata.find(procId)->second;
    }

    void storeprocs(const std::vector< process_equation >& procs)
    {
      for (const process_equation& e: procs)
      {
        insert_process_declaration(
          e.identifier(),
          e.formal_parameters(),
          e.expression(),
          unknown,false,false);
      }
    }

    bool searchProcDeclaration(
      const variable_list& parameters,
      const process_expression& body,
      const processstatustype s,
      const bool canterminate,
      const bool containstime,
      process_identifier& p) const
    {
      for(const std::pair<const aterm,objectdatatype>& d: objectdata)
      {
        if (d.second.object==proc &&
            d.second.parameters==parameters &&
            d.second.processbody==body &&
            d.second.canterminate==canterminate &&
            d.second.containstime==containstime &&
            d.second.processstatus==s)
        {
          p=process_identifier(d.second.objectname,d.second.parameters);
          return true;
        }
      }
      return false;
    }


    /************ storeinit *************************************************/

  public:
    process_identifier storeinit(const process_expression& init)
    {
      /* init is used as the name of the initial process,
         because it cannot occur as a string in the input */

      process_identifier initprocess(std::string("init"), variable_list());
      insert_process_declaration(initprocess,variable_list(),init,unknown,false,false);
      return initprocess;
    }

  private:

    /********** various functions on action and multi actions  ***************/



    action_list linMergeMultiActionList(const action_list& ma1, const action_list& ma2)
    {
      action_list result=ma2;
      for (const action& a: ma1)
      {
        result=insert(a,result);
      }
      return result;
    }


    action_list linMergeMultiActionListProcess(const process_expression& ma1, const process_expression& ma2)
    {
      return linMergeMultiActionList(to_action_list(ma1),to_action_list(ma2));
    }

    /************** determine_process_status ********************************/

    processstatustype determine_process_statusterm(
      const process_expression& body,  // intentionally not a reference.
      const processstatustype status)
    {
      /* In this procedure it is determined whether a process
         is of type mCRL, pCRL or a multiAction. pCRL processes
         occur strictly within mCRL processes, and multiActions
         occur strictly within pCRL processes. Processes that pass
         this procedure can be linearised. Bounded initialisation,
         the leftmerge and synchronization merge on the highest
         level are filtered out. */
      if (is_choice(body))
      {
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("Choice operator occurs in a multi-action in " + process::pp(body) + ". The lineariser cannot handle such a pattern. ");
        }
        const processstatustype s1=determine_process_statusterm(choice(body).left(),pCRL);
        const processstatustype s2=determine_process_statusterm(choice(body).right(),pCRL);
        if ((s1==mCRL)||(s2==mCRL))
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs within the scope of a choice operator in " + process::pp(body) +
                                      ". The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_seq(body))
      {
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("Sequential operator occurs in a multi-action in " + process::pp(body) +". The lineariser cannot handle such a pattern. ");
        }
        const processstatustype s1=determine_process_statusterm(seq(body).left(),pCRL);
        const processstatustype s2=determine_process_statusterm(seq(body).right(),pCRL);
        if ((s1==mCRL)||(s2==mCRL))
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of a sequential operator in " + process::pp(body) +". "
                                       + "The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_merge(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("The parallel operator occurs in the scope of recursion, or the condition, sequential or choice operator in " +
                                      process::pp(body) + ". The lineariser cannot handle such processes. ");
        }
        determine_process_statusterm(process::merge(body).left(),mCRL);
        determine_process_statusterm(process::merge(body).right(),mCRL);
        return mCRL;
      }

      if (is_left_merge(body))
      {
        throw mcrl2::runtime_error("Cannot linearise because the specification contains a leftmerge. ");
      }

      if (is_if_then(body))
      {
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("If-then occurs in a multi-action in " + process::pp(body) + ". The lineariser cannot linearise this. ");
        }
        const processstatustype s1=determine_process_statusterm(if_then(body).then_case(),pCRL);
        if (s1==mCRL)
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of the if-then operator in " + process::pp(body) + ". "
                                     + "The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_if_then_else(body))
      {
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("If-then-else occurs in a multi-action in " + process::pp(body) +".");
        }
        const processstatustype s1=determine_process_statusterm(if_then_else(body).then_case(),pCRL);
        const processstatustype s2=determine_process_statusterm(if_then_else(body).else_case(),pCRL);
        if ((s1==mCRL)||(s2==mCRL))
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of the if-then-else operator in " + process::pp(body) +
                                     ". " + "The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_sum(body))
      {
        /* insert the variable names of variables, to avoid
           that this variable name will be reused later on */
        insertvariables(sum(body).variables(),false);
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("The sum operator occurs within a multi-action in " + process::pp(body) + ". "
                                      + "The lineariser cannot handle such processes. ");
        }
        const processstatustype s1=determine_process_statusterm(sum(body).operand(),pCRL);
        if (s1==mCRL)
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of the sum operator in " + process::pp(body) + ". "
                                     "The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_stochastic_operator(body))
      {
        stochastic_operator_is_being_used=true;
        /* insert the variable names of variables, to avoid
           that this variable name will be reused later on */
        const stochastic_operator& sto=down_cast<const stochastic_operator>(body);
        insertvariables(sto.variables(),false);
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("Stochastic operator occurs within a multi-action in " + process::pp(body) +".");
        }
        const processstatustype s1=determine_process_statusterm(sto.operand(),pCRL);
        if (s1==mCRL)
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of the stochastic operator in " + process::pp(body) + ". "
                                     + "The lineariser cannot handle such processes. ");
        }
        return s1;
      }

      if (is_comm(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("The communication operator occurs in the scope of recursion, the condition, the sequential operation or the choice in "
                                      + process::pp(body) + ". The lineariser cannot linearise such processes. ");
        }
        determine_process_statusterm(comm(body).operand(),mCRL);
        return mCRL;
      }

      if (is_bounded_init(body))
      {
        throw mcrl2::runtime_error("Cannot linearise a specification with the bounded initialization operator.");
      }

      if (is_at(body))
      {
        timeIsBeingUsed = true;
        if (status==multiAction)
        {
          throw mcrl2::runtime_error("Time operator occurs in a multi-action in " + process::pp(body) +".");
        }
        const processstatustype s1=determine_process_statusterm(at(body).operand(),pCRL);
        if (s1==mCRL)
        {
          throw mcrl2::runtime_error("An operator ||, allow, hide, rename, or comm occurs in the scope of a time operator in " + process::pp(body) + "."
                                     + "The lineariser cannot handle such processes. ");
        }
        return pCRL;
      }

      if (is_sync(body))
      {
        const processstatustype s1=determine_process_statusterm(process::sync(body).left(),pCRL);
        const processstatustype s2=determine_process_statusterm(process::sync(body).right(),pCRL);
        if ((s1!=multiAction)||(s2!=multiAction))
        {
          throw mcrl2::runtime_error("Other objects than multi-actions occur in the scope of a synch operator in " + process::pp(body) +".");
        }
        return multiAction;
      }

      if (is_action(body))
      {
        return multiAction;
      }

      if (is_process_instance_assignment(body))
      {
        determine_process_status(process_instance_assignment(body).identifier(),status);
        return status;
      }

      if (is_delta(body))
      {
        return pCRL;
      }

      if (is_tau(body))
      {
        return multiAction;
      }

      if (is_hide(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("A hide operator occurs in the scope of recursion, or a condition, choice or sequential operator in " + process::pp(body) +".");
        }
        determine_process_statusterm(hide(body).operand(),mCRL);
        return mCRL;
      }

      if (is_rename(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("A rename operator occurs in the scope of recursion, or a condition, choice or sequential operator in " + process::pp(body) +".");
        }
        determine_process_statusterm(process::rename(body).operand(),mCRL);
        return mCRL;
      }

      if (is_allow(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("An allow operator occurs in the scope of recursion, or a condition, choice or sequential operator in " + process::pp(body) +".");
        }
        determine_process_statusterm(allow(body).operand(),mCRL);
        return mCRL;
      }

      if (is_block(body))
      {
        if (status!=mCRL)
        {
          throw mcrl2::runtime_error("A block operator occurs in the scope of recursion, or a condition, choice or sequential operator in " + process::pp(body) +".");
        }
        determine_process_statusterm(block(body).operand(),mCRL);
        return mCRL;
      }

      throw mcrl2::runtime_error("Process has unexpected format (2) " + process::pp(body) +".");
      return error;
    }


    void determine_process_status(
      const process_identifier& procDecl,
      const processstatustype status)
    {
      objectdatatype& object=objectIndex(procDecl);
      processstatustype s=object.processstatus;

      if (s==unknown)
      {
        object.processstatus=status;
        if (status==pCRL)
        {
          determine_process_statusterm(object.processbody,pCRL);
          return;
        }
        /* status==mCRL */
        s=determine_process_statusterm(object.processbody,mCRL);
        if (s!=status)
        {
          /* s==pCRL and status==mCRL */
          object.processstatus=s;
          determine_process_statusterm(object.processbody,pCRL);
        }
      }
      if (s==mCRL)
      {
        if (status==pCRL)
        {
          object.processstatus=pCRL;
          determine_process_statusterm(object.processbody,pCRL);
        }
      }
    }

    /***********  collect pcrlprocessen **********************************/

    void collectPcrlProcesses_term(const process_expression& body,  // Intentionally not a reference.
                                   std::vector <process_identifier>& pcrlprocesses,
                                   std::set <process_identifier>& visited)
    {
      if (is_if_then(body))
      {
        collectPcrlProcesses_term(if_then(body).then_case(),pcrlprocesses,visited);
        return;
      }

      if (is_if_then_else(body))
      {
        collectPcrlProcesses_term(if_then_else(body).then_case(),pcrlprocesses,visited);
        collectPcrlProcesses_term(if_then_else(body).else_case(),pcrlprocesses,visited);
        return;
      }

      if (is_choice(body))
      {
        collectPcrlProcesses_term(choice(body).left(),pcrlprocesses,visited);
        collectPcrlProcesses_term(choice(body).right(),pcrlprocesses,visited);
        return ;
      }

      if (is_seq(body))
      {
        collectPcrlProcesses_term(seq(body).left(),pcrlprocesses,visited);
        collectPcrlProcesses_term(seq(body).right(),pcrlprocesses,visited);
        return ;
      }

      if (is_merge(body))
      {
        collectPcrlProcesses_term(process::merge(body).left(),pcrlprocesses,visited);
        collectPcrlProcesses_term(process::merge(body).right(),pcrlprocesses,visited);
        return ;
      }

      if (is_sync(body))
      {
        collectPcrlProcesses_term(process::sync(body).left(),pcrlprocesses,visited);
        collectPcrlProcesses_term(process::sync(body).right(),pcrlprocesses,visited);
        return ;
      }

      if (is_sum(body))
      {
        collectPcrlProcesses_term(sum(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_stochastic_operator(body))
      {
        collectPcrlProcesses_term(stochastic_operator(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_at(body))
      {
        collectPcrlProcesses_term(at(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_process_instance_assignment(body))
      {
        collectPcrlProcesses(process_instance_assignment(body).identifier(),pcrlprocesses,visited);
        return;
      }

      if (is_hide(body))
      {
        collectPcrlProcesses_term(hide(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_rename(body))
      {
        collectPcrlProcesses_term(process::rename(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_allow(body))
      {
        collectPcrlProcesses_term(allow(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_block(body))
      {
        collectPcrlProcesses_term(block(body).operand(),pcrlprocesses,visited);
        return;
      }

      if (is_comm(body))
      {
        collectPcrlProcesses_term(comm(body).operand(),pcrlprocesses,visited);
        return;
      }

      if ((is_delta(body))||
          (is_tau(body))||
          (is_action(body)))
      {
        return;
      }

      throw mcrl2::runtime_error("process has unexpected format (1) " + process::pp(body) +".");
    }

    void collectPcrlProcesses(
      const process_identifier& procDecl,
      std::vector <process_identifier>& pcrlprocesses,
      std::set <process_identifier>& visited)
    {
      if (visited.count(procDecl)==0)
      {
        visited.insert(procDecl);
        objectdatatype& object=objectIndex(procDecl);
        if (object.processstatus==pCRL)
        {
          pcrlprocesses.push_back(procDecl);
        }
        collectPcrlProcesses_term(object.processbody,pcrlprocesses,visited);
      }
    }

    void collectPcrlProcesses(
      const process_identifier& procDecl,
      std::vector <process_identifier>& pcrlprocesses)
    {
      std::set <process_identifier> visited;
      collectPcrlProcesses(procDecl, pcrlprocesses, visited);
    }

    /****************  occursinterm *** occursintermlist ***********/

    void filter_vars_by_term(
      const data_expression& t,
      const std::set < variable >& vars_set,
      std::set < variable >& vars_result_set)
    {
      if (is_variable(t))
      {
        const variable& v = atermpp::down_cast<variable>(t);
        if (vars_set.find(v)!=vars_set.end())
        {
          vars_result_set.insert(v);
        }
        return;
      }

      if (is_function_symbol(t))
      {
        return;
      }

      if (is_machine_number(t))
      {
        return;
      }

      if (is_abstraction(t))
      {
        // mCRL2log(mcrl2::log::warning) << "filtering of variables expression with binders" << std::endl;
        return;
      }

      if (is_where_clause(t))
      {
        // mCRL2log(mcrl2::log::warning) << "filtering of variables expression with where clause" << std::endl;
        return;
      }

      if (!is_application(t))
      {
        mCRL2log(mcrl2::log::error) << "term of unexpected type " << t << std::endl;
      }

      assert(is_application(t));

      const application& a=atermpp::down_cast<const application>(t);
      filter_vars_by_term(a.head(),vars_set,vars_result_set);
      filter_vars_by_termlist(a.begin(),a.end(),vars_set,vars_result_set);
    }

    static
    bool occursintermlist(const variable& var, const data_expression_list& r)
    {
      for (const data_expression& d: r)
      {
        if (occursinterm(d, var))
        {
          return true;
        }
      }
      return false;
    }

    bool occursintermlist(const variable& var, const assignment_list& r, const process_identifier& proc_name) const
    {
      std::set<variable> assigned_variables;
      for (const assignment& l: r)
      {
        if (occursinterm(l.rhs(), var))
        {
          return true;
        }
       assigned_variables.insert(l.lhs());
      }
      // Check whether x does not occur in the assignment list. Then variable x is assigned to
      // itself, and it occurs in the process.
      for (const variable& v: objectIndex(proc_name).parameters)
      {
        if (var==v)
        {
          if (assigned_variables.count(var)==0) // This variable is not assigned, so it does occur!
          {
            return true;
          }
        }
      }
      return false;
    }

    template <typename Iterator>
    void filter_vars_by_termlist(
      Iterator begin,
      const Iterator& end,
      const std::set < variable >& vars_set,
      std::set < variable >& vars_result_set)
    {
      for (; begin != end; ++begin)
      {
        filter_vars_by_term(*begin,vars_set,vars_result_set);
      }
    }

    void filter_vars_by_multiaction(
      const action_list& multiaction,
      const std::set < variable >& vars_set,
      std::set < variable >& vars_result_set)
    {
      for (const action& ma: multiaction)
      {
        filter_vars_by_termlist(ma.arguments().begin(), ma.arguments().end(),vars_set,vars_result_set);
      }
      return;
    }

    void filter_vars_by_assignmentlist(
      const assignment_list& assignments,
      const variable_list& parameters,
      const std::set < variable >& vars_set,
      std::set < variable >& vars_result_set)
    {
      const auto& l= parameters | std::views::transform([](const variable& v) { return atermpp::down_cast<data_expression>(v); });
      filter_vars_by_termlist(l.begin(),l.end(),vars_set,vars_result_set);
      for (const assignment& a: assignments)
      {
        filter_vars_by_term(a.rhs(),vars_set,vars_result_set);
      }
    }

    bool occursinpCRLterm(const variable& var,
                          const process_expression& p,
                          const bool strict)
    {
      if (is_choice(p))
      {
        return occursinpCRLterm(var,choice(p).left(),strict)||
               occursinpCRLterm(var,choice(p).right(),strict);
      }
      if (is_seq(p))
      {
        return occursinpCRLterm(var,seq(p).left(),strict)||
               occursinpCRLterm(var,seq(p).right(),strict);
      }
      if (is_if_then(p))
      {
        return occursinterm(if_then(p).condition(), var) ||
               occursinpCRLterm(var,if_then(p).then_case(),strict);
      }

      if (is_sum(p))
      {
        if (strict)
        {
          return occursintermlist(var,variable_list_to_data_expression_list(sum(p).variables())) ||
                 occursinpCRLterm(var,sum(p).operand(),strict);
        /* below appears better? , but leads
           to errors. Should be investigated. */
        }
        else
        {
          return (!occursintermlist(var, variable_list_to_data_expression_list(sum(p).variables())))
                 && occursinpCRLterm(var, sum(p).operand(), strict);
        }
      }
      if (is_stochastic_operator(p))
      {
        const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(p);
        if (strict)
        {
          return occursintermlist(var,variable_list_to_data_expression_list(sto.variables())) ||
                      occursinterm(sto.distribution(), var) ||
                      occursinpCRLterm(var,sto.operand(),strict);
        }
        else
        {
          return (!occursintermlist(var,variable_list_to_data_expression_list(sto.variables()))) &&
                      (occursinterm(sto.distribution(), var) ||
                       occursinpCRLterm(var,sto.operand(),strict));
        }
      }
      if (is_process_instance_assignment(p))
      {
        return occursintermlist(var,process_instance_assignment(p).assignments(),process_instance_assignment(p).identifier());
      }
      if (is_action(p))
      {
        return occursintermlist(var,action(p).arguments());
      }
      if (is_sync(p))
      {
        return occursinpCRLterm(var,process::sync(p).left(),strict)||
               occursinpCRLterm(var,process::sync(p).right(),strict);
      }
      if (is_at(p))
      {
        return occursinterm(at(p).time_stamp(), var) ||
               occursinpCRLterm(var,at(p).operand(),strict);
      }
      if (is_delta(p))
      {
        return false;
      }
      if (is_tau(p))
      {
        return false;
      }
      throw mcrl2::runtime_error("unexpected process format in occursinCRLterm " + process::pp(p));
      return false;
    }

    template <class MutableSubstitution>
    void alphaconvertprocess(
      variable_list& sumvars,
      MutableSubstitution& sigma,
      const process_expression& p)
    {
      /* This function replaces the variables in sumvars
         by unique ones if these variables occur in the process expression p.
         The substitution sigma contains the renaming. */
      variable_vector newsumvars;

      for (const variable& var: sumvars)
      {
        if (occursinpCRLterm(var,p,true))
        {
          const variable newvar=get_fresh_variable(var.name(),var.sort());
          newsumvars.push_back(newvar);
          sigma[var]=newvar;
        }
        else
        {
          newsumvars.push_back(var);
        }
      }
      sumvars=variable_list(newsumvars.begin(), newsumvars.end());
    }

    template <class MutableSubstitution>
    void alphaconvert(
      variable_list& sumvars,
      MutableSubstitution& sigma,
      const variable_list& occurvars,
      const data_expression_list& occurterms)
    {
      /* This function replaces the variables in sumvars
         by unique ones if these variables occur in occurvars
         or occurterms. It extends rename_vars and rename
         terms to rename the replaced variables to new ones. */
      variable_list newsumvars;

      for (const variable& var: sumvars)
      {
        if (occursintermlist(var,variable_list_to_data_expression_list(occurvars)) ||
            occursintermlist(var,occurterms))
        {
          const variable newvar=get_fresh_variable(var.name(),var.sort());
          newsumvars.push_front(newvar);
          sigma[var]=newvar;
        }
        else
        {
          newsumvars.push_front(var);
        }
      }
      sumvars=reverse(newsumvars);
    }

    std::set< variable > find_free_variables_process(const process_expression& p)
    {
      std::set<variable> free_variables_in_p_new;
      free_variables_in_p_new=process::find_free_variables(p); 
      return free_variables_in_p_new; 
    }

    /* Remove assignments that do not appear in the parameter list. */
    static
    assignment_list filter_assignments(const assignment_list& assignments, const variable_list& parameters)
    {
      assignment_vector result;
      for(const assignment& a: assignments)
      {
        if (std::find(parameters.begin(),parameters.end(),a.lhs())!=parameters.end())   // found.
        {
          result.push_back(a);
        }
      }
      return assignment_list(result.begin(),result.end());
    }

    /* Check whether the assignments occur in the same order in the list of parameters */
    static bool check_assignment_list(
              const assignment_list& assignments,
              const variable_list& parameters)
    {
      assignment_list::iterator i=assignments.begin();
      for(const variable& v: parameters)
      {
        if (i!=assignments.end() && v==i->lhs())
        {
          ++i;
        }
      }
      return i==assignments.end();

    }

    /******************* substitute *****************************************/


    template <class Substitution>
    assignment_list substitute_assignmentlist(
      const assignment_list& assignments,
      const variable_list& parameters,
      const bool replacelhs,
      const bool replacerhs,
      Substitution& sigma)
    {
      assert(check_assignment_list(assignments, parameters));
      /* precondition: the variables in the assignment occur in
         the same sequence as in the parameters, which stands for the
         total list of parameters.

         This function replaces the variables in vars by the terms in terms
         in the right hand side of the assignments if replacerhs holds, and
         in the lefthandside of an assignment if replacelhs holds. If for some variable
         occuring in the parameterlist no assignment is present, whereas
         this variable occurs in vars, an assignment for it is added.

         It is not possible to use standard substitution functions to replace substitute_assignmentlis
         because they do not take the parameters of processes into account.
         For instance consider a process P(b:D)=... of which an instance P() exists.
         If the substitution sigma(b)=t is applied to P() the result should be P(b=t).
         The standard substitutions do not take this parameterlist into account, as it stands.
      */

      assert(replacelhs || replacerhs);
      if (parameters.empty())
      {
        assert(assignments.empty());
        return assignments;
      }

      const variable& parameter=parameters.front();

      if (!assignments.empty())
      {
        const assignment& ass=assignments.front();
        variable lhs=ass.lhs();
        if (parameter==lhs)
        {
          /* The assignment refers to parameter par. Substitute its
             left and righthandside */
          data_expression rhs=ass.rhs();

          if (replacelhs)
          {
            lhs = atermpp::down_cast<variable>(sigma(lhs));
          }
          if (replacerhs)
          {
            rhs=replace_variables_capture_avoiding_alt(rhs,sigma);
          }

          assignment_list result=
                   substitute_assignmentlist(
                     assignments.tail(),
                     parameters.tail(),
                     replacelhs,
                     replacerhs,
                     sigma);
          result.push_front(assignment(lhs,rhs));
          return result;
        }
      }

      /* Here the first parameter is not equal to the first
         assignment. So, we must find out whether a value
         for this variable is substituted, that is different
         from the variable, in which case an assignment must
         be added. */

      variable lhs=parameter;
      data_expression rhs=parameter;

      if (replacelhs)
      {
        lhs = atermpp::down_cast<data::variable>(sigma(lhs));
      }
      if (replacerhs)
      {
        rhs=replace_variables_capture_avoiding_alt(rhs,sigma);
      }

      if (lhs==rhs)
      {
        return substitute_assignmentlist(
                 assignments,
                 parameters.tail(),
                 replacelhs,
                 replacerhs,
                 sigma);
      }
      assignment_list result=
               substitute_assignmentlist(
                 assignments,
                 parameters.tail(),
                 replacelhs,
                 replacerhs,
                 sigma);
      result.push_front(assignment(lhs,rhs));
      return result;
    }

    // Sort the assignments, such that they have the same order as the parameters
    static
    assignment_list sort_assignments(const assignment_list& ass, const variable_list& parameters)
    {
      std::map<variable,data_expression>assignment_map;
      for(const assignment& a: ass)
      {
        assignment_map[a.lhs()]=a.rhs();
      }

      assignment_vector result;
      for(const variable& v: parameters)
      {
        const std::map<variable,data_expression>::const_iterator j=assignment_map.find(v);
        if (j!=assignment_map.end()) // found
        {
          result.emplace_back(j->first, j->second);
        }
      }
      return assignment_list(result.begin(),result.end());
    }

    bool check_valid_process_instance_assignment(
               const process_identifier& id,
               const assignment_list& assignments)
    {
      objectdatatype& object=objectIndex(id);
      variable_list parameters=object.parameters;
      for(const assignment& a: assignments)
      {
        // Every assignment must occur in the parameter list, in the right sequence.

        variable v;
        do
        {
          if (parameters.empty())
          {
            return false;
          }
          v=parameters.front();
          parameters.pop_front();
        }
        while (v!=a.lhs());

      }
      return true;
    }

    /* The function below calculates sigma(p) and replaces
       all variables that are bound by a sum in p by unique
       variables */


    /* The function below cannot be replace by replace_variables_capture_avoiding although
     * the interfaces are the same. As yet it is unclear why, but the difference shows itself
     * for instance when linearising lift3-final.mcrl2 and lift3_init.mcrl2 */
    template <class Substitution>
    process_expression substitute_pCRLproc(
      const process_expression& p,
      Substitution& sigma)
    {
      // return process::replace_variables_capture_avoiding(p, sigma, fresh_identifier_generator);
      if (is_choice(p))
      {
        process_expression left=substitute_pCRLproc(choice(p).left(),sigma);
        if (left==delta_at_zero())
        {
          return substitute_pCRLproc(choice(p).right(),sigma);
        }
        process_expression right=substitute_pCRLproc(choice(p).right(),sigma);
        if (right==delta_at_zero())
        {
          return left;
        }
        return choice(left,right);
      }
      if (is_seq(p))
      {
        process_expression q=substitute_pCRLproc(seq(p).left(),sigma);
        if (q==delta_at_zero())
        {
          return q;
        }
        return seq(q, substitute_pCRLproc(seq(p).right(),sigma));
      }
      if (is_if_then(p))
      {
        data_expression condition=RewriteTerm(replace_variables_capture_avoiding_alt(if_then(p).condition(), sigma));
        if (condition==sort_bool::false_())
        {
          return delta_at_zero();
        }
        if (condition==sort_bool::true_())
        {
          return substitute_pCRLproc(if_then(p).then_case(),sigma);
        }
        return if_then(condition,substitute_pCRLproc(if_then(p).then_case(),sigma));
      }
      if (is_if_then_else(p))
      {
        data_expression condition=RewriteTerm(replace_variables_capture_avoiding_alt(if_then_else(p).condition(), sigma));
        if (condition==sort_bool::false_())
        {
          return substitute_pCRLproc(if_then_else(p).else_case(),sigma);
        }
        if (condition==sort_bool::true_())
        {
          return substitute_pCRLproc(if_then_else(p).then_case(),sigma);
        }
        return if_then_else(
                 condition,
                 substitute_pCRLproc(if_then_else(p).then_case(),sigma),
                 substitute_pCRLproc(if_then_else(p).else_case(),sigma));
      }

      if (is_sum(p))
      {
        variable_list sumargs=sum(p).variables();
        variable_list vars;
        data_expression_list terms;

        for( std::map < variable, data_expression >::const_iterator i=sigma.begin(); i!=sigma.end(); ++i)
        {
          vars=push_back(vars,i->first);
          terms=push_back(terms,i->second);
        }

        maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma=sigma;
        alphaconvert(sumargs,local_sigma,vars,terms);

        return sum(sumargs,
                   substitute_pCRLproc(sum(p).operand(),local_sigma));
      }

      if (is_stochastic_operator(p))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(p);
        variable_list sumargs=sto.variables();
        variable_list vars;
        data_expression_list terms;

        for( std::map < variable, data_expression >::const_iterator i=sigma.begin(); i!=sigma.end(); ++i)
        {
          vars=push_back(vars,i->first);
          terms=push_back(terms,i->second);
        }

        maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma=sigma;
        alphaconvert(sumargs,local_sigma,vars,terms);

        return stochastic_operator(
                            sumargs,
                            replace_variables_capture_avoiding_alt(sto.distribution(),local_sigma),
                            substitute_pCRLproc(sto.operand(),local_sigma));
      }

      if (is_process_instance_assignment(p))
      {
        const process_instance_assignment q(p);
        objectdatatype& object=objectIndex(q.identifier());
        const variable_list parameters=object.parameters;
        const assignment_list new_assignments=substitute_assignmentlist(q.assignments(),parameters,false,true,sigma);
        assert(check_valid_process_instance_assignment(q.identifier(),new_assignments));
        return process_instance_assignment(q.identifier(),new_assignments);
      }

      if (is_action(p))
      {
        return action(action(p).label(),
                      replace_variables_capture_avoiding_alt(action(p).arguments(), sigma));
      }

      if (is_at(p))
      {
        return at(substitute_pCRLproc(at(p).operand(),sigma),
                      replace_variables_capture_avoiding_alt(at(p).time_stamp(),sigma));
      }

      if (is_delta(p))
      {
        return p;
      }

      if (is_tau(p))
      {
        return p;
      }

      if (is_sync(p))
      {
        return process::sync(
                 substitute_pCRLproc(process::sync(p).left(),sigma),
                 substitute_pCRLproc(process::sync(p).right(),sigma));
      }

      throw mcrl2::runtime_error("Internal error: expect a pCRL process (2) " + process::pp(p));
      return process_expression(); 
    }


    // The function below transforms a ProcessAssignment to a Process, provided
    // that the process is defined in objectnames.

    process_instance_assignment transform_process_instance_to_process_instance_assignment(
              const process_instance& procId,
              const std::set<variable>& bound_variables=std::set<variable>())
    {
      objectdatatype& object=objectIndex(procId.identifier());
      const variable_list process_parameters=object.parameters;
      const data_expression_list& rhss=procId.actual_parameters();

      assignment_vector new_assignments;
      data_expression_list::const_iterator j=rhss.begin();
      for(variable_list::const_iterator i=process_parameters.begin(); i!=process_parameters.end(); ++i, ++j)
      {
        assert(j!=rhss.end());
        if (*i==*j)
        {
          if (bound_variables.count(*i)>0) // Now *j is a different variable than *i
          {
            new_assignments.emplace_back(*i, *j);
          }
        }
        else
        {
          new_assignments.emplace_back(*i, *j);
        }
      }
      assert(j==rhss.end());

      assert(check_valid_process_instance_assignment(procId.identifier(),assignment_list(new_assignments.begin(),new_assignments.end())));
      process_instance_assignment p(procId.identifier(), assignment_list(new_assignments.begin(),new_assignments.end()));
      return p;
    }

    /********************************************************************/
    /*                                                                  */
    /*   BELOW THE PROCEDURES ARE GIVEN TO TRANSFORM PROCESSES TO       */
    /*   LINEAR PROCESSES.                                              */
    /*                                                                  */
    /*                                                                  */
    /*                                                                  */
    /********************************************************************/

    using variableposition = enum
    {
      first,
      later
    };

    /****************  tovarheadGNF  *********************************/

    variable_list parameters_that_occur_in_body(
      const variable_list& parameters,
      const process_expression& body)
    {
      variable_vector result;
      for(const variable& p: parameters)
      {
        if (occursinpCRLterm(p,body,false))
        {
          result.push_back(p);
        }
      }
      return variable_list(result.begin(),result.end());
    }

    // The variable below is used to count the number of new processes that
    // are made. If this number is very high, it is likely that the regular
    // flag is used, and an unbounded number of new processes are generated.
    // In such a case a warning is printed suggesting to use regular2.

    process_identifier newprocess(
      const variable_list& parameters,
      const process_expression& body,
      const processstatustype ps,
      const bool canterminate,
      const bool containstime)
    {
      assert(canterminatebody(body)==canterminate);
      assert(containstimebody(body)==containstime);

      process_identifier p1;
      if (searchProcDeclaration(parameters,body,ps, canterminate,containstime,p1))
      {
        return p1;  // The process did already exist. No need to make a new one.
      }

      static std::size_t numberOfNewProcesses = 0;
      static std::size_t warningNumber = 25;
      numberOfNewProcesses++;
      if (numberOfNewProcesses == warningNumber)
      {
        mCRL2log(mcrl2::log::warning) << "Generated " << numberOfNewProcesses << " new internal processes.";

        if (options.lin_method==lmRegular)
        {
          mCRL2log(mcrl2::log::warning) << " A possible unbounded loop can be avoided by using `regular2' or `stack' as linearisation method." << std::endl;
        }
        else if (options.lin_method==lmRegular2)
        {
          mCRL2log(mcrl2::log::warning) << " A possible unbounded loop can be avoided by using `stack' as the linearisation method." << std::endl;
        }
        else
        {
          mCRL2log(mcrl2::log::warning) << std::endl;
        }
        warningNumber=warningNumber*5;
      }
      const variable_list parameters1=parameters_that_occur_in_body(parameters, body);
      const core::identifier_string s=fresh_identifier_generator("P");
      const process_identifier p(s, parameters1);
      assert(std::string(p.name()).size()>0);
      insert_process_declaration(
        p,
        parameters1,
        body,
        ps,
        canterminate,
        containstime);
      return p;
    }


    process_expression wraptime(
      const process_expression& body,
      const data_expression& time,
      const variable_list& freevars)
    {
      if (is_choice(body))
      {
        return choice(
                 wraptime(choice(body).left(),time,freevars),
                 wraptime(choice(body).right(),time,freevars));
      }

      if (is_sum(body))
      {
        variable_list sumvars=sum(body).variables();
        process_expression body1=sum(body).operand();

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        body1=substitute_pCRLproc(body1, sigma);
        maintain_variables_in_rhs< mutable_map_substitution<> >  sigma_aux(sigma);
        const data_expression time1=/*data::*/replace_variables_capture_avoiding_alt(time, sigma_aux);
        body1=wraptime(body1,time1,sumvars+freevars);
        return sum(sumvars,body1);
      }

      if (is_stochastic_operator(body))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(body);
        variable_list sumvars=sto.variables();
        process_expression body1=sto.operand();

        maintain_variables_in_rhs<mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        body1=substitute_pCRLproc(body1, sigma);
        const data_expression new_distribution=replace_variables_capture_avoiding_alt(sto.distribution(), sigma);
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma_aux(sigma);
        const data_expression time1=replace_variables_capture_avoiding_alt(time, sigma_aux);
        body1=wraptime(body1,time1,sumvars+freevars);
        return stochastic_operator(sumvars,new_distribution,body1);
      }

      if (is_if_then(body))
      {
        return if_then(if_then(body).condition(),wraptime(if_then(body).then_case(),time,freevars));
      }

      if (is_seq(body))
      {
        return seq(wraptime(seq(body).left(),time,freevars),seq(body).right());
      }

      if (is_at(body))
      {
        /* make a new process */
        const process_identifier newproc=newprocess(freevars,body,pCRL,
                                         canterminatebody(body),containstimebody(body));
        assert(check_valid_process_instance_assignment(newproc,assignment_list()));
        return at(process_instance_assignment(
                    newproc,
                    assignment_list()),  //data::data_expression_list(objectdata[objectIndex(newproc)].parameters)),
                  time);
      }

      if (// (is_process_instance(body))||
          (is_process_instance_assignment(body))||
          (is_sync(body))||
          (is_action(body))||
          (is_tau(body))||
          (is_delta(body)))
      {
        return at(body,time);
      }

      throw mcrl2::runtime_error("Internal error: expect a pCRL process in wraptime " + process::pp(body));
      return process_expression();
    }

    using state = enum
    {
      alt_state,
      sum_state, /* cond,*/
      seq_state,
      name_state,
      multiaction_state
    };

    variable get_fresh_variable(const std::string& s, const sort_expression& sort, const int reuse_index=-1)
    {
      /* If reuse_index is smaller than 0 (-1 is the default value), an unused variable name is returned,
         based on the string s with sort `sort'. If reuse_index is larger or equal to
         0 the reuse_index+1 generated variable is returned. If for a particular reuse_index
         this function is called for the first time, it is guaranteed that the returned
         variable is unique, and not used as variable elsewhere. Upon subsequent calls
         get_fresh_variable will return the same variable for the same s,sort and reuse_triple.
         This feature is added to make it possible to avoid generating too many different variables. */


      if (reuse_index<0)
      {
        variable v(fresh_identifier_generator(s),sort);
        insertvariable(v,true);
        return v;
      }
      else
      {
        static std::map < int , std::map < variable,variable > > generated_variables;
        variable table_index_term(s,sort);
        variable old_variable;
        if (generated_variables[reuse_index].count(table_index_term)>0)
        {
          old_variable=generated_variables[reuse_index][table_index_term];
        }
        else
        {
          /* A new variable must be generated */
          old_variable=get_fresh_variable(s,sort);
          generated_variables[reuse_index][table_index_term]=old_variable;
        }
        return old_variable;
      }
    }

    variable_list make_pars(const sort_expression_list& sortlist)
    {
      /* this function returns a list of variables,
         corresponding to the sorts in sortlist */

      if (sortlist.empty())
      {
        return variable_list();
      }

      const sort_expression& sort=sortlist.front();

      variable_list result=make_pars(sortlist.tail());
      result.push_front(get_fresh_variable("a",sort));
      return result;
    }

    process_expression distributeActionOverConditions(
      const process_expression& act, // This is a multi-action, actually.
      const data_expression& condition,
      const process_expression& restterm,
      const variable_list& freevars,
      const std::set<variable>& variables_bound_in_sum)
    {
      if (is_if_then(restterm))
      {
        /* Here we check whether the process body has the form
           a (c -> x). For state space generation it turns out
           to be beneficial to rewrite this to c-> a x + !c -> a.delta@0, as in
           certain cases this leads to a reduction of the number
           of states. In this code, we recursively check whether
           the action must be distributed over x. This optimisation
           was observed by Yaroslav Usenko, May 2006. Implemented by JFG.
           On industrial examples, it appears to reduce the state space
           with a factor up to 2.
           Before october 2008 this code was wrong, as it transformed
           an expression to c-> a x, ommitting the a.delta@0. */

        const data_expression& c=down_cast<if_then>(restterm).condition();
        const process_expression r=choice(
                                     distributeActionOverConditions(
                                       act,
                                       lazy::and_(condition,c),
                                       if_then(restterm).then_case(),
                                       freevars,variables_bound_in_sum),
                                     distributeActionOverConditions(
                                       act,
                                       lazy::and_(condition,lazy::not_(c)),
                                       delta_at_zero(),
                                       freevars,variables_bound_in_sum));
        return r;
      }
      if (is_if_then_else(restterm))
      {
        /* Here we check whether the process body has the form
           a (c -> x <> y). For state space generation it turns out
           to be beneficial to rewrite this to c-> a x + !c -> a y, as in
           certain cases this leads to a reduction of the number
           of states, despite the duplication of the a action. In this code,
           we recursively check whether the action must be distributed over
           x and y. This optimisation
           was observed by Yaroslav Usenko, May 2006. Implemented by JFG.
           On industrial examples, it appears to reduce the state space
           with a factor up to 2. */

        const data_expression& c=down_cast<if_then_else>(restterm).condition();
        const process_expression r=choice(
                                     distributeActionOverConditions(
                                       act,
                                       lazy::and_(condition,c),
                                       if_then_else(restterm).then_case(),
                                       freevars,variables_bound_in_sum),
                                     distributeActionOverConditions(
                                       act,
                                       lazy::and_(condition,lazy::not_(c)),
                                       if_then_else(restterm).else_case(),
                                       freevars,variables_bound_in_sum));
        return r;
      }
      const process_expression restterm1=bodytovarheadGNF(restterm,seq_state,freevars,later,variables_bound_in_sum);
      return if_then(condition,seq(act,restterm1));
    }


   static
   assignment_list parameters_to_assignment_list(const variable_list& parameters, const std::set<variable>& variables_bound_in_sum)
   {
     assignment_vector result;
     for(const variable& v: parameters)
     {
       if (variables_bound_in_sum.count(v)>0)
       {
         result.emplace_back(v, v); // rhs is another variable than the lhs!!
       }
     }
     return assignment_list(result.begin(),result.end());
   }


    process_expression bodytovarheadGNF(
      const process_expression& body,
      const state s,
      const variable_list& freevars,
      const variableposition v,
      const std::set<variable>& variables_bound_in_sum)
    {
      /* it is assumed that we only receive processes with
         operators alt, seq, sum_state, cond, name, delta, tau, sync, at and stochastic operator in it */

      if (is_choice(body))
      {
        if (alt_state>=s)
        {
          const process_expression body1=bodytovarheadGNF(choice(body).left(),alt_state,freevars,first,variables_bound_in_sum);
          const process_expression body2=bodytovarheadGNF(choice(body).right(),alt_state,freevars,first,variables_bound_in_sum);
          if (isDeltaAtZero(body1))
          {
            return body2;
          }
          if (isDeltaAtZero(body2))
          {
            return body1;
          }
          return choice(body1,body2);
        }
        const process_expression body1=bodytovarheadGNF(body,alt_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body1,pCRL,
                                         canterminatebody(body1),
                                         containstimebody(body1));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));
      }

      if (is_sum(body))
      {
        if (sum_state>=s)
        {
          variable_list sumvars=sum(body).variables();
          process_expression body1=sum(body).operand();

          maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
          alphaconvert(sumvars,sigma,freevars,data_expression_list());
          body1=substitute_pCRLproc(body1,sigma);
          std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
          variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
          body1=bodytovarheadGNF(body1,sum_state,sumvars+freevars,first,variables_bound_in_sum1);
          /* Due to the optimisation below, suggested by Yaroslav Usenko, bodytovarheadGNF(...,sum_state,...)
             can deliver a process of the form c -> x + !c -> y. In this case, the
             sumvars must be distributed over both summands. */
          if (is_choice(body1))
          {
            return choice(sum(sumvars,choice(body1).left()),
                          sum(sumvars,choice(body1).right()));
          }
          return sum(sumvars,body1);
        }
        const process_expression body1=bodytovarheadGNF(body,alt_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body1,pCRL,
                                         canterminatebody(body1),
                                         containstimebody(body1));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));
      }

      if (is_stochastic_operator(body))
      {
        if (seq_state>=s)
        {
          const stochastic_operator& sto=down_cast<const stochastic_operator>(body);
          variable_list sumvars=sto.variables();
          data_expression distribution=sto.distribution();
          process_expression body1=sto.operand();

          maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
          alphaconvert(sumvars,sigma,freevars,data_expression_list());
          distribution=/* data::*/replace_variables_capture_avoiding_alt(distribution,sigma);
          body1=substitute_pCRLproc(body1,sigma);
          std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
          variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
          body1=bodytovarheadGNF(body1,name_state,sumvars+freevars,v,variables_bound_in_sum1);
          /* Due to the optimisation below, suggested by Yaroslav Usenko, bodytovarheadGNF(...,sum_state,...)
             can deliver a process of the form c -> x + !c -> y. In this case, the
             sumvars must be distributed over both summands. */
          return stochastic_operator(sumvars,distribution,body1);
        }
        const process_expression body_=bodytovarheadGNF(body,seq_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body_,pCRL,
                                         canterminatebody(body_),
                                         containstimebody(body_));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));
      }

      if (is_if_then(body))
      {
        const data_expression& condition=down_cast<if_then>(body).condition();
        const process_expression& body1=down_cast<if_then>(body).then_case();

        if (s<=sum_state)
        {
          return if_then(
                   condition,
                   bodytovarheadGNF(body1,seq_state,freevars,first,variables_bound_in_sum));
        }
        const process_expression body2=bodytovarheadGNF(body,alt_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body2,pCRL,
                                         canterminatebody(body2),
                                         containstimebody(body2));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));

      }

      if (is_if_then_else(body))
      {
        const data_expression& condition=down_cast<if_then_else>(body).condition();
        const process_expression& body1=down_cast<if_then_else>(body).then_case();
        const process_expression& body2=down_cast<if_then_else>(body).else_case();

        if (isDeltaAtZero(body1) && isDeltaAtZero(body2))
        {
          return body1;
        }

        if ((s<=sum_state) && ((isDeltaAtZero(body1))||(isDeltaAtZero(body2))))
        {
          if (isDeltaAtZero(body2))
          {
            return if_then(
                     condition,
                     bodytovarheadGNF(body1,seq_state,freevars,first,variables_bound_in_sum));
          }
          /* body1=="Delta@0" */
          {
            return if_then(
                     lazy::not_(condition),
                     bodytovarheadGNF(body2,seq_state,freevars,first,variables_bound_in_sum));
          }
        }
        if (alt_state==s) /* body1!=Delta@0 and body2!=Delta@0 */
        {
          return
            choice(
              if_then(
                condition,
                bodytovarheadGNF(body1,seq_state,freevars,first,variables_bound_in_sum)),
              if_then(
                lazy::not_(condition),
                bodytovarheadGNF(body2,seq_state,freevars,first,variables_bound_in_sum)));
        }
        const process_expression body3=bodytovarheadGNF(body,alt_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body3,pCRL,
                                         canterminatebody(body3),
                                         containstimebody(body3));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));

      }

      if (is_seq(body))
      {
        process_expression body1=seq(body).left();
        process_expression body2=seq(body).right();

        if (s<=seq_state)
        {
          body1=bodytovarheadGNF(body1,name_state,freevars,v,variables_bound_in_sum);
          if (!canterminatebody(body1))
          {
            /* In this case there is no need to investigate body2, as it cannot be reached. */
            return body1;
          }
          if ((is_if_then(body2)) && (s<=sum_state))
          {
            /* Here we check whether the process body has the form
               a (c -> x) + !c -> delta@0. For state space generation it turns out
               to be beneficial to rewrite this to c-> a x, as in
               certain cases this leads to a reduction of the number
               of states. An extra change (24/12/2006) is that the
               conditions are distributed recursively over
               all conditions. The optimisation
               was observed by Yaroslav Usenko, May 2006. Implemented by JFG.
               On industrial examples, it appears to reduce the state space
               with a factor up to 2. Until 1/11/2008 this code was incorrect,
               because the summand a (!c -> delta@0) was not forgotten.*/

            const data_expression& c=down_cast<if_then>(body2).condition();

            const process_expression r= choice(
                                          distributeActionOverConditions(body1,c,if_then(body2).then_case(),freevars,variables_bound_in_sum),
                                          distributeActionOverConditions(body1,lazy::not_(c),delta_at_zero(),freevars,variables_bound_in_sum));
            return r;
          }
          if ((is_if_then_else(body2)) && (s<=sum_state))
          {

            /* Here we check whether the process body has the form
               a (c -> x <> y). For state space generation it turns out
               to be beneficial to rewrite this to c-> a x + !c -> a y, as in
               certain cases this leads to a reduction of the number
               of states, despite the duplication of the a action. An extra
               change (24/12/2006) is that the conditions are distributed recursively over
               all conditions. The optimisation
               was observed by Yaroslav Usenko, May 2006. Implemented by JFG.
               On industrial examples, it appears to reduce the state space
               with a factor up to 2. */


            const data_expression& c=down_cast<if_then_else>(body2).condition();
            const process_expression r= choice(
                                          distributeActionOverConditions(body1,c,if_then_else(body2).then_case(),freevars,variables_bound_in_sum),
                                          distributeActionOverConditions(body1,lazy::not_(c),if_then_else(body2).else_case(),freevars,variables_bound_in_sum));
            return r;
          }
          body2=bodytovarheadGNF(body2,seq_state,freevars,later,variables_bound_in_sum);
          return seq(body1,body2);
        }
        body1=bodytovarheadGNF(body,alt_state,freevars,first,variables_bound_in_sum);
        const process_identifier newproc=newprocess(freevars,body1,pCRL,canterminatebody(body1),
                                         containstimebody(body1));

        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));
      }

      if (is_action(body))
      {
        if ((s==multiaction_state)||(v==first))
        {
          return body;
        }

        bool isnew=false;
        objectdatatype& object=addMultiAction(action(body),isnew);

        if (object.process_representing_action==process_identifier())
        {
          /* this action does not yet have a corresponding process, which
             must be constructed. The resulting process is stored in
             the variable process_representing_action in objectdata. Tempvar below is
             needed as objectdata may be realloced as a side effect
             of newprocess */
          const process_identifier tempvar=newprocess(
                                             object.parameters,
                                             object.processbody,
                                             GNF,true,false);
          object.process_representing_action=tempvar;
        }
        return transform_process_instance_to_process_instance_assignment(
                          process_instance(object.process_representing_action,
                          action(body).arguments()));
      }

      if (is_sync(body))
      {
        bool isnew=false;
        const process_expression& body1=down_cast<process::sync>(body).left();
        const process_expression& body2=down_cast<process::sync>(body).right();
        const action_list ma=linMergeMultiActionListProcess(
                               bodytovarheadGNF(body1,multiaction_state,freevars,v,variables_bound_in_sum),
                               bodytovarheadGNF(body2,multiaction_state,freevars,v,variables_bound_in_sum));

        const process_expression mp=action_list_to_process(ma);
        if ((s==multiaction_state)||(v==first))
        {
          return mp;
        }

        objectdatatype& object=addMultiAction(mp,isnew);

        if (object.process_representing_action==process_identifier())
        {
          /* this action does not yet have a corresponding process, which
             must be constructed. The resulting process is stored in
             the variable process_representing_action in objectdata. Tempvar below is needed
             as objectdata may be realloced as a side effect of newprocess */
          process_identifier tempvar=newprocess(
                                       object.parameters,
                                       object.processbody,
                                       GNF,true,false);
          object.process_representing_action=tempvar;
        }
        return transform_process_instance_to_process_instance_assignment(
                      process_instance(
                           process_identifier(object.process_representing_action),
                           getarguments(ma)));
      }

      if (is_at(body))
      {
        process_expression body1=bodytovarheadGNF(
                                   at(body).operand(),
                                   s,
                                   freevars,
                                   first,variables_bound_in_sum);
        data_expression time=data_expression(at(body).time_stamp());
        /* put the time operator around the first action or process */
        body1=wraptime(body1,time,freevars);
        if (v==first)
        {
          return body1;
        }

        /* make a new process, containing this process */
        const process_identifier newproc=newprocess(freevars,body1,pCRL,
                                         canterminatebody(body1),
                                         containstimebody(body1));
        assert(check_valid_process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum)));
        return process_instance_assignment(newproc,parameters_to_assignment_list(objectIndex(newproc).parameters,variables_bound_in_sum));
      }

      if (is_process_instance_assignment(body))
      {
        // return transform_process_assignment_to_process(body);
        return body;
      }


      if (is_tau(body))
      {
        if (v==first)
        {
          return tau();
        }
        assert(check_valid_process_instance_assignment(tau_process,assignment_list()));
        return process_instance_assignment(tau_process,assignment_list());
      }

      if (is_delta(body))
      {
        if (v==first)
        {
          return body;
        }
        assert(check_valid_process_instance_assignment(delta_process,assignment_list()));
        return process_instance_assignment(delta_process,assignment_list());
      }

      throw mcrl2::runtime_error("unexpected process format in bodytovarheadGNF " + process::pp(body) +".");
      return process_expression();
    }

    void procstovarheadGNF(const std::vector < process_identifier>& procs)
    {
      /* transform the processes in procs into pre-Greibach Normal Form */
      for (const process_identifier& i: procs)
      {
        objectdatatype& object=objectIndex(i);

        // The intermediate variable result is needed here
        // because objectdata can be reallocated as a side
        // effect of bodytovarheadGNF.

        std::set<variable> variables_bound_in_sum;
        const process_expression result=
          bodytovarheadGNF(
            object.processbody,
            alt_state,
            object.parameters,
            first,
            variables_bound_in_sum);
        object.processbody=result;
      }
    }

    /**************** towards real GREIBACH normal form **************/

    using terminationstatus = enum
    {
      terminating,
      infinite
    };

    process_expression putbehind(const process_expression& body1, const process_expression& body2)
    {
      if (is_choice(body1))
      {
        return choice(
                 putbehind(choice(body1).left(),body2),
                 putbehind(choice(body1).right(),body2));
      }

      if (is_seq(body1))
      {
        return seq(seq(body1).left(), putbehind(seq(body1).right(),body2));
      }

      if (is_if_then(body1))
      {
        return if_then(if_then(body1).condition(),putbehind(if_then(body1).then_case(),body2));
      }

      if (is_sum(body1))
      {
        /* we must take care that no variables in body2 are
            inadvertently bound */
        variable_list sumvars=sum(body1).variables();

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvertprocess(sumvars,sigma,body2);
        return sum(sumvars,
                   putbehind(substitute_pCRLproc(sum(body1).operand(), sigma),
                             body2));
      }

      if (is_stochastic_operator(body1))
      {
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(body1);
        variable_list stochvars=sto.variables();

        // See explanation at sum operator above.
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvertprocess(stochvars,sigma,body2);
        return stochastic_operator(
                 stochvars,
                 sto.distribution(),
                 putbehind(
                       substitute_pCRLproc(sto.operand(),sigma),
                       body2));
      }


      if (is_action(body1))
      {
        return seq(body1,body2);
      }

      if (is_sync(body1))
      {
        return seq(body1,body2);
      }

      if (is_process_instance_assignment(body1))
      {
        return seq(body1,body2);
      }

      if (is_delta(body1))
      {
        return body1;
      }

      if (is_tau(body1))
      {
        return seq(body1,body2);
        // throw mcrl2::runtime_error("Expect only multiactions, not a tau.");
      }

      if (is_at(body1))
      {
        return seq(body1,body2);
      }

      throw mcrl2::runtime_error("Internal error. Unexpected process format in putbehind " + process::pp(body1) +".");
      return process_expression();
    }


    process_expression distribute_condition(
      const process_expression& body1,
      const data_expression& condition)
    {
      if (is_choice(body1))
      {
        return choice(
                 distribute_condition(choice(body1).left(),condition),
                 distribute_condition(choice(body1).right(),condition));
      }

      if (is_seq(body1))
      {
        return if_then(condition,body1);
      }

      if (is_if_then(body1))
      {
        return if_then(
                 lazy::and_(if_then(body1).condition(),condition),
                 if_then(body1).then_case());
      }

      if (is_sum(body1))
      {
        /* we must take care that no variables in condition are
            inadvertently bound */
        variable_list sumvars=sum(body1).variables();
        maintain_variables_in_rhs< mutable_map_substitution<> >  sigma;
        alphaconvert(sumvars,sigma,variable_list(), data_expression_list({ condition }));
        return sum(
                 sumvars,
                 distribute_condition(
                   substitute_pCRLproc(sum(body1).operand(),sigma),
                   condition));
      }

      if (is_stochastic_operator(body1))
      {
        /* we must take care that no variables in condition are
            inadvertently bound */
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(body1);
        variable_list stochvars=sto.variables();
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(stochvars,sigma,variable_list(), data_expression_list({ condition }));

        return stochastic_operator(
                 stochvars,
                 if_(condition,
                        replace_variables_capture_avoiding_alt(
                                               sto.distribution(),
                                               sigma),
                        if_(variables_are_equal_to_default_values(stochvars),real_one(),real_zero())),
                 distribute_condition(
                       substitute_pCRLproc(sto.operand(),sigma),
                       condition));
      }

      if (is_at(body1)||
          is_action(body1)||
          is_sync(body1)||
          is_process_instance_assignment(body1)||
          is_delta(body1)||
          is_tau(body1))
      {
        return if_then(condition,body1);
      }

      throw mcrl2::runtime_error("Internal error. Unexpected process format in distribute condition " + process::pp(body1) +".");
    }

    /* This process calculates the equivalent of sum sumvars dist stochvars[distribution] body
       where the distribution occurs in front. This can only be done under limited circumstances.
       Assume we can enumerate the values of sumvars by e1,...,en. Introduce n copies of stochvars,
       i.e., stochvars1,...,stochvarsn. The new process becomes

           dist stochvars1,stochvars2,...,stochvarsn
                [distribution(e1,stochvars1)*distribution(e2,stochvars2)*...*distribution(en,stochvarsn)].

           body(e1,stochvars1)+
           body(e2,stochvars2)+
           ...
           body(en,stochvarsn)
    */

    process_expression enumerate_distribution_and_sums(
                         const variable_list& sumvars,
                         const variable_list& stochvars,
                         const data_expression& distribution,
                         const process_expression& body)
    {
      if (options.norewrite)
      {
        throw mcrl2::runtime_error("The use of the rewriter must be allowed to distribute a sum operator over a distribution.");
      }


      std::vector < data_expression_vector > data_vector(1,data_expression_vector());
      for(const variable& v:sumvars)
      {
        std::vector < data_expression_vector > new_data_vector;

        if (!data.is_certainly_finite(v.sort()))
        {
          throw mcrl2::runtime_error("Cannot distribute a sum variable of non finite sort " + pp(v.sort()) + " over a distribution, which is required for linearisation.");
        }

        for(const data_expression& d: enumerate_expressions(v.sort(),data,rewr))
        {
          for(const data_expression_vector& dv: data_vector)
          {
            data_expression_vector new_dv=dv;
            new_dv.push_back(d);
            new_data_vector.push_back(new_dv);
          }

        }
        data_vector.swap(new_data_vector);

      }
      assert(!data_vector.empty());

      process_expression resulting_body;
      data_expression resulting_distribution;
      variable_list resulting_stochastic_variables;
      bool result_defined=false;
      for(const data_expression_vector& d: data_vector)
      {
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        variable_list vl = stochvars;
        alphaconvert(vl,sigma,stochvars,data_expression_list());
        data_expression_vector::const_iterator i=d.begin();
        for(const variable& v: sumvars)
        {
          sigma[v] = *i;
          ++i;
        }
        const process_expression d1=substitute_pCRLproc(body,sigma);
        const data_expression new_distribution=replace_variables_capture_avoiding_alt(distribution, sigma);

        if (result_defined)
        {
          resulting_body=choice(resulting_body, d1);
          resulting_distribution=sort_real::times(resulting_distribution,new_distribution);
          resulting_stochastic_variables=resulting_stochastic_variables + vl;
        }
        else
        {
          resulting_body=d1;
          resulting_distribution=new_distribution;
          resulting_stochastic_variables=vl;
          result_defined=true;
        }
      }
      /* Put the distribution in front. */

      return stochastic_operator(resulting_stochastic_variables, resulting_distribution, resulting_body);
    }

    process_expression distribute_sum_over_a_stochastic_operator(
                         const variable_list& sumvars,
                         const variable_list& stochastic_variables,
                         const data_expression& distribution,
                         const process_expression& body)
    {
      if (is_sum(body)||
          is_choice(body)||
          is_seq(body)||
          is_if_then(body)||
          is_sync(body)||
          is_action(body)||
          is_tau(body)||
          is_at(body)||
          is_process_instance_assignment(body)||
          isDeltaAtZero(body))
      {
        return enumerate_distribution_and_sums(sumvars,stochastic_variables,distribution,body);
      }

      if (is_delta(body)||
          is_tau(body))
      {
        return body;
      }

      if (is_stochastic_operator(body))
      {
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(body);
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        variable_list inner_stoch_vars=sto.variables();
        alphaconvert(inner_stoch_vars,sigma,sumvars,data_expression_list());
        const process_expression new_body=substitute_pCRLproc(sto.operand(), sigma);
        const data_expression new_distribution=replace_variables_capture_avoiding_alt(sto.distribution(), sigma);
        return distribute_sum_over_a_stochastic_operator(sumvars,
                                                         stochastic_variables + inner_stoch_vars,
                                                         sort_real::times(distribution,new_distribution), new_body);
      }

      throw mcrl2::runtime_error("Internal error. Unexpected process format in distribute_sum " + process::pp(body) +".");
      return process_expression();
    }

    process_expression distribute_sum(
      const variable_list& sumvars,
      const process_expression& body1)
    {
      if (is_choice(body1))
      {
        return choice(
                 distribute_sum(sumvars,choice(body1).left()),
                 distribute_sum(sumvars,choice(body1).right()));
      }

      if (is_seq(body1)||
          is_if_then(body1)||
          is_sync(body1)||
          is_action(body1)||
          is_tau(body1)||
          is_at(body1)||
          is_process_instance_assignment(body1)||
          isDeltaAtZero(body1))
      {
        return sum(sumvars,body1);
      }

      if (is_sum(body1))
      {
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        variable_list inner_sumvars=sum(body1).variables();
        alphaconvert(inner_sumvars,sigma,sumvars,data_expression_list());
        const process_expression new_body1=substitute_pCRLproc(sum(body1).operand(), sigma);
        return sum(sumvars+inner_sumvars,new_body1);
      }

      if (is_stochastic_operator(body1))
      {
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(body1);
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        variable_list inner_stoch_vars=sto.variables();
        alphaconvert(inner_stoch_vars,sigma,sumvars,data_expression_list());
        const process_expression new_body1=substitute_pCRLproc(sto.operand(), sigma);
        const data_expression new_distribution=replace_variables_capture_avoiding_alt(sto.distribution(), sigma);
        return distribute_sum_over_a_stochastic_operator(sumvars, inner_stoch_vars, new_distribution, new_body1);
      }

      if (is_delta(body1)||
          is_tau(body1))
      {
        return body1;
      }

      throw mcrl2::runtime_error("Internal error. Unexpected process format in distribute_sum " + process::pp(body1) +".");
      return process_expression();
    }

    static
    int match_sequence(
      const std::vector < process_instance_assignment >& s1,
      const std::vector < process_instance_assignment >& s2,
      const bool regular2)
    {
      /* s1 and s2 are sequences of typed variables of
         the form Process(ProcVarId("P2",[SortId("Bit"),
         SortId("Bit")]),[OpId("b1",SortId("Bit")),OpId("b1",SortId("Bit"))]).
         This function yields true if the names and types of
         the processes in s1 and s2 match. */

      std::vector < process_instance_assignment >::const_iterator i2=s2.begin();
      for (std::vector < process_instance_assignment >::const_iterator i1=s1.begin();
           i1!=s1.end(); ++i1,++i2)
      {
        if (i2==s2.end())
        {
          return false;
        }
        if (regular2)
        {
          if (i1->identifier()!=i2->identifier())
          {
            return false;
          }
        }
        else
        {
          if (*i1!=*i2)
          {
            return false;
          }
        }
      }
      if (i2!=s2.end())
      {
        return false;
      }
      return true;
    }

    bool exists_variable_for_sequence(
      const std::vector < process_instance_assignment >& process_names,
      process_identifier& result)
    {
      std::vector < std::vector < process_instance_assignment > >::const_iterator rwalker=representedprocesses.begin();
      for (std::vector < process_identifier >::const_iterator walker=seq_varnames.begin();
           walker!=seq_varnames.end(); ++walker,++rwalker)
      {
        assert(rwalker!=representedprocesses.end());
        const process_identifier process=*walker;
        if (match_sequence(process_names,*rwalker,options.lin_method==lmRegular2))
        {
          result=process;
          return true;
        }
      }
      assert(rwalker==representedprocesses.end());
      return false;
    }

    void extract_names(
      const process_expression& sequence,
      std::vector < process_instance_assignment >& result)
    {
      if (is_action(sequence)||is_process_instance_assignment(sequence))
      {
        result.push_back(atermpp::down_cast<process_instance_assignment>(sequence));
        return;
      }

      if (is_stochastic_operator(sequence))
      {
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(sequence);
        extract_names(sto.operand(),result);
        return;
      }

      if (is_seq(sequence))
      {
        const process_expression& first=down_cast<seq>(sequence).left();
        if (is_process_instance_assignment(first))
        {
          result.push_back(atermpp::down_cast<process_instance_assignment>(first));
          objectdatatype& object=objectIndex(atermpp::down_cast<process_instance_assignment>(first).identifier());
          if (object.canterminate)
          {
            extract_names(seq(sequence).right(),result);
          }
          return;
        }
      }

      throw mcrl2::runtime_error("Internal error. Expected sequence of process names (1) " + process::pp(sequence) + ".");
    }

    variable_list parscollect(const process_expression& oldbody, process_expression& newbody)
    {
      /* we expect that oldbody is a sequence of process instances */

      if (is_process_instance_assignment(oldbody))
      {
        const process_identifier procId=process_instance_assignment(oldbody).identifier();
        const variable_list parameters=objectIndex(procId).parameters;
        assert(check_valid_process_instance_assignment(procId,data::assignment_list()));
        newbody=process_instance_assignment(procId,data::assignment_list());
        return parameters;
      }

      if (is_seq(oldbody))
      {
        const process_expression& first=down_cast<seq>(oldbody).left();
        if (is_process_instance_assignment(first))
        {
          objectdatatype& object=objectIndex(process_instance_assignment(first).identifier());
          if (object.canterminate)
          {
            const process_identifier procId=process_instance_assignment(first).identifier();
            const variable_list pars=parscollect(seq(oldbody).right(),newbody);
            variable_list pars1;
            variable_list pars2;
            const variable_list new_pars=construct_renaming(pars,objectIndex(procId).parameters,pars1,pars2,false);
            assignment_vector new_assignment;
            for(variable_list::const_iterator i=pars2.begin(), j=new_pars.begin(); i!=pars2.end(); ++i,++j)
            {
              assert(j!=new_pars.end());
              new_assignment.emplace_back(*i, *j);
            }
            assert(check_valid_process_instance_assignment(procId,assignment_list(new_assignment.begin(),new_assignment.end())));
            newbody=seq(process_instance_assignment(procId,assignment_list(new_assignment.begin(),new_assignment.end())),newbody);
            const variable_list result=pars1+pars;
            assert(std::set<variable>(result.begin(),result.end()).size()==result.size()); // all elements in the result are unique.
            return result;
          }
          else
          {
            return parscollect(first,newbody);
          }
        }
      }

      throw mcrl2::runtime_error("Internal error. Expected a sequence of process names (2) " + process::pp(oldbody) +".");
      return variable_list();
    }

    assignment_list argscollect_regular(
                         const process_expression& t,
                         const variable_list& vl,
                         const std::set<variable>& variables_bound_in_sum)
    {
      assignment_vector result;
      for(const variable& v: vl)
      {
        if (variables_bound_in_sum.count(v)>0 && occursinpCRLterm(v,t,false))
        {
          result.emplace_back(v, v); // Here an identity assignment is used, as it is possible
                                     // that the *i at the lhs is not the same variable as *i at the rhs.
        }
      }
      return assignment_list(result.begin(),result.end());
    }

    assignment_list argscollect_regular2(const process_expression& t, variable_list& vl)
    {
      if (is_process_instance_assignment(t))
      {
        const process_instance_assignment p(t);
        objectdatatype& object=objectIndex(p.identifier());

        const variable_list pars=object.parameters; // These are the old parameters of the process.
        assert(pars.size()<=vl.size());

        std::map<variable,data_expression>sigma;
        for(const assignment& a: p.assignments())
        {
          sigma[a.lhs()]=a.rhs();
        }

        assignment_list result;
        for(variable_list::const_iterator i=pars.begin(); i!=pars.end(); ++i, vl.pop_front())
        {
          assert(!vl.empty());
          const data_expression new_rhs=make_map_substitution(sigma)(*i);
          result.push_front(assignment(vl.front(),new_rhs)); // Identity assignments are stored as lhs and rhs may
                                                             // refer to different variables.
        }
        return reverse(result);
      }

      if (is_seq(t))
      {
        const process_instance_assignment firstproc=atermpp::down_cast<process_instance_assignment>(seq(t).left());
        objectdatatype& object=objectIndex(firstproc.identifier());
        const assignment_list first_assignment=argscollect_regular2(firstproc,vl);
        if (object.canterminate)
        {
          return first_assignment + argscollect_regular2(seq(t).right(),vl);
        }
        return first_assignment;
      }

      throw mcrl2::runtime_error("Internal error. Expected a sequence of process names (3) " + process::pp(t) +".");
    }

    process_expression cut_off_unreachable_tail(const process_expression& t)
    {
      if (is_process_instance_assignment(t)||is_delta(t)||is_action(t)||is_tau(t)||is_sync(t))
      {
        return t;
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(t);
        return stochastic_operator(sto.variables(),sto.distribution(), cut_off_unreachable_tail(sto.operand()));
      }

      if (is_seq(t))
      {
        const process_expression& firstproc=down_cast<seq>(t).left();
        objectdatatype& object=objectIndex(process_instance_assignment(firstproc).identifier());
        if (object.canterminate)
        {
          return seq(firstproc,cut_off_unreachable_tail(seq(t).right()));
        }
        return firstproc;
      }

      throw mcrl2::runtime_error("Internal error. Expected a sequence of process names (4) " + process::pp(t) +".");
      return process_expression();
    }

    process_expression create_regular_invocation(
      process_expression sequence,
      std::vector <process_identifier>& todo,
      const variable_list& freevars,
      const std::set<variable>& variables_bound_in_sum)
    {
      process_identifier new_process;

      /* Sequence consists of a sequence of process references,
         concatenated with the sequential composition operator, which may be preceded by
         a stochastic operator.  */
      if (is_stochastic_operator(sequence))
      {
        /* Leave the stochastic operators in place */
        const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(sequence);
        std::set<variable> variables_bound_in_sum_new=variables_bound_in_sum;
        variables_bound_in_sum_new.insert(sto.variables().begin(),sto.variables().end());
        return stochastic_operator(sto.variables(),
                                   sto.distribution(),
                                   create_regular_invocation(sto.operand(),todo,freevars+sto.variables(),variables_bound_in_sum_new));
      }

      sequence=cut_off_unreachable_tail(sequence);
      sequence=pCRLrewrite(sequence);
      std::vector < process_instance_assignment > process_names;
      extract_names(sequence,process_names);
      assert(!process_names.empty());

      if (process_names.size()==1)
      {
        /* length of list equals 1 */
        if (is_process_instance_assignment(sequence))
        {
          return sequence;
        }
        if (is_seq(sequence))
        {
          return seq(sequence).left();
        }
        throw mcrl2::runtime_error("Internal error. Expected a sequence of process names " + process::pp(sequence) +".");
      }
      /* There is more than one process name in the sequence,
         so, we must replace them by a single name */

      /* We first start out by searching whether
         there is already a variable with a matching sequence
         of variables */
      if (!exists_variable_for_sequence(process_names,new_process))
      {
        /* There does not exist an appropriate variable,
           so, make it and return its index in n */
        process_expression newbody;
        if (options.lin_method==lmRegular2)
        {
          variable_list pars=parscollect(sequence,newbody);
          new_process=newprocess(pars,newbody,pCRL,
                                 canterminatebody(newbody),
                                 containstimebody(newbody));
          representedprocesses.push_back(process_names);
        }
        else
        {
          new_process=newprocess(freevars,sequence,pCRL,
                                 canterminatebody(sequence),containstimebody(sequence));
          representedprocesses.push_back(process_names);
        }
        seq_varnames.push_back(new_process);
        todo.push_back(new_process);
      }
      /* now we must construct arguments */
      variable_list parameters=objectIndex(new_process).parameters;
      if (options.lin_method==lmRegular2)
      {
        const assignment_list args=argscollect_regular2(sequence,parameters);
        assert(check_valid_process_instance_assignment(new_process,args));
        const process_expression p=process_instance_assignment(new_process,args);
        return p;
      }
      else
      {
        assert(check_valid_process_instance_assignment(new_process,argscollect_regular(sequence,parameters,variables_bound_in_sum)));
        return process_instance_assignment(new_process,argscollect_regular(sequence,parameters,variables_bound_in_sum));
      }
    }

    process_expression to_regular_form(
      const process_expression& t,
      std::vector <process_identifier>& todo,
      const variable_list& freevars,
      const std::set<variable>& variables_bound_in_sum)
    /* t has the form of the sum, and condition over actions
       each followed by a sequence of variables. We replace
       this variable by a single one, putting the new variable
       on the todo list, to be transformed to regular form also. */
    {
      if (is_choice(t))
      {
        const process_expression t1=to_regular_form(choice(t).left(),todo,freevars,variables_bound_in_sum);
        const process_expression t2=to_regular_form(choice(t).right(),todo,freevars,variables_bound_in_sum);
        return choice(t1,t2);
      }

      if (is_seq(t))
      {
        const process_expression& firstact=down_cast<seq>(t).left();
        assert(is_at(firstact)||is_tau(firstact)||is_action(firstact)||is_sync(firstact));
        /* the sequence of variables in
                   the second argument must be replaced */
        return seq(firstact,create_regular_invocation(seq(t).right(),todo,freevars,variables_bound_in_sum));
      }

      if (is_if_then(t))
      {
        return if_then(if_then(t).condition(),to_regular_form(if_then(t).then_case(),todo,freevars,variables_bound_in_sum));
      }

      if (is_sum(t))
      {
        variable_list sumvars=sum(t).variables();

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        const process_expression body=substitute_pCRLproc(sum(t).operand(), sigma);

        std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
        variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
        return sum(sumvars,
                   to_regular_form(
                     body,
                     todo,
                     sumvars+freevars,
                     variables_bound_in_sum1));
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(t);
        variable_list sumvars=sto.variables();

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        const data_expression distribution=replace_variables_capture_avoiding_alt(
                                               sto.distribution(),
                                               sigma);
        const process_expression body=substitute_pCRLproc(sto.operand(),sigma);

        std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
        variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
        return stochastic_operator(
                   sumvars,
                   distribution,
                   to_regular_form(
                     body,
                     todo,
                     sumvars+freevars,
                     variables_bound_in_sum1));
      }

      if (is_sync(t)||is_action(t)||is_delta(t)||is_tau(t)||is_at(t))
      {
        return t;
      }

      throw mcrl2::runtime_error("to regular form expects GNF " + process::pp(t) +".");
      return process_expression();
    }

    process_expression distributeTime(
      const process_expression& body,
      const data_expression& time,
      const variable_list& freevars,
      data_expression& timecondition)
    {
      if (is_choice(body))
      {
        return choice(
                 distributeTime(choice(body).left(),time,freevars,timecondition),
                 distributeTime(choice(body).right(),time,freevars,timecondition));
      }

      if (is_sum(body))
      {
        variable_list sumvars=sum(body).variables();
        process_expression body1=sum(body).operand();
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        body1=substitute_pCRLproc(body1, sigma);
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma_aux(sigma);
        const data_expression time1=replace_variables_capture_avoiding_alt(time,sigma_aux);
        body1=distributeTime(body1,time1,sumvars+freevars,timecondition);
        return sum(sumvars,body1);
      }

      if (is_stochastic_operator(body))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(body);
        variable_list sumvars=sto.variables();
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,freevars,data_expression_list());
        process_expression new_body=substitute_pCRLproc(sto.operand(), sigma);
        data_expression new_distribution=/*data::*/replace_variables_capture_avoiding_alt(sto.distribution(), sigma);
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma_aux(sigma);
        const data_expression time1=/*data::*/replace_variables_capture_avoiding_alt(time,sigma_aux);
        new_body=distributeTime(new_body,time1,sumvars+freevars,timecondition);
        return stochastic_operator(sumvars,new_distribution,new_body);
      }

      if (is_if_then(body))
      {
        data_expression timecondition=sort_bool::true_();
        process_expression body1=distributeTime(
                                   if_then(body).then_case(),
                                   time,
                                   freevars,
                                   timecondition);

        return if_then(
                 lazy::and_(data_expression(if_then(body).condition()),timecondition),
                 body1);
      }

      if (is_seq(body))
      {
        return seq(
                 distributeTime(seq(body).left(), time,freevars,timecondition),
                 seq(body).right());
      }

      if (is_at(body))
      {
        /* make a new process */
        timecondition=equal_to(time,data_expression(at(body).time_stamp()));
        return body;
      }

      if ((is_sync(body))||
          (is_action(body))||
          (is_tau(body))||
          (is_delta(body)))
      {
        return at(body,time);
      }

      throw mcrl2::runtime_error("Internal error: expect a pCRL process in distributeTime " + process::pp(body) +".");
      return process_expression();
    }

    process_expression procstorealGNFbody(
      const process_expression& body,
      variableposition v,
      std::vector <process_identifier>& todo,
      const bool regular,
      processstatustype mode,
      const variable_list& freevars,
      const std::set <variable>& variables_bound_in_sum)
    /* This process delivers the transformation of body
       to GNF with actions as a head symbol, or it
       delivers NULL if body is not a pCRL process.
       If regular=1, then an attempt is made to obtain a
       GNF where one action is always followed by a
       variable. */
    {
      if (is_at(body))
      {
        data_expression timecondition=sort_bool::true_();
        process_expression body1=procstorealGNFbody(
                                   at(body).operand(),
                                   first,
                                   todo,
                                   regular,
                                   mode,
                                   freevars,
                                   variables_bound_in_sum);
        return distributeTime(
                 body1,
                 at(body).time_stamp(),
                 freevars,
                 timecondition);
      }

      if (is_choice(body))
      {
        const process_expression body1=procstorealGNFbody(choice(body).left(),first,todo,
                                       regular,mode,freevars,variables_bound_in_sum);
        const process_expression body2=procstorealGNFbody(choice(body).right(),first,todo,
                                       regular,mode,freevars,variables_bound_in_sum);
        return choice(body1,body2);
      }

      if (is_seq(body))
      {
        const process_expression body1=procstorealGNFbody(seq(body).left(),v,
                                       todo,regular,mode,freevars,variables_bound_in_sum);
        const process_expression body2=procstorealGNFbody(seq(body).right(),later,
                                       todo,regular,mode,freevars,variables_bound_in_sum);
        process_expression t3=putbehind(body1,body2);
        if ((regular) && (v==first))
        {
          /* We must transform t3 to regular form */
          t3=to_regular_form(t3,todo,freevars,variables_bound_in_sum);
        }
        return t3;
      }

      if (is_if_then(body))
      {
        const process_expression r=distribute_condition(
                                     procstorealGNFbody(if_then(body).then_case(),first,
                                         todo,regular,mode,freevars,variables_bound_in_sum),
                                     if_then(body).condition());
        return r;
      }

      if (is_sum(body))
      {
        const variable_list sumvars=sum(body).variables();
        std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
        variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
        return distribute_sum(sumvars,
                              procstorealGNFbody(sum(body).operand(),first,
                                  todo,regular,mode,sumvars+freevars,variables_bound_in_sum1));
      }

      if (is_stochastic_operator(body))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(body);
        const variable_list& sumvars=sto.variables();
        std::set<variable> variables_bound_in_sum1=variables_bound_in_sum;
        variables_bound_in_sum1.insert(sumvars.begin(),sumvars.end());
        return stochastic_operator(
                              sumvars,
                              sto.distribution(),
                              procstorealGNFbody(sto.operand(),v,
                                  todo,regular,mode,sumvars+freevars,variables_bound_in_sum1));
      }

      if (is_action(body))
      {
        return body;
      }

      if (is_sync(body))
      {
        return body;
      }

      if (is_process_instance_assignment(body))
      {
        process_identifier t=process_instance_assignment(body).identifier();

        if (v==later)
        {
          if (regular)
          {
            mode=mCRL;
          }
          todo.push_back(t);
          /* if ((!regular)||(mode=mCRL))
              todo.push_back(t);
                / * single = in `mode=mCRL' is important, otherwise crash
                   I do not understand the reason for this at this moment
                   JFG (9/5/2000) */
          return body;
        }

        objectdatatype& object=objectIndex(t);
        if (object.processstatus==mCRL)
        {
          todo.push_back(t);
          return process_expression();
        }
        /* The variable is a pCRL process and v==first, so,
           we must now substitute */
        procstorealGNFrec(t,first,todo,regular);

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        for(const assignment& a: static_cast<const process_instance_assignment&>(body).assignments())
        {
          sigma[a.lhs()]=a.rhs();
        }
        process_expression t3=substitute_pCRLproc(object.processbody,sigma);
        if (regular)
        {
          t3=to_regular_form(t3,todo,freevars,variables_bound_in_sum);
        }

        return t3;
      }

      if (is_delta(body))
      {
        return body;
      }

      if (is_tau(body))
      {
        return body;
      }

      if (is_merge(body))
      {
        procstorealGNFbody(process::merge(body).left(),later,
                           todo,regular,mode,freevars,variables_bound_in_sum);
        procstorealGNFbody(process::merge(body).right(),later,
                           todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      if (is_hide(body))
      {
        procstorealGNFbody(hide(body).operand(),later,todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      if (is_rename(body))
      {
        procstorealGNFbody(process::rename(body).operand(),later,todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      if (is_allow(body))
      {
        procstorealGNFbody(allow(body).operand(),later,todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      if (is_block(body))
      {
        procstorealGNFbody(block(body).operand(),later,todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      if (is_comm(body))
      {
        procstorealGNFbody(comm(body).operand(),later,todo,regular,mode,freevars,variables_bound_in_sum);
        return process_expression();
      }

      // The exception below can especially trigger on a process_instance as these ought to have
      // been removed earlier in favour of process_instance_assignments.
      throw mcrl2::runtime_error("unexpected process format in procstorealGNF " + process::pp(body) +".");
      return process_expression();
    }


    void procstorealGNFrec(
      const process_identifier& procIdDecl,
      const variableposition v,
      std::vector <process_identifier>& todo,
      const bool regular)

    /* Do a depth first search on process variables and substitute
       for the headvariable of a pCRL process, in case it is a process,
       such that we obtain a Greibach Normal Form. All pCRL processes will
       be labelled with GNF to indicate that they are in
       Greibach Normal Form. */

    {
      objectdatatype& object=objectIndex(procIdDecl);
      if (object.processstatus==pCRL)
      {
        object.processstatus=GNFbusy;
        std::set<variable> variables_bound_in_sum;
        const process_expression t=procstorealGNFbody(object.processbody,first,
                                   todo,regular,pCRL,object.parameters,variables_bound_in_sum);
        if (object.processstatus!=GNFbusy)
        {
          throw mcrl2::runtime_error("There is something wrong with recursion.");
        }

        object.processbody=t;
        object.processstatus=GNF;
        return;
      }

      if (object.processstatus==mCRL)
      {
        object.processstatus=mCRLbusy;
        std::set<variable> variables_bound_in_sum;
        const process_expression t=procstorealGNFbody(object.processbody,first,todo,
                                   regular,mCRL,object.parameters,variables_bound_in_sum);
        /* if the last result is not equal to NULL,
           the body of this process is itself a processidentifier */

        object.processstatus=mCRLdone;
        return;
      }

      if ((object.processstatus==GNFbusy) && (v==first))
      {
        throw mcrl2::runtime_error("Unguarded recursion in process " + process::pp(procIdDecl) +".");
      }

      if ((object.processstatus==GNFbusy)||
          (object.processstatus==GNF)||
          (object.processstatus==mCRLdone)||
          (object.processstatus==multiAction))
      {
        return;
      }

      if (object.processstatus==mCRLbusy)
      {
        throw mcrl2::runtime_error("Unguarded recursion in process " + process::pp(procIdDecl) +".");
      }

      throw mcrl2::runtime_error("strange process type: " + std::to_string(object.processstatus));
    }

    void procstorealGNF(const process_identifier& procsIdDecl,
                        const bool regular)
    {
      std::vector <process_identifier> todo;
      todo.push_back(procsIdDecl);
      while (!todo.empty())
      {
        const process_identifier pi=todo.back();
        todo.pop_back();
        procstorealGNFrec(pi,first,todo,regular);
      }
    }


    /**************** GENERaTE LPE **********************************/
    /*                                                              */
    /*                                                              */
    /*                                                              */
    /*                                                              */
    /*                                                              */
    /*                                                              */
    /****************************************************************/


    /**************** Make pCRL procs  ******************************/

    /* make_pCRL_procs searches for all process identifiers reachable from id in pCRL equations. */

    void make_pCRL_procs(const process_identifier& id,
                         std::set<process_identifier>& reachable_process_identifiers)
    {
      if (reachable_process_identifiers.count(id)==0)  // not found
      {
        reachable_process_identifiers.insert(id);
        make_pCRL_procs(objectIndex(id).processbody,reachable_process_identifiers);
      }
      return;
    }

    void make_pCRL_procs(const process_expression& t,
                         std::set<process_identifier>& reachable_process_identifiers)
    {
      if (is_choice(t))
      {
        make_pCRL_procs(choice(t).left(),reachable_process_identifiers);
        make_pCRL_procs(choice(t).right(),reachable_process_identifiers);
        return;
      }

      if (is_seq(t))
      {
        make_pCRL_procs(seq(t).left(),reachable_process_identifiers);
        make_pCRL_procs(seq(t).right(),reachable_process_identifiers);
        return;
      }

      if (is_if_then(t))
      {
        make_pCRL_procs(if_then(t).then_case(),reachable_process_identifiers);
        return;
      }

      if (is_sum(t))
      {
        make_pCRL_procs(sum(t).operand(),reachable_process_identifiers);
        return;
      }

      if (is_stochastic_operator(t))
      {
        make_pCRL_procs(stochastic_operator(t).operand(),reachable_process_identifiers);
        return;
      }

      if (is_process_instance_assignment(t))
      {
        make_pCRL_procs(atermpp::down_cast<process_instance_assignment>(t).identifier(), reachable_process_identifiers);
        return;
      }

      if (is_sync(t)||is_action(t)||is_tau(t)||is_delta(t)||is_at(t))
      {
        return;
      }

      throw mcrl2::runtime_error("unexpected process format " + process::pp(t) + " in make_pCRL_procs");
    }

    /**************** minimize_set_of_reachable_process_identifiers  ******************************/

    /* The process identifiers in reachable_process_identifiers have bodies and parameters lists
     * that can be the same. If that is the case, they are identified. The right hand sides are substituted to
     * reflect the identified processes. Example: P1 = a.P2, P3=a.P4, P2=b.Q, P4=b.Q. This reduces to P1 = a.P2, P2 = b.Q.
     * The result is a map from remaining process identifiers to their adapted process bodies.
     * The initial process was part of the reachable_process_identifiers and must be part of the map.
     */

    /* set_proc_identifier_map is an auxiliary procedure that effectively sets identfier_identifier_map[id1]:=... or
       identifier_identifier_map[id1]:=... in such a way that they map to the same element. In this way identifier_identifier_map
       acts as a representation of the equivalence set induced by equations id1==id2. A constraint is that
       initial_process may never occur as a left hand side. An invariant of identifier_identifier_map is
       that no id1 that occurs at the right of the map can occur at the left.
    */

    static
    process_identifier get_last(const process_identifier& id, const std::map< process_identifier, process_identifier >& identifier_identifier_map)
    {
      process_identifier target_id=id;
      bool ready=false;
      do
      {
        const std::map< process_identifier, process_identifier >::const_iterator i=identifier_identifier_map.find(target_id);
        if (i==identifier_identifier_map.end())
        {
          ready=true;
        }
        else
        {
          target_id=i->second;
        }
      }
      while (!ready);
      return target_id;
    }


    static
    void set_proc_identifier_map(
            std::map< process_identifier, process_identifier >& identifier_identifier_map,
            const process_identifier& id1_,
            const process_identifier& id2_,
            const process_identifier& initial_process)
    {
      assert(id1_!=id2_);

      /* Take care that id1 is the last identifier or that id2 is arger than id1. This guarantees
         that there will be no loops in the mapping of identifiers. */
      // make_substitution sigma(identifier_identifier_map);
      process_identifier id1= get_last(id1_,identifier_identifier_map);
      process_identifier id2= get_last(id2_,identifier_identifier_map);
      if (id1==initial_process)
      {
        id1.swap(id2);
      }
      if (id1!=id2)
      {
        identifier_identifier_map[id1]=id2;
      }
    }

    static
    void complete_proc_identifier_map(std::map< process_identifier, process_identifier >& identifier_identifier_map)
    {
      std::map< process_identifier, process_identifier > new_identifier_identifier_map;
      for(const std::pair<const process_identifier, process_identifier >& p: identifier_identifier_map)
      {
        new_identifier_identifier_map[p.first]=get_last(p.second,identifier_identifier_map);
      }
      identifier_identifier_map.swap(new_identifier_identifier_map);
#ifndef NDEBUG
      /* In the result no right hand side occurs as the left hand side of identifier_identifier_map */
      using identifier_identifier_pair = std::pair<const process_identifier, process_identifier>;
      for(const identifier_identifier_pair& p: identifier_identifier_map)
      {
        assert(identifier_identifier_map.count(p.second)==0);
      }
#endif
    }

    struct make_substitution
    {
      using result_type = process_identifier;
      using argument_type = process_identifier;
      const std::map< process_identifier, process_identifier >& m_map;

      make_substitution(const std::map< process_identifier, process_identifier >& map)
        : m_map(map)
      {}

      process_identifier operator()(const process_identifier& id) const
      {
        const std::map< process_identifier, process_identifier >::const_iterator i=m_map.find(id);
        if (i==m_map.end()) // Not found
        {
          return id;
        }
        return i->second;
      }
    };

    std::set< process_identifier >
          minimize_set_of_reachable_process_identifiers(const std::set<process_identifier>& reachable_process_identifiers,
                                                        const process_identifier& initial_process)
    {
      assert(reachable_process_identifiers.count(initial_process)>0);
      using parameters_process_pair = std::pair<variable_list, process_expression>;
      using mapping_type = std::map<std::pair<variable_list, process_expression>, process_identifier>;
      using mapping_type_pair = std::pair<const std::pair<variable_list, process_expression>, process_identifier>;

      /* First put the identifiers in reachable_process_identifiers in the process mapping */
      mapping_type process_mapping;
      std::map< process_identifier, process_identifier > identifier_identifier_map;
      for(const process_identifier& id: reachable_process_identifiers)
      {
        objectdatatype& object=objectIndex(id);
        const parameters_process_pair p(object.parameters,object.processbody);
        mapping_type::const_iterator i=process_mapping.find(p);
        if (i==process_mapping.end())   // Not found.
        {
          process_mapping[p]=id;
        }
        else
        {
          set_proc_identifier_map(identifier_identifier_map,id,i->second,initial_process);
        }
      }
      complete_proc_identifier_map(identifier_identifier_map);

      while (!identifier_identifier_map.empty())
      {
        /* Move the elements of process_mapping into new_process_mapping, under application
         * of the identifier_identifier_map. If this yields new process_identifiers to be identified store
         * these in new_proc_identifier_map. Do this until no process_identifiers are mapped onto
         * each other anymore. */
        std::map< process_identifier, process_identifier > new_identifier_identifier_map;
        mapping_type new_process_mapping;
        for(const mapping_type_pair& p: process_mapping)
        {
          const make_substitution sigma(identifier_identifier_map);
          const parameters_process_pair p_new(p.first.first, replace_process_identifiers( p.first.second, sigma));
          mapping_type::const_iterator i=new_process_mapping.find(p_new);
          if (i==new_process_mapping.end())   // Not found.
          {
            new_process_mapping[p_new]=sigma(p.second);
          }
          else
          {
            set_proc_identifier_map(new_identifier_identifier_map,sigma(p.second),sigma(i->second),initial_process);
          }
        }
        complete_proc_identifier_map(new_identifier_identifier_map);
        identifier_identifier_map.swap(new_identifier_identifier_map);
        process_mapping.swap(new_process_mapping);
      }

      /* Store the pairs from process mapping into result */
      /* std::map< process_identifier, process_expression > result;
      for(const mapping_type_pair& p: process_mapping)
      {
        result[p.second]=p.first.second;
      }
      assert(result.count(initial_process)>0);
      return result; */
      std::set< process_identifier > result;
      for(const mapping_type_pair& p: process_mapping)
      {
        result.insert(p.second);
        objectdatatype& object=objectIndex(p.second);
        object.processbody=p.first.second;
      }
      assert(result.count(initial_process)>0);
      return result;
    }

    /****************   remove_stochastic_operators_from_front **********************/

   class process_pid_pair
   {
     protected:
       process_expression m_process_body;
       process_identifier m_pid;

     public:
       process_pid_pair(const process_expression& process_body, const process_identifier& pid)
         : m_process_body(process_body), m_pid(pid)
       {}

       process_pid_pair(const process_pid_pair& other) = default;
       process_pid_pair(process_pid_pair&& other) = default;
       process_pid_pair& operator=(const process_pid_pair& other) = default;
       process_pid_pair& operator=(process_pid_pair&& other) = default;

       const process_expression& process_body() const
       {
         return m_process_body;
       }

       const process_identifier& process_id() const
       {
         return m_pid;
       }
    };  // end process_pid_pair

    /* Take the process expressions from process_equations and transform these such that
       no process equation has an initial distribution. The changed equations are reported
       back, including adapted parameter lists because sometimes the variables bound in the
       stochastic operators need to be added.
       The initial distribution of the initial process is reported back in
       initial_stochastic_distribution, including if necessary adapted sets of arguments. */

    std::set< process_identifier >
              remove_stochastic_operators_from_front(
                   const std::set< process_identifier >& reachable_process_identifiers,
                   process_identifier& initial_process_id, // If parameters change, another initial_process_id will be substituted.
                   stochastic_distribution& initial_stochastic_distribution)
    {
      /* First obtain the stochastic distribution for each process variable. */
      std::map< process_identifier, process_pid_pair > processes_with_stochastic_distribution_first;
      std::set< process_identifier > result;
      for(const process_identifier& p: reachable_process_identifiers)
      {
        objectdatatype& object = objectIndex(p);
        process_expression proc_=obtain_initial_distribution_term(object.processbody);
        if (!is_stochastic_operator(proc_))
        {
          processes_with_stochastic_distribution_first.insert(std::pair< process_identifier, process_pid_pair >(p, process_pid_pair(proc_,p)));
          result.insert(p);
        }
        else
        {
          const stochastic_operator& proc=down_cast<const stochastic_operator>(proc_);
          assert(!is_process_instance_assignment(proc.operand()));
          objectdatatype& object=objectIndex(p);
          maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
          variable_list vars=proc.variables();
          alphaconvert(vars,local_sigma, vars + object.parameters, data_expression_list());

          const process_identifier newproc=newprocess(
                                           vars + object.parameters,
                                           process::replace_variables_capture_avoiding_with_an_identifier_generator(proc.operand(),
                                                                                       local_sigma,
                                                                                       fresh_identifier_generator),
                                           pCRL,
                                           canterminatebody(proc.operand()),
                                           containstimebody(proc.operand()));
          // calculate the substitution to be applied on proc.distribution() which is moved outside the
          // body the process.
          processes_with_stochastic_distribution_first.insert(
                        std::pair< process_identifier, process_pid_pair >
                                 (p,
                                  process_pid_pair(stochastic_operator(
                                                        vars,
                                                        replace_variables_capture_avoiding_alt(proc.distribution(), local_sigma),
                                                        proc.operand()),
                                                   newproc)));
          result.insert(newproc);
        }
      }

      for(const process_identifier& p: reachable_process_identifiers)
      {
        objectdatatype& object=objectIndex(processes_with_stochastic_distribution_first.at(p).process_id());
        assert(!is_stochastic_operator(object.processbody));
        object.processbody=transform_initial_distribution_term(object.processbody,processes_with_stochastic_distribution_first);
        assert(!is_stochastic_operator(object.processbody));
      }

      // Adapt the initial process
      process_expression initial_distribution=processes_with_stochastic_distribution_first.at(initial_process_id).process_body();
      if (is_stochastic_operator(initial_distribution))
      {
        const stochastic_operator sto=atermpp::down_cast<stochastic_operator>(initial_distribution);
        initial_stochastic_distribution = stochastic_distribution(sto.variables(), sto.distribution());
        initial_process_id=processes_with_stochastic_distribution_first.at(initial_process_id).process_id();
      }
      else
      {
        initial_stochastic_distribution = stochastic_distribution(variable_list(), sort_real::real_one());
      }

      return result;
    }


    /**************** Collectparameterlist ******************************/

    bool alreadypresent(variable& var,
                        const variable_list& vl, 
                        mutable_indexed_substitution<>& parameter_renaming)
    {
      /* Note: variables can be different, although they have the
         same string, due to different types. If they have the
         same string, but different types, the conflict must
         be resolved by renaming the name of the variable.
         Note that var is changed to the new variable name. */

      for(const variable& v: vl)
      {
        if (v.name()==var.name())
        {
          if (v.sort()==var.sort())
          {
            return true; // The variable is present. 
          }
          else 
          {
            const variable new_var(fresh_identifier_generator(var.name()),var.sort());
            parameter_renaming[var]=new_var;
            var=new_var;
            return false; // variable var is renamed, as name clashed with a variable with a different sort. Not present.
          }
        }
      }
      return false; // The variable var is not present. 
    }

    variable_list joinparameters(const variable_list& par1,
                                 const variable_list& par2,
                                 mutable_indexed_substitution<>& parameter_renaming)
    {
      if (par2.empty())
      {
        return par1;
      }

      variable var2=par2.front();
      assert(is_variable(var2));

      variable_list result=joinparameters(par1,par2.tail(), parameter_renaming);
      if (alreadypresent(var2,par1, parameter_renaming))
      {
        return result;
      }

      result.push_front(var2);
      return result;
    }

    variable_list collectparameterlist(const std::set< process_identifier >& pCRLprocs, 
                                       mutable_indexed_substitution<>& parameter_renaming)
    {
      variable_list parameters;
      for (const process_identifier& p: pCRLprocs)
      {
        const objectdatatype& object=objectIndex(p);
        parameters=joinparameters(parameters,object.parameters, parameter_renaming);
      }
      return parameters;
    }

    /****************  Declare local datatypes  ******************************/

    void declare_control_state(
      const std::set < process_identifier >& pCRLprocs)
    {
      create_enumeratedtype(pCRLprocs.size());
    }

    class stackoperations
    {
      public:
        variable_list parameter_list;
        sort_expression stacksort;
        sort_expression_list sorts;
        function_symbol_list get;
        data::function_symbol push;
        data::function_symbol emptystack;
        data::function_symbol empty;
        data::function_symbol pop;
        data::function_symbol getstate;
        stackoperations* next;

        // Stack operations are not supposed to be copied.
        stackoperations(const stackoperations& )=delete;
        stackoperations& operator=(const stackoperations& )=delete;

        stackoperations(const variable_list& pl,
                        specification_basic_type& spec)
        {
          parameter_list=pl;
          next=spec.stack_operations_list;
          spec.stack_operations_list=this;

          //create structured sort
          //  Stack = struct emptystack?is_empty
          //               | push(getstate: Pos, getx1: S1, ..., getxn: Sn, pop: Stack)
          //               ;

          basic_sort stack_sort_alias(spec.fresh_identifier_generator("Stack"));
          structured_sort_constructor_argument_vector sp_push_arguments;
          for (const variable& v: pl)
          {
            sp_push_arguments.emplace_back(spec.fresh_identifier_generator("get" + std::string(v.name())), v.sort());
            sorts.push_front(v.sort());
          }
          sp_push_arguments.emplace_back(spec.fresh_identifier_generator("pop"), stack_sort_alias);
          sorts=reverse(sorts);
          structured_sort_constructor sc_push(spec.fresh_identifier_generator("push"), sp_push_arguments);
          structured_sort_constructor sc_emptystack(spec.fresh_identifier_generator("emptystack"),spec.fresh_identifier_generator("isempty"));

          structured_sort_constructor_vector constructors(1,sc_push);
          constructors.push_back(sc_emptystack);
          //add data declarations for structured sort
          spec.data.add_alias(alias(stack_sort_alias,structured_sort(constructors)));
          stacksort=data::normalize_sorts(stack_sort_alias,spec.data);
          push=sc_push.constructor_function(stack_sort_alias);
          emptystack=sc_emptystack.constructor_function(stack_sort_alias);
          empty=sc_emptystack.recogniser_function(stack_sort_alias);
          const std::vector< data::function_symbol > projection_functions =
            sc_push.projection_functions(stack_sort_alias);
          pop=projection_functions.back();
          getstate=projection_functions.front();
          get=function_symbol_list(projection_functions.begin()+1,projection_functions.end()-1);
        }

        ~stackoperations() = default;
    };

    class stacklisttype
    {
      public:
        stackoperations* opns;
        variable_list parameters;
        variable stackvar;
        std::size_t no_of_states;
        /* the boolean state variables occur in reverse
           order, i.e. the least significant first, whereas
           in parameter lists, the order is reversed. */
        variable_list booleanStateVariables;


        /* All datatypes for different stacks that are being generated
           are stored in the following list, such that it can be investigated
           whether a suitable stacktype already exist, before generating a new
           one */

        stacklisttype(const stacklisttype& )=delete;
        stacklisttype& operator=(const stacklisttype& )=delete;

        static
        stackoperations* find_suitable_stack_operations(
          const variable_list& parameters,
          stackoperations* stack_operations_list)
        {
          if (stack_operations_list==nullptr)
          {
            return nullptr;
          }
          if (parameters==stack_operations_list->parameter_list)
          {
            return stack_operations_list;
          }
          return find_suitable_stack_operations(parameters,stack_operations_list->next);
        }

        /// \brief Constructor
        stacklisttype(const variable_list& parlist,
                      specification_basic_type& spec,
                      const bool regular,
                      const std::set < process_identifier >& pCRLprocs,
                      const bool singlecontrolstate)
        {
          assert(pCRLprocs.size()>0);
          parameters=parlist;

          no_of_states=pCRLprocs.size();
          process_identifier last= *pCRLprocs.begin();
          const std::string s3((spec.options.statenames)?std::string(last.name()):std::string("s"));
          if ((spec.options.binary) && (spec.options.newstate))
          {
            std::size_t i=spec.upperpowerof2(no_of_states);
            for (; i>0 ; i--)
            {
              variable name(spec.fresh_identifier_generator("bst"),sort_bool::bool_());
              spec.insertvariable(name,true);
              booleanStateVariables.push_front(name);
            }
          }

          if (regular)
          {
            opns=nullptr;
            if (spec.options.newstate)
            {
              if (!spec.options.binary)
              {
                if (!singlecontrolstate)
                {
                  const std::size_t e=spec.create_enumeratedtype(no_of_states);
                  stackvar=variable(spec.fresh_identifier_generator(s3), spec.enumeratedtypes[e].sortId);
                }
                else
                {
                  /* Generate a stackvariable that is never used */
                  stackvar=variable(spec.fresh_identifier_generator("Never_used"), sort_bool::bool_());
                }
              }
              else
              {
                stackvar=variable(spec.fresh_identifier_generator(s3),sort_bool::bool_());
              }
            }
            else
            {
              stackvar=variable(spec.fresh_identifier_generator(s3), sort_pos::pos());
            }
            spec.insertvariable(stackvar,true);
          }
          else
          {
            if (spec.options.newstate)
            {
              throw mcrl2::runtime_error("cannot combine stacks with " +
                                         (spec.options.binary?std::string("binary"):std::string("an enumerated type")));
            }
            opns=find_suitable_stack_operations(parlist,spec.stack_operations_list);

            if (opns!=nullptr)
            {
              stackvar=variable(spec.fresh_identifier_generator(s3),opns->stacksort);
              spec.insertvariable(stackvar,true);
            }
            else
            {
              variable_list temp=parlist;
              temp.push_front(variable("state",sort_pos::pos()));
              opns=(stackoperations*) new stackoperations(temp,spec);
              stackvar = variable(spec.fresh_identifier_generator(s3), opns->stacksort);
              spec.insertvariable(stackvar,true);
            }
          }
        }

        ~stacklisttype() = default;
    };

    bool is_global_variable(const data_expression& d) const
    {
      return is_variable(d) && global_variables.count(atermpp::down_cast<variable>(d)) > 0;
    }

    data_expression getvar(const variable& var,
                           const stacklisttype& stack) const
    {
      /* first search whether the variable is a free process variable */
      if (global_variables.count(var)>0)
      {
        return var;
      }

      /* otherwise find out whether the variable matches a parameter */
      function_symbol_list::const_iterator getmappings=stack.opns->get.begin();
      for (variable_list::const_iterator walker=stack.parameters.begin() ;
           walker!=stack.parameters.end() ; ++walker,++getmappings)
      {
        if (*walker==var)
        {
          return application(*getmappings,stack.stackvar);
        }
        assert(getmappings!=stack.opns->get.end());
      }
      assert(0); /* We cannot end up here, because that means that we
                    are looking for in non-existing variable */
      return var;
    }

    assignment_list processencoding(
      std::size_t i,
      const assignment_list& t1,
      const stacklisttype& stack)
    {
      assert(i>0);
      assignment_list t(t1);
      if (!options.newstate)
      {
        assignment_list result=t;
        result.push_front(assignment(stack.stackvar,sort_pos::pos(i)));
        return result;
      }

      i=i-1; /* below we count from 0 instead from 1 as done in the
                first version of the prover */

      if (!options.binary)
      {
        const std::size_t e=create_enumeratedtype(stack.no_of_states);
        data_expression_list l=enumeratedtypes[e].elementnames;
        for (; i>0 ; i--)
        {
          assert(l.size()>0);
          l.pop_front();
        }
        assignment_list result=t;
        assert(l.size()>0);
        result.push_front(assignment(stack.stackvar,l.front()));
        return result;
      }
      /* else a sequence of boolean values needs to be generated,
         representing the value i, when there are l->n elements */
      {
        std::size_t k=upperpowerof2(stack.no_of_states);
        variable_list::const_iterator boolean_state_variables=stack.booleanStateVariables.begin();
        for (; k>0 ; k--, ++boolean_state_variables)
        {
          if ((i % 2)==0)
          {
            t.push_front(assignment(*boolean_state_variables,sort_bool::false_()));
            i=i/2;
          }
          else
          {
            t.push_front(assignment(*boolean_state_variables,sort_bool::true_()));
            i=(i-1)/2;
          }
        }
        return t;
      }
    }

    data_expression_list processencoding(
      std::size_t i,
      const data_expression_list& t1,
      const stacklisttype& stack)
    {
      data_expression_list t(t1);
      if (!options.newstate)
      {
        data_expression_list result=t;
        result.push_front(sort_pos::pos(i));
        return result;
      }

      i=i-1; /* below we count from 0 instead from 1 as done in the
                first version of the linearizer */

      if (!options.binary)
      {
        const std::size_t e=create_enumeratedtype(stack.no_of_states);
        data_expression_list l(enumeratedtypes[e].elementnames);
        for (; i>0 ; i--)
        {
          l.pop_front();
        }
        data_expression_list result=t;
        result.push_front(l.front());
        return result;
      }
      /* else a sequence of boolean values needs to be generated,
         representing the value i, when there are l->n elements */
      {
        std::size_t k=upperpowerof2(stack.no_of_states);
        variable_list::const_iterator boolean_state_variables=stack.booleanStateVariables.begin();
        for (; k>0 ; k--, ++boolean_state_variables)
        {
          if ((i % 2)==0)
          {
            t.push_front(sort_bool::false_());
            i=i/2;
          }
          else
          {
            t.push_front(sort_bool::true_());
            i=(i-1)/2;
          }
        }
        return t;
      }
    }

    data_expression correctstatecond(
      const process_identifier& procId,
      const std::set < process_identifier >& pCRLproc,
      const stacklisttype& stack,
      int regular)
    {
      std::size_t i=1;
      for (const process_identifier& p: pCRLproc)
      {
        if (p==procId)
        {
          break;
        }
        ++i;
      }
      /* i is the index of the current process */

      if (!options.newstate)
      {
        if (regular)
        {
          return equal_to(stack.stackvar, processencoding(i,assignment_list(),stack).front().rhs());
        }
        return equal_to(
                 application(stack.opns->getstate,stack.stackvar),
                 processencoding(i,assignment_list(),stack).front().rhs());
      }

      if (!options.binary) /* Here a state encoding using enumerated types
                              must be declared */
      {
        create_enumeratedtype(stack.no_of_states);
        if (regular)
        {
          return equal_to(stack.stackvar,
                          processencoding(i,assignment_list(),stack).front().rhs());
        }
        return equal_to(
                 application(stack.opns->getstate, stack.stackvar),
                 processencoding(i,assignment_list(),stack).front().rhs());
      }

      /* in this case we must encode the condition using
         boolean variables */

      const variable_list vars=stack.booleanStateVariables;

      i=i-1; /* start counting from 0, instead from 1 */
      data_expression t3(sort_bool::true_());
      for (const variable& v: vars)
      {
        if ((i % 2)==0)
        {
          t3=lazy::and_(lazy::not_(v),t3);
          i=i/2;
        }
        else
        {
          t3=lazy::and_(v,t3);
          i=(i-1)/2;
        }

      }
      assert(i==0);
      return t3;
    }

    data_expression adapt_term_to_stack(
      const data_expression& t,
      const stacklisttype& stack,
      const variable_list& vars,
      const variable_list& stochastic_variables)
    {
      if (is_function_symbol(t))
      {
        return t;
      }

      if (is_machine_number(t))
      {
        return t;
      }

      if (is_variable(t))
      {
        if (std::find(vars.begin(),vars.end(),t)!=vars.end())
        {
          /* t occurs in vars, so, t does not have to be reconstructed
             from the stack */
          return t;
        }
        else
        if (std::find(stochastic_variables.begin(),stochastic_variables.end(),t)!=stochastic_variables.end())
        {
          /* t occurs in stochastic_variables, so, t does not have to be reconstructed
             from the stack */
          return t;
        }
        else
        {
          return getvar(atermpp::down_cast<variable>(t), stack);
        }
      }

      if (is_application(t))
      {
        const application&a=atermpp::down_cast<application>(t);
        return application(
                 adapt_term_to_stack(a.head(),stack,vars,stochastic_variables),
                 adapt_termlist_to_stack(a.begin(),a.end(),stack,vars,stochastic_variables));
      }

      if (is_abstraction(t))
      {
        const abstraction& abs_t=down_cast<abstraction>(t);
        return abstraction(
                 abs_t.binding_operator(),
                 abs_t.variables(),
                 adapt_term_to_stack(abs_t.body(),stack,abs_t.variables() + vars,stochastic_variables));
      }

      if (is_where_clause(t))
      {
        const where_clause where_t(t);
        const assignment_list old_assignments=reverse(where_t.assignments());
        variable_list new_vars=vars;
        assignment_list new_assignments;
        for (const assignment& a: old_assignments)
        {
          new_vars.push_front(a.lhs());
          new_assignments.push_front(
                             assignment(
                               a.lhs(),
                               adapt_term_to_stack(a.rhs(),stack,vars,stochastic_variables)));

        }
        return where_clause(
                 adapt_term_to_stack(where_t,stack,new_vars,stochastic_variables),
                 new_assignments);

      }

      assert(0);  // expected a term;
      return t;   // in case of non-debug mode, try to return something as decent as possible.
    }

    template <typename Iterator>
    data_expression_vector adapt_termlist_to_stack(
      Iterator begin,
      const Iterator& end,
      const stacklisttype& stack,
      const variable_list& vars,
      const variable_list& stochastic_variables)
    {
      data_expression_vector result;
      for (; begin != end; ++begin)
      {
        result.push_back(adapt_term_to_stack(*begin,stack, vars,stochastic_variables));
      }
      return result;
    }

    /// Adapt multiactions to stack by traversing the list.
    /// Note: the result is sorted w.r.t action_compare().
    action_list
    adapt_multiaction_to_stack(const action_list& multiAction, const stacklisttype& stack, const variable_list& vars)
    {
      action_vector result;

      for (const action& act : multiAction)
      {
        const data_expression_vector vec(
            adapt_termlist_to_stack(act.arguments().begin(), act.arguments().end(), stack, vars, variable_list()));
        result.emplace_back(act.label(), data_expression_list(vec.begin(), vec.end()));
      }

      // As the arguments change, sort order w.r.t. action_label is preserved, but order w.r.t. arguments is changed.
      // resort the list.
      std::sort(result.begin(), result.end(), action_compare());

      return action_list(result.begin(), result.end());
    }

    data_expression representative_generator_internal(const sort_expression& s, const bool allow_dont_care_var=true)
    {
      if ((!options.noglobalvars) && allow_dont_care_var)
      {
        const variable newVariable(fresh_identifier_generator("dc"),s);
        insertvariable(newVariable,true);
        global_variables.insert(newVariable);
        return newVariable;
      }
      return representative_generator(data)(s);
    }

    data_expression find_(
      const variable& s,
      const assignment_list& args,
      const stacklisttype& stack,
      const variable_list& vars,
      const std::set<variable>& free_variables_in_body,
      const variable_list& stochastic_variables)
    {
      /* We generate the value for variable s in the list of
         the parameters of the process. If s is equal to some
         variable in pars, it is an argument of the current
         process, and it must be replaced by the corresponding
         argument in args.
           If s does not occur in pars, it must be replaced
         by a dummy value.
      */
      for (const assignment& a: args)
      {
        if (s==a.lhs())
        {
          return adapt_term_to_stack(a.rhs(),stack,vars,stochastic_variables);
        }
      }

      if (free_variables_in_body.find(s)==free_variables_in_body.end())
      {
        const data_expression result=representative_generator_internal(s.sort());
        return adapt_term_to_stack(result,stack,vars,stochastic_variables);
      }
      return adapt_term_to_stack(s,stack,vars,stochastic_variables);
    }


    data_expression_list findarguments(
      const variable_list& pars,
      const variable_list& parlist,
      const assignment_list& args,
      const data_expression_list& t2,
      const stacklisttype& stack,
      const variable_list& vars,
      const std::set<variable>& free_variables_in_body,
      const variable_list& stochastic_variables)
    {
      if (parlist.empty())
      {
        return t2;
      }
      data_expression_list result=findarguments(pars,parlist.tail(),args,t2,stack,vars,free_variables_in_body,stochastic_variables);
      data_expression rhs=find_(parlist.front(),args,stack,vars,free_variables_in_body,stochastic_variables);

      result.push_front(rhs);
      return result;
    }

    assignment_list find_dummy_arguments(
      const variable_list& parlist,   // The list of all parameters.
      const assignment_list& args,
      const std::set<variable>& free_variables_in_body,
      const variable_list& stochastic_variables)
    {
      std::map<variable,data_expression> assignment_map;
      for(const assignment& a: args)
      {
        assignment_map[a.lhs()]=a.rhs();
      }

      assignment_vector result;
      for(const variable& v: parlist)
      {
        if (std::find(stochastic_variables.begin(),stochastic_variables.end(),v)!=stochastic_variables.end())
        {
          // v is a stochastic variable. Insert the identity assignment.
          result.emplace_back(v, v);
        }
        else if (free_variables_in_body.find(v)==free_variables_in_body.end())
        {
          {
            // The variable *i must get a default value.
            const data_expression rhs=representative_generator_internal(v.sort());
            result.emplace_back(v, rhs);
          }
        }
        else
        {
          const std::map<variable,data_expression>::iterator k=assignment_map.find(v);
          if (k!=assignment_map.end())  // There is assignment for v. Use it.
          {
            result.emplace_back(k->first, k->second);
            assignment_map.erase(k);
          }
        }

      }
      return assignment_list(result.begin(), result.end());
    }



    assignment_list push_regular(
      const process_identifier& procId,
      const assignment_list& args,
      const stacklisttype& stack,
      const std::set < process_identifier >& pCRLprocs,
      bool singlestate,
      const variable_list& stochastic_variables)
    {
      objectdatatype& object=objectIndex(procId);
      const assignment_list t=find_dummy_arguments(stack.parameters,args,get_free_variables(object),stochastic_variables);

      if (singlestate)
      {
        return t;
      }

      std::size_t i=1;
      for (const process_identifier& p: pCRLprocs)
      {
        if (p==procId)
        {
          break;
        }
        ++i;
      }
      return processencoding(i,t,stack);
    }


    assignment_list make_procargs_regular(
      const process_expression& t,
      const stacklisttype& stack,
      const std::set < process_identifier >& pcrlprcs,
      const bool singlestate,
      const variable_list& stochastic_variables)
    {
      /* t is a sequential composition of process variables */

      if (is_seq(t))
      {
        throw mcrl2::runtime_error("Process is not regular, as it has stacking vars: " + process::pp(t) + ".");
      }

      if (is_process_instance_assignment(t))
      {
        const process_identifier& procId=atermpp::down_cast<process_instance_assignment>(t).identifier();
        const assignment_list& t1=atermpp::down_cast<process_instance_assignment>(t).assignments();
        return push_regular(procId,
                            t1,
                            stack,
                            pcrlprcs,
                            singlestate,
                            stochastic_variables);
      }

      throw mcrl2::runtime_error("Expect seq or name constructing a recursive invocation: " + process::pp(t) + ".");
    }

    data_expression push_stack(
      const process_identifier& procId,
      const assignment_list& args,
      const data_expression_list& t2,
      const stacklisttype& stack,
      const std::set < process_identifier >& pCRLprocs,
      const variable_list& vars,
      const variable_list& stochastic_variables)
    {
      objectdatatype& object=objectIndex(procId);
      const data_expression_list t=findarguments(object.parameters,
                                                 stack.parameters,
                                                 args,t2,stack,vars,
                                                 get_free_variables(object),
                                                 stochastic_variables);

      std::size_t i=1;
      for (const process_identifier& p: pCRLprocs)
      {
        if (p==procId)
        {
          break;
        }
        ++i;
      }
      const data_expression_list l=processencoding(i,t,stack);
      assert(l.size()==function_sort(stack.opns->push.sort()).domain().size());
      return application(stack.opns->push,l);
    }

    data_expression make_procargs_stack(
      const process_expression& t,
      const stacklisttype& stack,
      const std::set < process_identifier >& pcrlprcs,
      const variable_list& vars,
      const variable_list& stochastic_variables)
    {
      /* t is a sequential composition of process variables */

      if (is_seq(t))
      {
        const process_instance_assignment process=atermpp::down_cast<process_instance_assignment>(seq(t).left());
        const process_expression& process2=down_cast<seq>(t).right();
        const process_identifier& procId=process.identifier();
        const assignment_list& t1=process.assignments();

        if (objectIndex(procId).canterminate)
        {
          const data_expression stackframe=make_procargs_stack(process2,stack,pcrlprcs, vars,stochastic_variables);
          return push_stack(procId,t1, data_expression_list({ stackframe }),stack,pcrlprcs,vars,stochastic_variables);
        }

        return push_stack(procId,t1, data_expression_list({ data_expression(stack.opns->emptystack) }),
                                           stack,pcrlprcs,vars,stochastic_variables);
      }

      if (is_process_instance_assignment(t))
      {
        const process_identifier procId=process_instance_assignment(t).identifier();
        const assignment_list t1=process_instance_assignment(t).assignments();

        if (objectIndex(procId).canterminate)
        {
          return push_stack(procId,
                            t1,
                            data_expression_list({ data_expression(application(stack.opns->pop,stack.stackvar)) }),
                            stack,
                            pcrlprcs,
                            vars,
                            stochastic_variables);
        }
        return push_stack(procId,
                          t1,
                          data_expression_list({ data_expression(stack.opns->emptystack) }),
                          stack,
                          pcrlprcs,
                          vars,
                          stochastic_variables);
      }

      throw mcrl2::runtime_error("Expect seq or name putting processes on a stack: " + process::pp(t) +".");
    }

    assignment_list make_procargs(
      const process_expression& t,
      const stacklisttype& stack,
      const std::set < process_identifier >& pcrlprcs,
      const variable_list& vars,
      const bool regular,
      const bool singlestate,
      const variable_list& stochastic_variables)
    {
      if (regular)
      {
        return make_procargs_regular(t,stack,pcrlprcs,singlestate,stochastic_variables);
      }
      /* return a stackframe */
      data_expression sf=make_procargs_stack(t,stack,pcrlprcs,vars,stochastic_variables);
      return assignment_list({ assignment(stack.stackvar,sf) });
    }

    static
    bool occursin(const variable& name,
                  const variable_list& pars)
    {
      assert(is_variable(name));
      for (const variable& v: pars)
      {
        if (name.name()==v.name())
        {
          return true;
        }
      }
      return false;
    }


    data_expression_list pushdummy_regular_data_expressions(
      const variable_list& pars,
      const stacklisttype& stack)
    {
      /* stack.parameters is the total list of parameters of the
         aggregated pCRL process. The variable pars contains
         the list of all variables occuring in the initial
         process. */

      data_expression_vector result;
      for(const variable& par: stack.parameters)
      {
        /* // Check whether it is a stochastic variable.
        if (std::find(stochastic_variables.begin(),stochastic_variables.end(),par)!=pars.end())
        {
          result.push_back(assignment(par,par));
        }
        // Otherwise, check that is is an ordinary parameter.
        else */
        if (std::find(pars.begin(),pars.end(),par)!=pars.end())
        {
          result.push_back(par);
        }
        /* otherwise the value of this argument is irrelevant, so
           make it a don't care variable. */
        else
        {
          result.push_back(representative_generator_internal(par.sort()));
        }
      }
      return data_expression_list(result.begin(), result.end());
    }

    assignment_list pushdummy_regular(
      const variable_list& pars,
      const stacklisttype& stack,
      const variable_list& stochastic_variables)
    {
      /* stack.parameters is the total list of parameters of the
         aggregated pCRL process. The variable pars contains
         the list of all variables occuring in the initial
         process. */

      assignment_vector result;
      for(const variable& par: stack.parameters)
      {
        // Check whether it is a stochastic variable.
        if (std::find(stochastic_variables.begin(),stochastic_variables.end(),par)!=pars.end())
        {
          result.emplace_back(par, par);
        }
        // Otherwise, check that is is an ordinary parameter.
        else if (std::find(pars.begin(),pars.end(),par)!=pars.end())
        {
        }
        /* otherwise the value of this argument is irrelevant, so
           make it a don't care variable. */
        else
        {
          result.emplace_back(par, representative_generator_internal(par.sort()));
        }
      }
      return assignment_list(result.begin(), result.end());
    }

    data_expression_list pushdummyrec_stack(
      const variable_list& totalpars,
      const variable_list& pars,
      const stacklisttype& stack,
      const variable_list& stochastic_variables)
    {
      /* totalpars is the total list of parameters of the
         aggregated pCRL process. The variable pars contains
         the list of all variables occuring in the initial
         process. */

      if (totalpars.empty())
      {
        return data_expression_list({ data_expression(stack.opns->emptystack) });
      }

      const variable& par=totalpars.front();
      if (std::find(pars.begin(),pars.end(),par)!=pars.end())
      {
        data_expression_list result=pushdummyrec_stack(totalpars.tail(),pars,stack,stochastic_variables);
        result.push_front(par);
        return result;
      }
      /* Check whether the parameter par refers to a stochastic variables */
      if (std::find(stochastic_variables.begin(),stochastic_variables.end(),par)!=pars.end())
      {
        data_expression_list result=pushdummyrec_stack(totalpars.tail(),pars,stack,stochastic_variables);
        result.push_front(par);
        return result;
      }
      /* otherwise the value of this argument is irrelevant, so
         make it Nil, if a regular translation is made. If a translation
         with stacks is made, then yield a default `unique' term. */
      data_expression_list result=pushdummyrec_stack(totalpars.tail(),pars,stack,stochastic_variables);
      result.push_front(representative_generator_internal(par.sort()));
      return result;
    }

    data_expression_list pushdummy_stack(
      const variable_list& parameters,
      const stacklisttype& stack,
      const variable_list& stochastic_variables)
    {
      return pushdummyrec_stack(stack.parameters,parameters,stack,stochastic_variables);
    }

    data_expression_list make_initialstate(
      const process_identifier& initialProcId,
      const stacklisttype& stack,
      const std::set < process_identifier >& pcrlprcs,
      const bool regular,
      const bool singlecontrolstate,
      const stochastic_distribution& initial_stochastic_distribution)
    {
      std::size_t i=1;
      for (const process_identifier& p: pcrlprcs)
      {
        if (p==initialProcId)
        {
          break;
        }
        ++i;
      }
      /* i is the index of the initial state */

      if (regular)
      {
        data_expression_list result=
          pushdummy_regular_data_expressions(objectIndex(initialProcId).parameters,
                                             stack);
        if (!singlecontrolstate)
        {
          return processencoding(i,result,stack);
        }
        return result;
      }
      else
      {
        data_expression_list result=
                pushdummy_stack(objectIndex(initialProcId).parameters,
                                stack,
                                initial_stochastic_distribution.variables());
        const data_expression_list l=processencoding(i,result,stack);
        assert(l.size()==function_sort(stack.opns->push.sort()).domain().size());
        return data_expression_list({ application(stack.opns->push,l) });
      }
    }

    /*************************  Routines for summands  **************************/

    assignment_list dummyparameterlist(const stacklisttype& stack,
                                       const bool singlestate)
    {
      if (singlestate)
      {
        return assignment_list();
      }

      return processencoding(1,assignment_list(),stack); /* Take 1 as dummy indicator */
    }


    void insert_summand(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const variable_list& sumvars,
      const data_expression& condition,
      const action_list& multiAction,
      const data_expression& actTime,
      const stochastic_distribution& distribution,
      const assignment_list& procargs,
      const bool has_time,
      const bool is_deadlock_summand)
    {
      assert(distribution != stochastic_distribution());
      assert(std::is_sorted(multiAction.begin(), multiAction.end(), action_compare()));

      const data_expression rewritten_condition=RewriteTerm(condition);
      if (rewritten_condition==sort_bool::false_())
      {
        return;
      }

      if (is_deadlock_summand)
      {
        insert_timed_delta_summand(action_summands,
                                   deadlock_summands,
                                   deadlock_summand(sumvars,
                                                    rewritten_condition,
                                                    has_time?deadlock(actTime):deadlock()),options.ignore_time);
      }
      else
      {
        action_summands.emplace_back(sumvars,
            rewritten_condition,
            has_time ? multi_action(multiAction, actTime) : multi_action(multiAction),
            procargs,
            stochastic_distribution(distribution.variables(), RewriteTerm(distribution.distribution())));
      }
    }


    void add_summands(
      const process_identifier& procId,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      process_expression summandterm,
      const std::set < process_identifier >& pCRLprocs,
      const stacklisttype& stack,
      const bool regular,
      const bool singlestate,
      const variable_list& process_parameters)
      
    {
      data_expression atTime;
      action_list multiAction;
      bool is_delta_summand=false;
      bool has_time=false;

      if (isDeltaAtZero(summandterm))
      {
        // delta@0 does not need to be added.
        return;
      }

      /* remove the sum operators; collect the sum variables in the
         list sumvars */

      variable_list sumvars;
      while (is_sum(summandterm))
      {
        sumvars=sum(summandterm).variables() + sumvars;
        summandterm=sum(summandterm).operand();
      }

      // Check whether the variables in sumvars clash with the parameter list,
      // and rename the summandterm accordingly.
      maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
      alphaconvert(sumvars,local_sigma,process_parameters,data_expression_list());
      summandterm=substitute_pCRLproc(summandterm, local_sigma);


      /* translate the condition */

      data_expression condition1;
      if (regular && singlestate)
      {
        condition1=sort_bool::true_();
      }
      else
      {
        condition1=correctstatecond(procId,pCRLprocs,stack,regular);
      }

      stochastic_distribution cumulative_distribution(variable_list(),sort_real::real_one());

      /* The conditions are collected for use. The stochastic operators before the action are ignored */
      while ((is_if_then(summandterm)||is_stochastic_operator(summandterm)))
      {
        if (is_if_then(summandterm))
        {
          const data_expression& localcondition=down_cast<if_then>(summandterm).condition();
          if (!(regular && singlestate))
          {
            condition1=lazy::and_(
                         condition1,
                         ((regular)?localcondition:
                          adapt_term_to_stack(
                            localcondition,
                            stack,
                            sumvars,
                            cumulative_distribution.variables())));
          }
          else
          {
            /* regular and singlestate */
            condition1=lazy::and_(localcondition,condition1);
          }
          summandterm=if_then(summandterm).then_case();
        }
        else
        {
          const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(summandterm);
          cumulative_distribution=stochastic_distribution(
                                    cumulative_distribution.variables()+sto.variables(),
                                    real_times_optimized(cumulative_distribution.distribution(),
                                                         sto.distribution()));
          summandterm=sto.operand();
        }
      }

      if (!cumulative_distribution.variables().empty() && !sumvars.empty())
      {
        throw mcrl2::runtime_error("Cannot permute sum operator and stochastic distribution in " + data::pp(summandterm));
      }

      if (is_seq(summandterm))
      {
        /* only one summand is needed */
        process_expression t1=seq(summandterm).left();
        process_expression t2=seq(summandterm).right();
        if (is_at(t1))
        {
          has_time=true;
          atTime=at(t1).time_stamp();
          t1=at(t1).operand();
        }

        if (t1==delta())
        {
          is_delta_summand=true;
        }
        else
        {
          assert(is_tau(t1)||is_action(t1)||is_sync(t1));
          multiAction=to_sorted_action_list(t1);
        }

        stochastic_distribution distribution(variable_list(),sort_real::real_one());
        process_expression t2_new=t2;
        if (is_stochastic_operator(t2))
        {
          const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(t2);
          distribution=stochastic_distribution(sto.variables(),sto.distribution());
          t2_new=sto.operand();
        }
        const assignment_list procargs=make_procargs(t2_new,stack,pCRLprocs,sumvars,regular,singlestate,distribution.variables());

        if (!regular)
        {
          if (!is_delta_summand)
          {
            multiAction=adapt_multiaction_to_stack(
                          multiAction,stack,sumvars);
          }
          if (has_time)
          {
            atTime=adapt_term_to_stack(
                     atTime,stack,sumvars,variable_list());
          }
          distribution=stochastic_distribution(
                                  distribution.variables(),
                                  adapt_term_to_stack(distribution.distribution(),stack,sumvars,distribution.variables()));
        }
        insert_summand(action_summands,deadlock_summands,
                       sumvars,condition1,multiAction,
                       atTime,distribution,procargs,
                       has_time,is_delta_summand);
        return;
      }

      /* There is a single initial multiaction or deadlock, possibly timed*/
      if (is_at(summandterm))
      {
        atTime=at(summandterm).time_stamp();
        summandterm=at(summandterm).operand();
        has_time=true;
      }
      else
      {
        // do nothing
      }

      if (is_delta(summandterm))
      {
        is_delta_summand=true;
      }
      else if (is_tau(summandterm))
      {
        // multiAction is already empty.
      }
      else if (is_action(summandterm))
      {
        assert(multiAction.empty()); // so the result must be sorted.
        multiAction.push_front(action(summandterm));
      }
      else if (is_sync(summandterm))
      {
        multiAction=to_sorted_action_list(summandterm);
      }
      else
      {
        throw mcrl2::runtime_error("expected multiaction " + process::pp(summandterm) +".");
      }

      if (regular)
      {
        if (!is_delta_summand)
          /* As termination has been replaced by an explicit action terminated, followed
           * by delta, a single terminating action cannot exist for regular processes. */
        {
          throw mcrl2::runtime_error("terminating processes should not exist when using the regular flag");
        }
        insert_summand(action_summands,
                       deadlock_summands,
                       sumvars,
                       condition1,
                       multiAction,
                       atTime,
                       stochastic_distribution(variable_list(),sort_real::real_one()),
                       dummyparameterlist(stack,singlestate),
                       has_time,
                       is_delta_summand);
        return;
      }

      multiAction=adapt_multiaction_to_stack(multiAction,stack,sumvars);
      assignment_list procargs ({ assignment(stack.stackvar,application(stack.opns->pop,stack.stackvar)) });

      insert_summand(
                action_summands,
                deadlock_summands,
                sumvars,
                condition1,
                multiAction,
                atTime,
                stochastic_distribution(variable_list(),sort_real::real_one()),   // TODO: UNLIKELY THAT THIS IS CORRECT.
                procargs,
                has_time,
                is_delta_summand);

      return;
    }


    void collectsumlistterm(
      const process_identifier& procId,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const process_expression& body,
      const variable_list& pars,
      const stacklisttype& stack,
      const bool regular,
      const bool singlestate,
      const std::set < process_identifier >& pCRLprocs,
      mutable_indexed_substitution<>& parameter_renaming)
    {
      if (is_choice(body))
      {
        const process_expression& t1=down_cast<choice>(body).left();
        const process_expression& t2=down_cast<choice>(body).right();

        collectsumlistterm(procId,action_summands,deadlock_summands,t1,pars,stack,
                           regular,singlestate,pCRLprocs,parameter_renaming);
        collectsumlistterm(procId,action_summands,deadlock_summands,t2,pars,stack,
                           regular,singlestate,pCRLprocs,parameter_renaming);
        return;
      }
      if (is_stochastic_operator(body))
      {
        /* Remove leading stochastic operators. They will not become
           part of this summand */
        const stochastic_operator& sto=atermpp::down_cast<const stochastic_operator>(body);
        collectsumlistterm(procId,action_summands,deadlock_summands,sto.operand(),pars,stack,
                           regular,singlestate,pCRLprocs,parameter_renaming);
      }
      else
      {
        add_summands(procId,
                     action_summands,
                     deadlock_summands,
                     process::replace_variables(body,parameter_renaming), // This renaming renames variables with the same name but 
                                                                          // conflicting types, e.g., x:Bool and x:Nat. 
                     pCRLprocs,
                     stack,
                     regular,
                     singlestate,
                     pars);
      }
    }

    void collectsumlist(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const std::set < process_identifier >& pCRLprocs,
      const variable_list& pars,
      const stacklisttype& stack,
      bool regular,
      bool singlestate,
      mutable_indexed_substitution<>& parameter_renaming)
    {
      for (const process_identifier& p: pCRLprocs)
      {
        objectdatatype& object=objectIndex(p);

        collectsumlistterm(
          p,
          action_summands,
          deadlock_summands,
          object.processbody,
          pars,
          stack,
          regular,
          singlestate,
          pCRLprocs,
          parameter_renaming);
      }
    }


    /**************** Enumtype and enumeratedtype **********************************/

    class enumeratedtype
    {
      public:
        std::size_t size;
        sort_expression sortId;
        data_expression_list elementnames;
        function_symbol_list functions;

        enumeratedtype(const std::size_t n,
                       specification_basic_type& spec)
        {
          size=n;
          if (n==2)
          {
            sortId = sort_bool::bool_();
            elementnames = data_expression_list({ sort_bool::false_(), sort_bool::true_()});
          }
          else
          {
            //create new sort identifier
            basic_sort sort_id(spec.fresh_identifier_generator("Enum" + std::to_string(n)));
            sortId=sort_id;
            //create structured sort
            //  Enumi = struct en_i | ... | e0_i;
            structured_sort_constructor_list struct_conss;
            for (std::size_t j=0 ; (j<n) ; j++)
            {
              //create constructor declaration of the structured sort
              const identifier_string s=spec.fresh_identifier_generator("e" + std::to_string(j) + "_" + std::to_string(n));
              const structured_sort_constructor struct_cons(s,"");

              struct_conss.push_front(struct_cons);
            }
            structured_sort sort_struct(struct_conss);

            //add declaration of standard functions
            spec.data.add_alias(alias(sort_id, sort_struct));

            //store new declarations in return value w
            sortId = sort_id;
            const function_symbol_vector& constructors=spec.data.constructors(sort_id);
            assert(constructors.size()==n);
            elementnames = data::data_expression_list(constructors.begin(),constructors.end());
          }
        }

        enumeratedtype(const enumeratedtype& e)
        {
          size=e.size;
          sortId=e.sortId;
          elementnames=e.elementnames;
          functions=e.functions;
        }

        void operator=(const enumeratedtype& e)
        {
          size=e.size;
          sortId=e.sortId;
          elementnames=e.elementnames;
          functions=e.functions;
        }

        ~enumeratedtype() = default;
    };

    std::size_t create_enumeratedtype(const std::size_t n)
    {
      std::size_t w;

      for (w=0; ((w<enumeratedtypes.size()) && (enumeratedtypes[w].size!=n)); ++w) {};

      if (w==enumeratedtypes.size()) // There is no enumeratedtype of arity n.
      {
        enumeratedtypes.emplace_back(n, *this);
      }
      return w;
    }

    data::function_symbol find_case_function(std::size_t index, const sort_expression& sort) const
    {
      for (const data::function_symbol& f: enumeratedtypes[index].functions)
      {
        sort_expression_list domain = function_sort(f.sort()).domain();
        assert(domain.size() >= 2);
        if (*(++domain.begin())==sort)
        {
          return f;
        }
      };

      throw mcrl2::runtime_error("searching for a nonexisting case function on sort " + data::pp(sort) +".");
      return data::function_symbol();
    }

    void define_equations_for_case_function(
      const std::size_t index,
      const data::function_symbol& functionname,
      const sort_expression& sort)
    {
      variable_list vars;
      data_expression_list args;
      data_expression_list xxxterm;

      const variable v1=get_fresh_variable("x",sort);
      const std::size_t n=enumeratedtypes[index].size;
      for (std::size_t j=0; (j<n); j++)
      {
        const variable v=get_fresh_variable("y",sort);
        vars.push_front(v);
        args.push_front(v);
        xxxterm.push_front(v1);
      }

      /* I generate here an equation of the form
         C(e,x,x,x,...x)=x for a variable x. */
      const sort_expression& s=enumeratedtypes[index].sortId;
      const variable v=get_fresh_variable("e",s);

      // we add e in front of xxxterm; note that xxxterm is not used afterwards
      // anymore, so we don't need to create a temporary copy for it here.
      xxxterm.push_front(data_expression(v));
      data.add_equation(
        data_equation(
          variable_list({ v1, v }),
          application(functionname,xxxterm),
          v1));
      fresh_equation_added=true;

      variable_list auxvars=vars;

      const data_expression_list& elementnames=enumeratedtypes[index].elementnames;
      for (const data_expression& w: elementnames)
      {
        assert(auxvars.size()>0);
        args.push_front(w);
        data.add_equation(data_equation(
                          vars,
                          application(functionname,args),
                          auxvars.front()));
        //fresh_equation_added=true; already set before the loop, omitting here
        args.pop_front();
        auxvars.pop_front();
      }
    }

    void create_case_function_on_enumeratedtype(
      const sort_expression& sort,
      const std::size_t enumeratedtype_index)
    {
      assert(enumeratedtype_index<enumeratedtypes.size());
      /* first find out whether the function exists already, in which
         case nothing needs to be done */

      const function_symbol_list& functions=enumeratedtypes[enumeratedtype_index].functions;
      const function_symbol_list::const_iterator i = std::find_if(functions.begin(), functions.end(), [&sort](const data::function_symbol& w){
        const function_sort& w1sort(down_cast<function_sort>(w.sort()));
        assert(function_sort(w1sort).domain().size()>1);
        // w matches if the second sort of the case function equals sort
        return *(++(w1sort.domain().begin())) == sort;
      });
      if(i != functions.end()) {
        return;
      }
      /* The function does not exist;
         Create a new function of enumeratedtype e, on sort */

      if (enumeratedtypes[enumeratedtype_index].sortId==sort_bool::bool_())
      {
        /* take the if function on sort 'sort' */
        enumeratedtypes[enumeratedtype_index].functions.push_front(if_(sort));
        return;
      }
      // else
      sort_expression_list newsortlist;
      const std::size_t n=enumeratedtypes[enumeratedtype_index].size;
      for (std::size_t j=0; j<n ; ++j)
      {
        newsortlist.push_front(sort);
      }
      const sort_expression& sid=enumeratedtypes[enumeratedtype_index].sortId;
      newsortlist.push_front(sid);

      const function_sort newsort(newsortlist,sort);
      const data::function_symbol casefunction(
        fresh_identifier_generator("C" + std::to_string(n) + "_" +
                         (!is_basic_sort(newsort)?"":std::string(basic_sort(sort).name()))), newsort);
      data.add_mapping(casefunction);
      enumeratedtypes[enumeratedtype_index].functions.push_front(casefunction);

      define_equations_for_case_function(enumeratedtype_index,casefunction,sort);
      return;
    }

    class enumtype
    {
      public:
        // enumtypes are not supposed to be copied.
        enumtype(const enumtype&)=delete;
        enumtype& operator =(const enumtype& )=delete;

        std::size_t enumeratedtype_index;
        variable var;

        enumtype(std::size_t n,
                 const sort_expression_list& fsorts,
                 const sort_expression_list& gsorts,
                 specification_basic_type& spec)
        {
          enumeratedtype_index=spec.create_enumeratedtype(n);

          var=variable(spec.fresh_identifier_generator("e"),spec.enumeratedtypes[enumeratedtype_index].sortId);
          spec.insertvariable(var,true);

          for (const sort_expression& f: fsorts)
          {
            spec.create_case_function_on_enumeratedtype(f,enumeratedtype_index);
          }

          for (const sort_expression& f: gsorts)
          {
            spec.create_case_function_on_enumeratedtype(f,enumeratedtype_index);
          }

          spec.create_case_function_on_enumeratedtype(sort_bool::bool_(),enumeratedtype_index);

          if (spec.timeIsBeingUsed || spec.stochastic_operator_is_being_used)
          {
            spec.create_case_function_on_enumeratedtype(sort_real::real_(),enumeratedtype_index);
          }
        }

        ~enumtype() = default;
    };
    /************** Merge summands using enumerated type ***********************/

    /* The function below returns true if the variable var could be mapped
       on an existing variable v' in matchinglist. The pars and args form pair
       form a substitution that will be extended with the pair [var,v']. i
       It returns false if the variable is new.

       If var is added (and not mapped on some other variable in the matchinglist/aka v)
       it is checked whether var occurs in  v or in the process_parameters,
       in which case var is renamed to a fresh variable. The renaming is added
       to the substitution encoded in pars/args.
    */


    bool mergeoccursin(
      variable& var,
      const variable_list& v,
      variable_list& matchinglist,
      variable_list& pars,
      data_expression_list& args,
      const variable_list& process_parameters)
    {
      variable_list auxmatchinglist;

      /* First find out whether var:sort can be matched to a
         term in the matching list */

      /* first find out whether the variable occurs in the matching
         list, so, they can be joined */

      for (variable_list::const_iterator i=matchinglist.begin();
           i!=matchinglist.end(); ++i)
      {
        variable var1=*i;
        if (var.sort()==var1.sort())
        {
          /* sorts match, so, we join the variables */
          if (var!=var1)
          {
            pars.push_front(var);
            args.push_front(data_expression(var1));
          }
          /* copy remainder of matching list */
          for (++i ; i!=matchinglist.end(); ++i)
          {
            auxmatchinglist.push_front(*i);
          }
          matchinglist=reverse(auxmatchinglist);
          return true;
        }
        else
        {
          auxmatchinglist.push_front(var1);
        }
      }

      /* turn auxmatchinglist back in normal order, and put the result
         in *matchinglist */

      matchinglist=reverse(auxmatchinglist);

      /* in this case no matching argument has been found.
      So, we must find out whether *var is an allowed variable, not
      occuring in the variablelist v.
      But if so, we must replace it by a new one. */
      for (const variable& var1 : v)
      {
        if (var.name() == var1.name())
        {
          pars.push_front(var);
          var=get_fresh_variable(var.name(),var.sort());
          args.push_front(data_expression(var));
          return false;
        }
      }

      /* Check whether the variable occurs in the prcoess parameter list, in which case
         it also needs to be renamed */
      for (const variable& var1: process_parameters)
      {
        if (var.name()==var1.name())
        {
          pars.push_front(var);
          var=get_fresh_variable(var.name(),var.sort());
          args.push_front(data_expression(var));
          return false;
        }
      }

      return false;
    }

    static
    data_expression_list extend(const data_expression& c, const data_expression_list& cl)
    {
      return data_expression_list(cl.begin(),
                                  cl.end(),
                                  [&c](const data_expression& e)->data_expression { return lazy::and_(c,e);} );

    }

    data_expression variables_are_equal_to_default_values(const variable_list& vl)
    {
      data_expression result=sort_bool::true_();
      for(const variable& v: vl)
      {
        const data_expression unique=representative_generator_internal(v.sort(),false);
        result=lazy::and_(result,equal_to(v,unique));
      }
      return result;
    }


    data_expression_list extend_conditions(
      const variable& var,
      const data_expression_list& conditionlist)
    {
      try
      {
        const data_expression unique=representative_generator_internal(var.sort(),false);
        const data_expression newcondition=equal_to(var,unique);
        return extend(newcondition,conditionlist);
      }
      catch (mcrl2::runtime_error&)
      {
        // The representative generator failed to find a term of var.sort().
        // No condition is generated, meaning that var will be unrestrained. This
        // is correct, and as the var.sort() has no concrete value, this will most
        // likely not effect matters as state space generation.
        return conditionlist;
      }
    }


    data_expression transform_matching_list(const variable_list& matchinglist)
    {
      data_expression result=sort_bool::true_();
      for(const variable& v: matchinglist)
      {
        try
        {
          data_expression unique=representative_generator_internal(v.sort(),false);
          result=lazy::and_(result, equal_to(v,unique));
        }
        catch (mcrl2::runtime_error&)
        {
          // No representant for sort v.sort() could be found. No condition is added.
        }
      }
      return result;
    }


    data_expression_list addcondition(
      const variable_list& matchinglist,
      const data_expression_list& conditionlist)
    {
      data_expression_list result=conditionlist;
      result.push_front(transform_matching_list(matchinglist));
      return result;
    }

/* Join the variables of v1 to v2 and rename the variables in v1
 * if needed. The conditionlist gives conditions to restrain variables
 * that did not occur in the other list. renaming pars and args give
 * renamings to be applied if variables in v1 had to be renamed. It
 * is not allowed to rename to names already occurring in the parameter
 * list. */

   variable_list merge_var(
      const variable_list& v1,
      const variable_list& v2,
      std::vector < variable_list>& renamings_pars,
      std::vector < data_expression_list>& renamings_args,
      data_expression_list& conditionlist,
      const variable_list& process_parameters)
    {
      data_expression_list renamingargs;
      variable_list renamingpars;
      variable_list matchinglist=v2;

      /* If the sequence of sum variables is reversed,
       * the variables are merged in the same sequence for all
       * summands (due to a suggestion of Muck van Weerdenburg) */

      variable_list v1h=(v2.empty()?reverse(v1):v1);

      variable_list result=v2;
      for (const variable& v_: v1h)
      {
        variable v=v_;
        if (!mergeoccursin(v,v2,
                           matchinglist,renamingpars,renamingargs,process_parameters))
        {
          result.push_front(v);
          conditionlist=extend_conditions(v,conditionlist);
        }
      }
      conditionlist=addcondition(matchinglist,conditionlist);
      renamings_pars.push_back(renamingpars);
      renamings_args.push_back(renamingargs);
      return result;
    }

    variable_list make_binary_sums(
      std::size_t n,
      const sort_expression& enumtypename,
      data_expression& condition,
      const variable_list& tail)
    {
      variable_list result;
      assert(n>1);
      condition=sort_bool::true_();

      n=n-1;
      for (result=tail ; (n>0) ; n=n/2)
      {
        variable sumvar=get_fresh_variable("e",enumtypename);
        result.push_front(sumvar);
        if ((n % 2)==0)
        {
          condition=lazy::and_(sumvar,condition);
        }
        else
        {
          condition=lazy::or_(sumvar,condition);
        }
      }
      return result;
    }

    data_expression construct_binary_case_tree_rec(
      std::size_t n,
      const variable_list& sums,
      data_expression_list& terms,
      const sort_expression& termsort,
      const enumtype& e)
    {
      assert(!terms.empty());

      if (n<=0)
      {
        assert(!terms.empty());
        const data_expression t=terms.front();
        terms.pop_front();
        return t;
      }

      assert(!sums.empty());
      const variable& casevar=sums.front();

      const data_expression t=construct_binary_case_tree_rec(n / 2,sums.tail(),terms,termsort,e);

      if (terms.empty())
      {
        return t;
      }

      const data_expression t1=construct_binary_case_tree_rec(n / 2,sums.tail(),terms,termsort,e);

      if (t==t1)
      {
        return t;
      }
      return application(find_case_function(e.enumeratedtype_index, termsort), casevar, t, t1);
    }

    template <class T>
    bool all_equal(const atermpp::term_list<T>& l)
    {
      if (l.empty())
      {
        return true;
      }
      typename atermpp::term_list<T>::const_iterator i=l.begin();
      const T& first=*i;
      for(++i ; i!=l.end(); ++i)
      {
        if (*i!=first)
        {
          return false;
        }
      }
      return true;
    }

    data_expression construct_binary_case_tree(
      std::size_t n,
      const variable_list& sums,
      data_expression_list terms,
      const sort_expression& termsort,
      const enumtype& e)
    {
      return construct_binary_case_tree_rec(n-1,sums,terms,termsort,e);
    }

    static
    bool summandsCanBeClustered(
      const stochastic_action_summand& summand1,
      const stochastic_action_summand& summand2)
    {
      if (summand1.has_time()!= summand2.has_time())
      {
        return false;
      }

      /* Here the multiactions are proper multi actions,
         both with or without a time indication */

      /* The actions can be clustered if they contain
         the same actions, with the same sorts for the
         actions. We assume that the multiactions are
         ordered.
      */

      const action_list multiactionlist1=summand1.multi_action().actions();
      const action_list multiactionlist2=summand2.multi_action().actions();
      action_list::const_iterator i2=multiactionlist2.begin();
      for (action_list::const_iterator i1=multiactionlist1.begin(); i1!=multiactionlist1.end(); ++i1,++i2)
      {
        if (i2==multiactionlist2.end())
        {
          return false;
        }
        if (i1->label()!=i2->label())
        {
          return false;
        }
      }
      return i2 == multiactionlist2.end();
    }

    static
    data_expression getRHSassignment(const variable& var, const assignment_list& as)
    {
      for (const assignment& a: as)
      {
        if (a.lhs()==var)
        {
          return a.rhs();
        }
      }
      return  var;
    }

    stochastic_action_summand collect_sum_arg_arg_cond(
      const enumtype& e,
      std::size_t n,
      const stochastic_action_summand_vector& action_summands,
      const variable_list& parameters)
    {
      /* This function gets a list of summands, with
         the same multiaction and time
         status. It yields a single clustered summand
         by introducing an auxiliary sum operator, with
         a variable of enumtype. In case binary is used,
         a sequence of variables are introduced of sort Bool */

      variable_list resultsum;
      data_expression resultcondition;
      action_list resultmultiaction;
      data_expression resulttime;

      std::vector < variable_list > rename_list_pars;
      std::vector < data_expression_list > rename_list_args;
      /* rename list is a list of pairs of variable and term lists */
      data_expression_list conditionlist;
      data_expression binarysumcondition;
      int equaluptillnow=1;

      std::set<variable> all_sum_variables;
      for (const stochastic_action_summand& smmnd: action_summands)
      {
        const variable_list sumvars=smmnd.summation_variables();
        resultsum=merge_var(sumvars,resultsum,rename_list_pars,rename_list_args,conditionlist,parameters);
        all_sum_variables.insert(sumvars.begin(),sumvars.end());
      }

      if (options.binary)
      {
        resultsum=make_binary_sums(
                    n,
                    enumeratedtypes[e.enumeratedtype_index].sortId,
                    binarysumcondition,
                    resultsum);
      }
      else
      {
        resultsum.push_front(e.var);
      }

      /* we construct the resulting condition */
      data_expression_list auxresult;
      std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
      std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();

      data_expression equalterm;
      equaluptillnow=1;
      for (const stochastic_action_summand& smmnd: action_summands)
      {
        const data_expression& condition=smmnd.condition();
        assert(auxrename_list_pars!=rename_list_pars.end());
        assert(auxrename_list_args!=rename_list_args.end());
        const variable_list auxpars= *auxrename_list_pars;
        ++auxrename_list_pars;
        const data_expression_list auxargs= *auxrename_list_args;
        ++auxrename_list_args;
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        data_expression_list::const_iterator j=auxargs.begin();
        for (variable_list::const_iterator i=auxpars.begin();
             i!=auxpars.end(); ++i, ++j)
        {
          /* Substitutions are carried out from left to right. The first applicable substitution counts */
          if (sigma(*i)==*i)  // sigma *i is undefined.
          {
            sigma[*i]=*j;
          }
        }
        maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);

        const data_expression auxresult1=replace_variables_capture_avoiding_alt(condition,mutable_sigma);
        if (equalterm==data_expression()||is_global_variable(equalterm))
        {
          equalterm=auxresult1;
        }
        else
        {
          if (equaluptillnow)
          {
            equaluptillnow=((auxresult1==equalterm)||is_global_variable(auxresult1));
          }
        }
        auxresult.push_front(auxresult1);
      }
      if (options.binary)
      {
        resultcondition=construct_binary_case_tree(n,
                        resultsum,auxresult,sort_bool::bool_(),e);
        resultcondition=lazy::and_(binarysumcondition,resultcondition);
        resultcondition=lazy::and_(
                          construct_binary_case_tree(n,
                              resultsum,conditionlist,sort_bool::bool_(),e),
                          resultcondition);
      }
      else
      {
        if (equaluptillnow)
        {
          if (all_equal(conditionlist))
          {
            resultcondition=lazy::and_(conditionlist.front(), equalterm);
          }
          else
          {
            data_expression_list tempconditionlist=conditionlist;
            tempconditionlist.push_front(data_expression(e.var));
            resultcondition=lazy::and_(
                              application(
                                find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                                tempconditionlist),
                              equalterm);
          }
        }
        else
        {
          data_expression_list tempauxresult=auxresult;
          tempauxresult.push_front(data_expression(e.var));
          resultcondition=application(
                            find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                            tempauxresult);
          if (all_equal(conditionlist))
          {
           resultcondition=lazy::and_(conditionlist.front(),resultcondition);
          }
          else
          {
           data_expression_list tempconditionlist=conditionlist;
           tempconditionlist.push_front(data_expression(e.var));
           resultcondition=lazy::and_(
                            application(
                              find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                              tempconditionlist),
                            resultcondition);
          }
        }
      }

      /* now we construct the arguments of the action */
      /* First we collect all multi-actions in a separate vector
         of multiactions */
      std::vector < action_list > multiActionList;

      for (const stochastic_action_summand& smmnd: action_summands)
      {
        multiActionList.push_back(smmnd.multi_action().actions());
      }

      action_list resultmultiactionlist;
      std::size_t multiactioncount= multiActionList[0].size(); // The number of multi actions.
      for (; multiactioncount>0 ; multiactioncount--)
      {
        data_expression_list resultf;
        // fcnt is the arity of the action with index multiactioncount-1;
        // const action a= *(multiActionList[0].begin()+(multiactioncount-1));
        action_list::const_iterator a=multiActionList[0].begin();
        for (std::size_t i=1 ; i<multiactioncount ; ++i,++a) {}
        // const action a= *((multiActionList[0]).begin()+(multiactioncount-1));
        std::size_t fcnt=(a->arguments()).size();
        data_expression f;

        for (; fcnt>0 ; fcnt--)
        {
          data_expression_list auxresult;
          data_expression equalterm;
          bool equaluptillnow=true;
          std::vector < variable_list >  ::const_iterator auxrename_list_pars=rename_list_pars.begin();
          std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();
          std::vector<action_list>::const_iterator  multiactionwalker=multiActionList.begin();
          for (stochastic_action_summand_vector::const_iterator walker=action_summands.begin(); walker!=action_summands.end();
               ++walker,++multiactionwalker)
          {
            assert(auxrename_list_pars!=rename_list_pars.end());
            assert(auxrename_list_args!=rename_list_args.end());
            const variable_list auxpars= *auxrename_list_pars;
            ++auxrename_list_pars;
            const data_expression_list auxargs= *auxrename_list_args;
            ++auxrename_list_args;
            // f is the fcnt-th argument of the multiactioncount-th action in the list
            action_list::const_iterator a1=multiactionwalker->begin();
            for (std::size_t i=1; i<multiactioncount; ++i, ++a1) {};
            data_expression_list::const_iterator d1=(a1->arguments()).begin();
            for (std::size_t i=1; i<fcnt; ++i, ++d1) {};
            f= *d1;
            maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
            data_expression_list::const_iterator j=auxargs.begin();
            for (variable_list::const_iterator i=auxpars.begin();
                 i!=auxpars.end(); ++i, ++j)
            {
              /* Substitutions are carried out from left to right. The first applicable substitution counts */
              if (sigma(*i)==*i)  // *i is not assigned in sigma.
              {
                sigma[*i]=*j;
              }
            }

            maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
            const data_expression auxresult1=replace_variables_capture_avoiding_alt(f,mutable_sigma);

            if (equalterm==data_expression()||is_global_variable(equalterm))
            {
              equalterm=auxresult1;
            }
            else
            {
              if (equaluptillnow)
              {
                equaluptillnow=((equalterm==auxresult1)||is_global_variable(auxresult1));
              }
            }
            auxresult.push_front(auxresult1);
          }
          if (equaluptillnow)
          {
            resultf.push_front(equalterm);
          }
          else
          {
            if (!options.binary)
            {
              data_expression_list tempauxresult=auxresult;
              tempauxresult.push_front(data_expression(e.var));
              resultf.push_front(data_expression(application(
                                                   find_case_function(e.enumeratedtype_index,f.sort()),
                                                   tempauxresult)));
            }
            else
            {
              data_expression temp=construct_binary_case_tree(
                                     n,
                                     resultsum,
                                     auxresult,
                                     f.sort(),
                                     e);
              resultf.push_front(temp);
            }
          }
        }
        a=multiActionList[0].begin();
        for (std::size_t i=1 ; i<multiactioncount ; ++i,++a) {}
        resultmultiactionlist.push_front(action(a->label(),resultf));
      }

      /* Construct resulttime, the time of the action ... */

      equaluptillnow=true;
      equalterm=data_expression();
      bool some_summand_has_time=false;
      bool all_summands_have_time=true;

      // first find out whether there is a summand with explicit time.
      for (const stochastic_action_summand& smmnd: action_summands)
      {
        if (smmnd.has_time())
        {
          some_summand_has_time=true;
        }
        else
        {
          all_summands_have_time=false;
        }
      }

      if (some_summand_has_time)
      {
        variable dummy_time_variable;
        if (!all_summands_have_time)
        {
          // Generate a fresh dummy variable, and add it to the summand variables
          dummy_time_variable=get_fresh_variable("dt",sort_real::real_());
          resultsum.push_front(dummy_time_variable);
        }
        std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();
        data_expression_list auxresult;
        for (const stochastic_action_summand& smmnd: action_summands)
        {
          if (smmnd.has_time())
          {
            const data_expression& actiontime=smmnd.multi_action().time();

            assert(auxrename_list_pars!=rename_list_pars.end());
            assert(auxrename_list_args!=rename_list_args.end());
            const variable_list auxpars= *auxrename_list_pars;
            ++auxrename_list_pars;
            const data_expression_list auxargs= *auxrename_list_args;
            ++auxrename_list_args;

            maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
            data_expression_list::const_iterator j=auxargs.begin();
            for (variable_list::const_iterator i=auxpars.begin();
                 i!=auxpars.end(); ++i, ++j)
            {
              /* Substitutions are carried out from left to right. The first applicable substitution counts */
              if (sigma(*i)==*i)   // sigma is not defined for *i.
              {
                sigma[*i]=*j;
              }
            }

            maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
            const data_expression auxresult1=replace_variables_capture_avoiding_alt(actiontime, mutable_sigma);
            if (equalterm==data_expression()||is_global_variable(equalterm))
            {
              equalterm=auxresult1;
            }
            else
            {
              if (equaluptillnow)
              {
                equaluptillnow=((auxresult1==equalterm)||is_global_variable(auxresult1));
              }
            }
            auxresult.push_front(auxresult1);
          }
          else
          {
            // this summand does not have time. But some summands have.
            auxresult.push_front(data_expression(dummy_time_variable));
            equaluptillnow=false;
          }
        }
        if (options.binary==1)
        {
          resulttime=construct_binary_case_tree(n,
                                                resultsum,auxresult,sort_real::real_(),e);
        }
        else
        {
          if (equaluptillnow)
          {
            resulttime=equalterm;
          }
          else
          {
            data_expression_list tempauxresult=auxresult;
            tempauxresult.push_front(data_expression(e.var));
            resulttime=application(
                         find_case_function(e.enumeratedtype_index,sort_real::real_()),
                         tempauxresult);
          }
        }
      }

      /* we construct the resulting distribution --------------------------------------------------------- */


      variable_list resulting_stochastic_variables;
      std::vector < variable_list > stochastic_rename_list_pars;
      std::vector < data_expression_list > stochastic_rename_list_args;
      data_expression resulting_distribution=sort_real::real_one();
      {
        data_expression_list stochastic_conditionlist;
        for (const stochastic_action_summand& smmnd: action_summands)
        {
          const variable_list stochastic_vars=smmnd.distribution().variables();
          resulting_stochastic_variables=merge_var(stochastic_vars,resulting_stochastic_variables,
                              stochastic_rename_list_pars,stochastic_rename_list_args,stochastic_conditionlist,
                              parameters + variable_list(all_sum_variables.begin(),all_sum_variables.end()));
        }

        std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();

        std::vector < variable_list >::const_iterator stochastic_auxrename_list_pars=stochastic_rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator stochastic_auxrename_list_args=stochastic_rename_list_args.begin();

        data_expression_list aux_result;
        data_expression equalterm;
        equaluptillnow=1;
        for (const stochastic_action_summand& smmnd: action_summands)
        {
          const data_expression& dist=smmnd.distribution().distribution();
          const variable_list stochastic_variables=smmnd.distribution().variables();
          assert(auxrename_list_pars!=rename_list_pars.end());
          assert(auxrename_list_args!=rename_list_args.end());
          assert(stochastic_auxrename_list_pars!=stochastic_rename_list_pars.end());
          assert(stochastic_auxrename_list_args!=stochastic_rename_list_args.end());

          const variable_list auxpars= *auxrename_list_pars;
          ++auxrename_list_pars;
          const data_expression_list auxargs= *auxrename_list_args;
          ++auxrename_list_args;
          const variable_list stochastic_auxpars= *stochastic_auxrename_list_pars;
          ++stochastic_auxrename_list_pars;
          const data_expression_list stochastic_auxargs= *stochastic_auxrename_list_args;
          ++stochastic_auxrename_list_args;
          maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
          data_expression_list::const_iterator j=stochastic_auxargs.begin();
          for (variable_list::const_iterator i=stochastic_auxpars.begin();
               i!=stochastic_auxpars.end(); ++i, ++j)
          {
            /* Substitutions are carried out from left to right. The first applicable substitution counts */
            if (sigma(*i)==*i)  // sigma is not defined for *i.
            {
              sigma[*i]=*j;
            }
          }
          j=auxargs.begin();
          for (variable_list::const_iterator i=auxpars.begin();
               i!=auxpars.end(); ++i, ++j)
          {
            /* Substitutions are carried out from left to right. The first applicable substitution counts */
            if (sigma(*i)==*i) // sigma is not defined for *i.
            {
              sigma[*i]=*j;
            }
          }
          maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);

          const data_expression auxresult1=replace_variables_capture_avoiding_alt(dist,mutable_sigma);
          if (equalterm==data_expression()||is_global_variable(equalterm))
          {
            equalterm=auxresult1;
          }
          else
          {
            if (equaluptillnow)
            {
              equaluptillnow=((auxresult1==equalterm)||is_global_variable(auxresult1));
            }
          }
          aux_result.push_front(auxresult1);
        }
        if (options.binary)
        {
          resulting_distribution=construct_binary_case_tree(n,
                          resultsum,aux_result,sort_real::real_(),e);
          resulting_distribution=lazy::and_(
                            construct_binary_case_tree(n,
                                resultsum,stochastic_conditionlist,sort_bool::bool_(),e),
                            resulting_distribution);
        }
        else
        {
          if (equaluptillnow)
          {
            if (all_equal(stochastic_conditionlist))
            {
              if (stochastic_conditionlist.front()==sort_bool::true_())
              {
                resulting_distribution=equalterm;
              }
              else
              {
                resulting_distribution=real_times_optimized(
                                                  if_(stochastic_conditionlist.front(),sort_real::real_one(),sort_real::real_zero()),
                                                  equalterm);
              }
            }
            else
            {
              data_expression_list temp_stochastic_conditionlist=stochastic_conditionlist;
              temp_stochastic_conditionlist.push_front(data_expression(e.var));
              resulting_distribution=real_times_optimized(
                                if_(application(
                                      find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                                      temp_stochastic_conditionlist),
                                      sort_real::real_one(),
                                      sort_real::real_zero()),
                                equalterm);
            }
          }
          else
          {
            data_expression_list tempauxresult=aux_result;
            tempauxresult.push_front(data_expression(e.var));
            resulting_distribution=application(
                              find_case_function(e.enumeratedtype_index,sort_real::real_()),
                              tempauxresult);
            if (all_equal(stochastic_conditionlist))
            {
              if (stochastic_conditionlist.front()==sort_bool::true_())
              {
               resulting_distribution=real_times_optimized(
                                                 if_(stochastic_conditionlist.front(),sort_real::real_one(),sort_real::real_zero()),
                                                 resulting_distribution);
              }
              else
              {
               resulting_distribution=real_times_optimized(
                                                 if_(stochastic_conditionlist.front(),sort_real::real_one(),sort_real::real_zero()),
                                                 resulting_distribution);
              }
            }
            else
            {
             data_expression_list temp_stochastic_conditionlist=stochastic_conditionlist;
             temp_stochastic_conditionlist.push_front(data_expression(e.var));
             resulting_distribution=real_times_optimized(
                              if_(application(
                                     find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                                     temp_stochastic_conditionlist),
                                  sort_real::real_one(),
                                  sort_real::real_zero()),
                              resulting_distribution);
            }
          }
        }
      }
      /* now we construct the arguments of the invoked -----------------------------------------------------
         process, i.e. the new function g */
      data_expression_list resultnextstate;

      for (const variable& var: parameters)
      {
        equalterm=data_expression();
        equaluptillnow=1;
        std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();
        std::vector < variable_list >::const_iterator stochastic_auxrename_list_pars=stochastic_rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator stochastic_auxrename_list_args=stochastic_rename_list_args.begin();
        data_expression_list auxresult;
        for (const stochastic_action_summand& smmnd: action_summands)
        {
          const assignment_list nextstate=smmnd.assignments();
          assert(auxrename_list_pars!=rename_list_pars.end());
          assert(auxrename_list_args!=rename_list_args.end());
          assert(stochastic_auxrename_list_pars!=stochastic_rename_list_pars.end());
          assert(stochastic_auxrename_list_args!=stochastic_rename_list_args.end());
          const variable_list auxpars= *auxrename_list_pars;
          ++auxrename_list_pars;
          const data_expression_list auxargs= *auxrename_list_args;
          ++auxrename_list_args;
          const variable_list stochastic_auxpars= *stochastic_auxrename_list_pars;
          ++stochastic_auxrename_list_pars;
          const data_expression_list stochastic_auxargs= *stochastic_auxrename_list_args;
          ++stochastic_auxrename_list_args;

          data_expression nextstateparameter;
          nextstateparameter=getRHSassignment(var,nextstate);

          maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
          data_expression_list::const_iterator j=stochastic_auxargs.begin();
          for (variable_list::const_iterator i=stochastic_auxpars.begin();
                 i!=stochastic_auxpars.end(); ++i, ++j)
          {
            /* Substitutions are carried out from left to right. The first applicable substitution counts */
            if (sigma(*i)==*i)  // sigma is not defined for *i.
            {
              sigma[*i]=*j;
            }
          }
          j=auxargs.begin();
          for (variable_list::const_iterator i=auxpars.begin();
                 i!=auxpars.end(); ++i, ++j)
          {
            /* Substitutions are carried out from left to right. The first applicable substitution counts */
            if (sigma(*i)==*i)  // sigma is not defined for *i.
            {
              sigma[*i]=*j;
            }
          }

          maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
          data_expression auxresult1=replace_variables_capture_avoiding_alt(nextstateparameter, mutable_sigma);
          if (equalterm==data_expression()) // ||is_global_variable(equalterm)) Adding this last part leads to smaller case functions,
                                            // but bigger and less structured state spaces, as parameters are less often put to default
                                            // values. Constant elimination can less often be applied as constant eliminiation does not
                                            // have subtle knowledge of case functions.
          {
            equalterm=auxresult1;
          }
          else if (equaluptillnow)
          {
            // set equaluptillnow if the arguments of this case function are all equal,
            // or are all equal to global variables. Setting a case function
            // C(e,v,dc1) to the value v, where dc1 is a global variable can result in
            // the growth of the state space, as dc1 is not set to a default value, but
            // keeps the value v.
            equaluptillnow=((equalterm==auxresult1)||
                               ((equalterm==data_expression()||is_global_variable(equalterm)) &&
                                                    is_global_variable(auxresult1)));
          }

          auxresult.push_front(auxresult1);
        }
        if (equaluptillnow)
        {
          resultnextstate.push_front(equalterm);
        }
        else
        {
          if (!options.binary)
          {
            data_expression_list tempauxresult=auxresult;
            tempauxresult.push_front(data_expression(e.var));
            resultnextstate.push_front(
                         data_expression(application(
                                           find_case_function(
                                             e.enumeratedtype_index,
                                             var.sort()),
                                           tempauxresult)));
          }
          else
          {
            data_expression temp=construct_binary_case_tree(
                                   n,
                                   resultsum,
                                   auxresult,
                                   var.sort(),
                                   e);
            resultnextstate.push_front(temp);
          }
        }
      }
      /* Now turn *resultg around */
      resultnextstate=reverse(resultnextstate);
      /* The list of arguments in nextstate are now in a sequential form, and
           must be transformed back to a list of assignments */
      const assignment_list final_resultnextstate=make_assignment_list(parameters,resultnextstate);
      return stochastic_action_summand(
                            resultsum,
                            resultcondition,
                            some_summand_has_time?multi_action(resultmultiactionlist,resulttime):multi_action(resultmultiactionlist),
                            final_resultnextstate,
                            stochastic_distribution(resulting_stochastic_variables,resulting_distribution));
    }

    deadlock_summand collect_sum_arg_arg_cond(
      const enumtype& e,
      std::size_t n,
      const deadlock_summand_vector& deadlock_summands,
      const variable_list& parameters)
    {
      /* This function gets a list of summands, with
         the same multiaction and time
         status. It yields a single clustered summand
         by introducing an auxiliary sum operator, with
         a variable of enumtype. In case binary is used,
         a sequence of variables are introduced of sort Bool */

      variable_list resultsum;
      data_expression resultcondition;
      action_list resultmultiaction;
      data_expression resulttime;

      std::vector < variable_list > rename_list_pars;
      std::vector < data_expression_list > rename_list_args;
      /* rename list is a list of pairs of variable and term lists */
      data_expression_list conditionlist;
      data_expression binarysumcondition;
      int equaluptillnow=1;

      for (const deadlock_summand& smmnd: deadlock_summands)
      {
        const variable_list sumvars=smmnd.summation_variables();
        resultsum=merge_var(sumvars,resultsum,rename_list_pars,rename_list_args,conditionlist,parameters);
      }

      if (options.binary)
      {
        resultsum=make_binary_sums(
                    n,
                    enumeratedtypes[e.enumeratedtype_index].sortId,
                    binarysumcondition,
                    resultsum);
      }
      else
      {
        resultsum.push_front(e.var);
      }

      /* we construct the resulting condition */
      data_expression_list auxresult;
      std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
      std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();

      data_expression equalterm;
      equaluptillnow=1;
      for (const deadlock_summand& smmnd: deadlock_summands)
      {
        const data_expression& condition=smmnd.condition();
        assert(auxrename_list_pars!=rename_list_pars.end());
        assert(auxrename_list_args!=rename_list_args.end());
        const variable_list auxpars= *auxrename_list_pars;
        ++auxrename_list_pars;
        const data_expression_list auxargs= *auxrename_list_args;
        ++auxrename_list_args;
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        data_expression_list::const_iterator j=auxargs.begin();
        for (variable_list::const_iterator i=auxpars.begin();
             i!=auxpars.end(); ++i, ++j)
        {
          /* Substitutions are carried out from left to right. The first applicable substitution counts */
          if (sigma(*i)== *i) // sigma is not defined for *i.
          {
            sigma[*i]=*j;
          }
        }

        maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
        const data_expression auxresult1=replace_variables_capture_avoiding_alt(condition, mutable_sigma);
        if (equalterm==data_expression()||is_global_variable(equalterm))
        {
          equalterm=auxresult1;
        }
        else
        {
          if (equaluptillnow)
          {
            equaluptillnow=((auxresult1==equalterm)||is_global_variable(auxresult1));
          }
        }
        auxresult.push_front(auxresult1);
      }
      if (options.binary)
      {
        resultcondition=construct_binary_case_tree(n,
                        resultsum,auxresult,sort_bool::bool_(),e);
        resultcondition=lazy::and_(binarysumcondition,resultcondition);
        resultcondition=lazy::and_(
                          construct_binary_case_tree(n,
                              resultsum,conditionlist,sort_bool::bool_(),e),
                          resultcondition);
      }
      else
      {
        if (equaluptillnow)
        {
          if (all_equal(conditionlist))
          {
            resultcondition=lazy::and_(conditionlist.front(),equalterm);
          }
          else
          {
            data_expression_list tempconditionlist=conditionlist;
            tempconditionlist.push_front(data_expression(e.var));
            resultcondition=lazy::and_(
                            application(
                              find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                              tempconditionlist),
                            equalterm);
          }
        }
        else
        {
          data_expression_list tempauxresult=auxresult;
          tempauxresult.push_front(data_expression(e.var));
          resultcondition=application(
                            find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                            tempauxresult);
          if (all_equal(conditionlist))
          {
            resultcondition=lazy::and_(conditionlist.front(),resultcondition);
          }
          else
          {
            data_expression_list tempconditionlist=conditionlist;
            tempconditionlist.push_front(data_expression(e.var));
            resultcondition=lazy::and_(
                            application(
                              find_case_function(e.enumeratedtype_index,sort_bool::bool_()),
                              tempconditionlist),
                            resultcondition);
          }
        }
      }

      /* now we construct the arguments of the action */
      /* First we collect all multi-actions in a separate vector
         of multiactions */
      std::vector < action_list > multiActionList;

      action_list resultmultiactionlist;

      /* Construct resulttime, the time of the action ... */

      equaluptillnow=true;
      equalterm=data_expression();
      bool some_summand_has_time=false;
      bool all_summands_have_time=true;

      // first find out whether there is a summand with explicit time.
      for (const deadlock_summand& smmnd: deadlock_summands)
      {
        if (smmnd.deadlock().has_time())
        {
          some_summand_has_time=true;
        }
        else
        {
          all_summands_have_time=false;
        }
      }

      if ((some_summand_has_time))
      {
        variable dummy_time_variable;
        if (!all_summands_have_time)
        {
          // Generate a fresh dummy variable, and add it to the summand variables
          dummy_time_variable=get_fresh_variable("dt",sort_real::real_());
          resultsum.push_front(dummy_time_variable);
        }
        std::vector < variable_list >::const_iterator auxrename_list_pars=rename_list_pars.begin();
        std::vector < data_expression_list >::const_iterator auxrename_list_args=rename_list_args.begin();
        data_expression_list auxresult;
        for (const deadlock_summand& smmnd: deadlock_summands)
        {
          if (smmnd.deadlock().has_time())
          {
            const data_expression& actiontime=smmnd.deadlock().time();

            assert(auxrename_list_pars!=rename_list_pars.end());
            assert(auxrename_list_args!=rename_list_args.end());
            const variable_list auxpars= *auxrename_list_pars;
            ++auxrename_list_pars;
            const data_expression_list auxargs= *auxrename_list_args;
            ++auxrename_list_args;

            maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
            data_expression_list::const_iterator j=auxargs.begin();
            for (variable_list::const_iterator i=auxpars.begin();
                 i!=auxpars.end(); ++i, ++j)
            {
              /* Substitutions are carried out from left to right. The first applicable substitution counts */
              if (sigma(*i)==*i)  // sigma is not defined for *i.
              {
                sigma[*i]=*j;
              }
            }

            maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
            const data_expression auxresult1=replace_variables_capture_avoiding_alt(actiontime, mutable_sigma);
            if (equalterm==data_expression()||is_global_variable(equalterm))
            {
              equalterm=auxresult1;
            }
            else
            {
              if (equaluptillnow)
              {
                equaluptillnow=((auxresult1==equalterm)||is_global_variable(auxresult1));
              }
            }
            auxresult.push_front(auxresult1);
          }
          else
          {
            // this summand does not have time. But some summands have.
            auxresult.push_front(data_expression(dummy_time_variable));
            equaluptillnow=false;
          }
        }
        if (options.binary==1)
        {
          resulttime=construct_binary_case_tree(n,
                                                resultsum,auxresult,sort_real::real_(),e);
        }
        else
        {
          if (equaluptillnow)
          {
            resulttime=equalterm;
          }
          else
          {
            data_expression_list tempauxresult=auxresult;
            tempauxresult.push_front(data_expression(e.var));
            resulttime=application(
                         find_case_function(e.enumeratedtype_index,sort_real::real_()),
                         tempauxresult);
          }
        }
      }

      return deadlock_summand(resultsum,
                              resultcondition,
                              some_summand_has_time?deadlock(resulttime):deadlock());
    }

    static
    sort_expression_list getActionSorts(const action_list& actionlist)
    {
      sort_expression_list resultsorts;

      for (const action& act: actionlist)
      {
        resultsorts=act.label().sorts()+resultsorts;
      }
      return resultsorts;
    }

    void cluster_actions(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const variable_list& pars)
    {
      {
        /* We cluster first the action summands with the action
            occurring in the first summand of sums.
            These summands are first stored in w1. */

        stochastic_action_summand_vector result;
        stochastic_action_summand_vector reducible_sumlist=action_summands;

        for (const stochastic_action_summand& summand1: action_summands)
        {
          stochastic_action_summand_vector w1;
          stochastic_action_summand_vector w2;
          for (const stochastic_action_summand& summand2: reducible_sumlist)
          {
            if (summandsCanBeClustered(summand1,summand2))
            {
              w1.push_back(summand2);
            }
            else
            {
              w2.push_back(summand2);
            }
          }
          reducible_sumlist.swap(w2);

          /* In w1 we now find all the summands labelled with
             similar multiactions, actiontime and terminationstatus.
             We must now construct its clustered form. */
          std::size_t n=w1.size();

          if (n>0)
          {
            if (n>1)
            {
              const action_list multiaction=w1.front().multi_action().actions();
              sort_expression_list actionsorts;
              actionsorts=getActionSorts(multiaction);

              const enumtype enumeratedtype_(options.binary?2:n,actionsorts,get_sorts(pars),*this);

              result.push_back(collect_sum_arg_arg_cond(enumeratedtype_,n,w1,pars));
            }
            else
            {
              // result=w1 + result;
              for(const stochastic_action_summand& summand: result)
              {
                w1.push_back(summand);
              }
              w1.swap(result);
            }
          }
        }
        action_summands=result;
      }

      // Now the delta summands are clustered.
      deadlock_summand_vector result;
      deadlock_summand_vector reducible_sumlist=deadlock_summands;

      for (const deadlock_summand& summand1: deadlock_summands)
      {
        deadlock_summand_vector w1;
        deadlock_summand_vector w2;
        for (const deadlock_summand& summand2: reducible_sumlist)
        {
          if (summand1.deadlock().has_time()==summand2.deadlock().has_time())
          {
            w1.push_back(summand2);
          }
          else
          {
            w2.push_back(summand2);
          }
        }
        reducible_sumlist.swap(w2);

        /* In w1 we now find all the summands labelled with
           either no or an action time.
           We must now construct its clustered form. */
        const std::size_t n=w1.size();

        if (n>0)
        {
          if (n>1)
          {
            sort_expression_list actionsorts;

            const enumtype enumeratedtype_(options.binary?2:n,actionsorts,get_sorts(pars),*this);

            result.push_back(collect_sum_arg_arg_cond(enumeratedtype_,n,w1,pars));
          }
          else
          {
            assert(w1.size()==1);
            result.push_back(w1.front());
          }
        }
      }
      assert(reducible_sumlist.empty());
      deadlock_summands.swap(result);
    }


    /**************** GENERaTE LPEpCRL **********************************/


    /* The variable regular indicates that we are interested in generating
       a LPE assuming the pCRL term under consideration is regular */

    void generateLPEpCRL(stochastic_action_summand_vector& action_summands,
                         deadlock_summand_vector& deadlock_summands,
                         const process_identifier& procId,
                         const bool containstime,
                         const bool regular,
                         variable_list& parameters,
                         data_expression_list& init,
                         stochastic_distribution& initial_stochastic_distribution)
    /* A pair of initial state and linear process must be extracted
       from the underlying GNF */
    {
      // We use action_summands and deadlock_summands as an output.
      assert(action_summands.empty());
      assert(deadlock_summands.empty());

      bool singlecontrolstate=false;
      std::set< process_identifier > reachable_process_identifiers;
      make_pCRL_procs(procId, reachable_process_identifiers);

      /* now reachable process identifiers contains a list of all process id's in this pCRL process.
         Some of these processes have equal parameter lists and right hand sides. A typical example is
         P1 = a.P2, P3=a.P4, P2=b.Q, P4=b.Q. This reduces to P1 = a.P2, P2 = b.Q.
         We reduce the set of process_identifiers in reachable_process_identifiers and perform the
         appropriate substitution in the right hand side.
         The result is a map from process_identifiers to adapted process bodies.
      */
      const std::set< process_identifier > reduced_set_of_process_identifiers =
          minimize_set_of_reachable_process_identifiers(reachable_process_identifiers,procId);

      /* The equations can still contain stochastic operators that occur before the first action.
         We translate the equations such that the stochastic operators do not occur in front anymore.
         Moving the stochastic operators to the front means that the number of parameters of each
         process can increase. In this case we introduce a new process identifiers. The resulting
         ids of processes occur in stochastic_normalized_process_identifiers. */

      process_identifier initial_proc_id=procId;  // the initial process id may be renamed.
      const std::set< process_identifier > stochastic_normalized_process_identifiers =
                  remove_stochastic_operators_from_front(reduced_set_of_process_identifiers, initial_proc_id, initial_stochastic_distribution);

      /* collect the parameters, but assume that variables
         have a unique sort */
      if (stochastic_normalized_process_identifiers.size()==1)
      {
        singlecontrolstate=true;
      }
      mutable_indexed_substitution<> parameter_renaming;
      parameters=collectparameterlist(stochastic_normalized_process_identifiers,parameter_renaming);

      if ((!regular)||((!singlecontrolstate) && (options.newstate) && (!options.binary)))
      {
        declare_control_state(stochastic_normalized_process_identifiers);
      }
      stacklisttype stack(parameters,*this,regular,stochastic_normalized_process_identifiers,singlecontrolstate);

      if (regular)
      {
        if ((options.binary) && (options.newstate))
        {
          parameters=stack.parameters;
          if (!singlecontrolstate)
          {
            parameters=reverse(stack.booleanStateVariables) + parameters;
          }
        }
        else  /* !binary or oldstate */
        {
          variable_list tempparameters=stack.parameters;
          tempparameters.push_front(stack.stackvar);
          parameters=
            ((!singlecontrolstate)?tempparameters:stack.parameters);
        }
      }
      else /* not regular, use a stack */
      {
        parameters = variable_list({ stack.stackvar });
      }
      init=make_initialstate(initial_proc_id,stack,stochastic_normalized_process_identifiers,regular,singlecontrolstate,initial_stochastic_distribution);
      assert(init.size()==parameters.size());
      collectsumlist(action_summands,deadlock_summands,stochastic_normalized_process_identifiers,parameters,stack,regular,singlecontrolstate,parameter_renaming);

      if (!options.no_intermediate_cluster)
      {
        cluster_actions(action_summands,deadlock_summands,parameters);
      }

      if ((!containstime) || options.ignore_time)
      {
        /* We add a delta summand to each process, if the flag
           ignore_time is set, or if the process does not contain any
           explicit reference to time. This affects the behaviour of each
           process in the sense that each process can idle
           indefinitely. It has the advantage that large numbers
           numbers of timed delta summands are subsumed by
           this delta. As the removal of timed delta's
           is time consuming in the linearisation, the use
           of this flag, can speed up linearisation considerably */
        deadlock_summands.emplace_back(variable_list(), sort_bool::true_(), deadlock());
      }
    }


    /**************** hiding *****************************************/

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
    static
    std::string log_hide_application(const lps_statistics_t& lps_statistics_before,
                                     const lps_statistics_t& lps_statistics_after,
                                     const std::size_t num_hidden_actions,
                                     size_t indent = 0)
    {
      std::string indent_str(indent, ' ');
      std::ostringstream os;

      os << indent_str << "- operator: hide" << std::endl;

      indent += 2;
      indent_str = std::string(indent, ' ');
      os << indent_str << "number of hidden actions: " << num_hidden_actions << std::endl
         << indent_str << "before:" << std::endl << print(lps_statistics_before, indent+2)
         << indent_str << "after:" << std::endl << print(lps_statistics_after, indent+2);

      return os.str();
    }
#endif // MCRL2_LOG_LPS_LINEARISE_STATISTICS

    static
    action_list hide_(const identifier_string_list& hidelist, const action_list& multiaction)
    {
      action_list resultactionlist;

      for (const action& a: multiaction)
      {
        if (std::find(hidelist.begin(),hidelist.end(),a.label().name())==hidelist.end())
        {
          resultactionlist.push_front(a);
        }
      }

      /* reverse the actionlist to maintain the ordering */
      resultactionlist = reverse(resultactionlist);
      return resultactionlist;
    }

    static
    void hidecomposition(const identifier_string_list& hidelist, stochastic_action_summand_vector& action_summands)
    {
#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_before = get_statistics(action_summands);
#endif
      for (stochastic_action_summand_vector::iterator i=action_summands.begin(); i!=action_summands.end() ; ++i)
      {
        const action_list acts=hide_(hidelist,i->multi_action().actions());
        *i=stochastic_action_summand(i->summation_variables(),
                          i->condition(),
                          i->has_time()?multi_action(acts,i->multi_action().time()):multi_action(acts),
                          i->assignments(),
                          i->distribution());
      }

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_after = get_statistics(action_summands);
      std::cout << log_hide_application(lps_statistics_before, lps_statistics_after, hidelist.size());
#endif
    }

    /**************** equalargs ****************************************/

    static
    bool occursinvarandremove(const variable& var, variable_list& vl)
    {
      bool result=false;

      if (vl.empty())
      {
        return false;
      }
      vl.pop_front();
      const variable var1=vl.front();
      if (var==var1)
      {
        return true;
      }

      // Names of variables cannot be the same, even if they have different types.
      if (var.name()==var1.name())
      {
        throw mcrl2::runtime_error("variable conflict " + data::pp(var) + ":" + data::pp(var.sort()) + " versus " +
                                   data::pp(var1) + ":" + data::pp(var1.sort()) + ".");
      }

      result=occursinvarandremove(var,vl);
      vl.push_front(var1);
      return result;
    }

    /********************** construct renaming **************************/

    variable_list construct_renaming(
      const variable_list& pars1,
      const variable_list& pars2,
      variable_list& pars3,
      variable_list& pars4,
      const bool unique=true)
    {
      /* check whether the variables in pars2 are unique,
         wrt to those in pars1, and deliver:
         - in pars3 a list of renamed parameters pars2, such that
           pars3 is unique with respect to pars1;
         - in pars4 a list of parameters that need to be renamed;
         - as a return result, new values for the parameters in pars4.
           This allows using substitute_data(list) to rename
           action and process arguments and conditions to adapt
           to the new parameter names.
         The variable unique indicates whether the generated variables
         are unique, and not occurring elsewhere. If unique is false,
         it is attempted to reuse previously generated variable names,
         as long as they do not occur in pars1. The default value for
         unique is true.
       */

      variable_list t;
      variable_list t1;
      variable_list t2;

      if (pars2.empty())
      {
        pars3=variable_list();
        pars4=variable_list();
      }
      else
      {
        const variable& var2=pars2.front();
        variable var3=var2;
        for (std::size_t i=0 ; occursin(var3,pars1)||occursin(var3,pars2) ; ++i)
        {
          var3=get_fresh_variable(var2.name(),var2.sort(),(unique?-1:i));
        }
        if (var3!=var2)
        {
          t1=construct_renaming(pars1,pars2.tail(),t,t2,unique);
          t1.push_front(var3);

          pars4=t2;
          pars4.push_front(var2);
          pars3=t;
          pars3.push_front(var3);
        }
        else
        {
          t1=construct_renaming(pars1,pars2.tail(),t,pars4,unique);
          pars3=t;
          pars3.push_front(var2);
        }

      }
      return t1;
    }

    /**************** communication operator composition ****************/

    template <typename List>
    static sort_expression_list get_sorts(const List& l)
    {
      return sort_expression_list(l.begin(), l.end(), [](const typename List::value_type& d) -> sort_expression {return d.sort();});
    }

    static
    bool check_real_variable_occurrence(
      const variable_list& sumvars,
      const data_expression& actiontime,
      const data_expression& condition)
    {
      /* Check whether actiontime is an expression
         of the form t1 +...+ tn, where one of the
         ti is a variable in sumvars that does not occur in condition */

      if (is_variable(actiontime))
      {
        const variable& t = atermpp::down_cast<variable>(actiontime);
        if (occursintermlist(t, variable_list_to_data_expression_list(sumvars)) && !occursinterm(condition, t))
        {
          return true;
        }
      }

      if (sort_real::is_plus_application(actiontime))
      {
        return (check_real_variable_occurrence(sumvars,data::binary_left(application(actiontime)),condition) ||
                check_real_variable_occurrence(sumvars,data::binary_right(application(actiontime)),condition));
      }

      return false;
    }

    data_expression makesingleultimatedelaycondition(
      const variable_list& sumvars,
      const variable_list& freevars,
      const data_expression& condition,
      const bool has_time,
      const variable& timevariable,
      const data_expression& actiontime,
      variable_list& used_sumvars)
    {
      /* Generate a condition of the form:

           exists sumvars. condition && timevariable<actiontime

         where the sumvars are added to the existentially quantified
         variables, and the resulting expression is

           condition && timevariable<actiontime

         The comments below refer to old code, where an explicit
         existential quantor was generated.

         OLD:
         Currently, the existential quantifier must use an equation,
         which represents a higher order function. The existential
         quantifier is namely of type exists:sorts1->Bool, where sorts1
         are the sorts of the quantified variables.

         If the sum variables do not occur in the expression, they
         are not quantified.

         If the actiontime is of the form t1+t2+...+tn where one
         of the ti is a quantified real variable in sumvars, and this
         variable does not occur in the condition, then the expression
         of the form timevariable < actiontime is omitted.
      */

      assert(used_sumvars.empty());
      data_expression result;
      variable_list variables;
      if (!has_time || (check_real_variable_occurrence(sumvars,actiontime,condition)))
      {
        result=condition;
      }
      else
      {
        result=RewriteTerm(
                 lazy::and_(
                   condition,
                   data::less(timevariable,actiontime)));
        variables.push_front(timevariable);
      }
      for (const variable& v: freevars)
      {
        if (occursinterm(result, v))
        {
          variables.push_front(v);
        }
      }

      for (const variable& v: global_variables)
      {
        if (occursinterm(result, v))
        {
          variables.push_front(v);
        }
      }

      for (const variable& v: sumvars)
      {
        if (occursinterm(result, v))
        {
          used_sumvars.push_front(v);
        }
      }
      used_sumvars = reverse(used_sumvars);

      return result;
    }

    lps::detail::ultimate_delay getUltimateDelayCondition(
      const stochastic_action_summand_vector& action_summands,
      const deadlock_summand_vector& deadlock_summands,
      const variable_list& freevars)
    {
      /* First walk through the summands to see whether
         a summand with condition true that does not refer
         to time exists. In that case the ultimate delay
         condition is true */

      variable time_variable=get_fresh_variable("timevar",sort_real::real_());
      for (const deadlock_summand& summand: deadlock_summands)
      {
        if ((!summand.deadlock().has_time()) && (summand.condition()==sort_bool::true_()))
        {
          return lps::detail::ultimate_delay(time_variable);
        }
      }

      for (const stochastic_action_summand& summand: action_summands)
      {
        if ((!summand.multi_action().has_time()) && (summand.condition()==sort_bool::true_()))
        {
          return lps::detail::ultimate_delay(time_variable);
        }
      }

      /* Unfortunately, no ultimate delay condition true can
         be generated. So, we must now traverse all conditions
         to generate a non trivial ultimate delay condition */

      data_expression_list results;
      data_expression_list condition_list;
      std::vector < variable_list> renamings_pars;
      std::vector < data_expression_list> renamings_args;
      variable_list existentially_quantified_variables;
      for (const deadlock_summand& s: deadlock_summands)
      {
        variable_list new_existentially_quantified_variables;
        const data_expression ult_del_condition=
             makesingleultimatedelaycondition(
                             s.summation_variables(),
                             freevars,
                             s.condition(),
                             s.deadlock().has_time(),
                             time_variable,
                             s.deadlock().time(),
                             new_existentially_quantified_variables);
        existentially_quantified_variables=merge_var(
                                             new_existentially_quantified_variables,
                                             existentially_quantified_variables,
                                             renamings_pars,
                                             renamings_args,
                                             condition_list,
                                             variable_list());
        results.push_front(ult_del_condition);
      }

      for (const stochastic_action_summand& s: action_summands)
      {
        variable_list new_existentially_quantified_variables;
        const data_expression ult_del_condition=
             makesingleultimatedelaycondition(
                             s.summation_variables(),
                             freevars,
                             s.condition(),
                             s.multi_action().has_time(),
                             time_variable,
                             s.multi_action().time(),
                             new_existentially_quantified_variables);
        existentially_quantified_variables=merge_var(
                                             new_existentially_quantified_variables,
                                             existentially_quantified_variables,
                                             renamings_pars,
                                             renamings_args,
                                             condition_list,
                                             variable_list());
        results.push_front(ult_del_condition);
      }

      data_expression result=sort_bool::false_();

      assert(results.size()==condition_list.size());
      assert(results.size()==renamings_pars.size());
      assert(results.size()==renamings_args.size());

      std::vector < variable_list>::const_iterator renamings_par=renamings_pars.begin();
      std::vector < data_expression_list>::const_iterator renamings_arg=renamings_args.begin();
      condition_list=reverse(condition_list);
      results=reverse(results);
      data_expression_list::const_iterator j=condition_list.begin();
      for(data_expression_list::const_iterator i=results.begin();
              i!=results.end(); ++i,++j,++renamings_par,++renamings_arg)
      {
        const variable_list& auxpars=*renamings_par;
        const data_expression_list& auxargs=*renamings_arg;

        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        data_expression_list::const_iterator j1=auxargs.begin();
        for (variable_list::const_iterator i1=auxpars.begin();
             i1!=auxpars.end(); ++i1, ++j1)
        {
          /* Substitutions are carried out from left to right. The first applicable substitution counts */
          if (sigma(*i1)==*i1)
          {
            sigma[*i1]=*j1;

          }
        }

        maintain_variables_in_rhs< mutable_map_substitution<> > mutable_sigma(sigma);
        result=lazy::or_(result,replace_variables_capture_avoiding_alt(lazy::and_(*i,*j), mutable_sigma));
      }
      return lps::detail::ultimate_delay(time_variable, existentially_quantified_variables, RewriteTerm(result));
    }


    /******** make_unique_variables **********************/

    data::maintain_variables_in_rhs< data::mutable_map_substitution<> >  make_unique_variables(
      const variable_list& var_list,
      const std::string& hint)
    {
      /* This function generates a list of variables with the same sorts
         as in variable_list, where all names are unique */

      data::maintain_variables_in_rhs< data::mutable_map_substitution<> >  sigma;

      for(const variable& var: var_list)
      {
        const data::variable v =
              get_fresh_variable(std::string(var.name()) + ((hint.empty())?"":"_") + hint, var.sort());
        sigma[var] = v;
      }
      return sigma;
    }

    /******** make_parameters_and_variables_unique **********************/

    void make_parameters_and_sum_variables_unique(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      variable_list& pars,
      lps::detail::ultimate_delay& ultimate_delay_condition,
      const std::string& hint="")
    {
      stochastic_action_summand_vector result_action_summands;

      data::maintain_variables_in_rhs<data::mutable_map_substitution<> > sigma=make_unique_variables(pars, hint);
      const variable_list unique_pars=data::replace_variables(pars, sigma);

      if (!options.ignore_time)
      {
        // Remove variables locally bound in the ultimate_delay_condition
        maintain_variables_in_rhs< data::mutable_map_substitution<> > local_sigma=sigma;
        for(const variable& v: ultimate_delay_condition.variables())
        {
          local_sigma[v]=v;
        }
        ultimate_delay_condition.constraint()=replace_variables_capture_avoiding_alt(
                                                       ultimate_delay_condition.constraint(),
                                                       local_sigma);  // Only substitute the variables in the lhs.
      }
      for (const stochastic_action_summand& smmnd: action_summands)
      {
        const variable_list& sumvars=smmnd.summation_variables();
        data::maintain_variables_in_rhs<data::mutable_map_substitution<> > sigma_sumvars=make_unique_variables(sumvars,hint);
        const variable_list unique_sumvars=data::replace_variables(sumvars, sigma_sumvars);

        stochastic_distribution distribution=smmnd.distribution();
        const variable_list stochastic_vars=distribution.variables();
        data::maintain_variables_in_rhs<data::mutable_map_substitution<> > sigma_stochastic_vars=
                                          make_unique_variables(stochastic_vars,hint);
        const variable_list unique_stochastic_vars=data::replace_variables(stochastic_vars, sigma_stochastic_vars);

        data_expression condition=smmnd.condition();
        action_list multiaction=smmnd.multi_action().actions();
        data_expression actiontime=smmnd.multi_action().time();
        assignment_list nextstate=smmnd.assignments();

        condition=replace_variables_capture_avoiding_alt(condition, sigma_sumvars);
        condition=replace_variables_capture_avoiding_alt(condition, sigma);

        actiontime=replace_variables_capture_avoiding_alt(actiontime, sigma_sumvars);
        actiontime=replace_variables_capture_avoiding_alt(actiontime, sigma);
        multiaction=lps::replace_variables_capture_avoiding_with_an_identifier_generator(multiaction, sigma_sumvars, fresh_identifier_generator);
        multiaction=lps::replace_variables_capture_avoiding_with_an_identifier_generator(multiaction, sigma, fresh_identifier_generator);

        distribution=stochastic_distribution(
                       unique_stochastic_vars,
                       replace_variables_capture_avoiding_alt(distribution.distribution(), sigma_sumvars));
        distribution=stochastic_distribution(
                       unique_stochastic_vars,
                       replace_variables_capture_avoiding_alt(distribution.distribution(), sigma_stochastic_vars));
        distribution=stochastic_distribution(
                       distribution.variables(),
                       replace_variables_capture_avoiding_alt(distribution.distribution(), sigma));

        nextstate=substitute_assignmentlist(nextstate,pars,false,true,sigma_sumvars);
        nextstate=substitute_assignmentlist(nextstate,pars,false,true,sigma_stochastic_vars);
        nextstate=substitute_assignmentlist(nextstate,pars,true,true,sigma);

        result_action_summands.emplace_back(unique_sumvars,
            condition,
            smmnd.multi_action().has_time() ? multi_action(multiaction, actiontime) : multi_action(multiaction),
            nextstate,
            distribution);
      }
      pars=unique_pars;
      action_summands.swap(result_action_summands);

      deadlock_summand_vector result_deadlock_summands;

      assert(unique_pars.size()==pars.size());

      for (const deadlock_summand& smmnd: deadlock_summands)
      {
        const variable_list& sumvars=smmnd.summation_variables();
        maintain_variables_in_rhs<data::mutable_map_substitution<> > sigma_sumvars=make_unique_variables(sumvars,hint);
        const variable_list unique_sumvars=data::replace_variables(sumvars, sigma_sumvars);

        assert(unique_sumvars.size()==sumvars.size());
        data_expression condition=smmnd.condition();
        data_expression actiontime=smmnd.deadlock().time();

        condition=replace_variables_capture_avoiding_alt(condition, sigma_sumvars);
        condition=replace_variables_capture_avoiding_alt(condition, sigma);

        actiontime=replace_variables_capture_avoiding_alt(actiontime, sigma_sumvars);
        actiontime=replace_variables_capture_avoiding_alt(actiontime, sigma);

        result_deadlock_summands.emplace_back(unique_sumvars,
            condition,
            smmnd.deadlock().has_time() ? deadlock(actiontime) : deadlock());
      }
      pars=unique_pars;
      result_deadlock_summands.swap(deadlock_summands);
    }





    /**************** parallel composition ******************************/



    /// \brief Returns the conjunction of the two delay conditions and the join of the variables, where
    ///        the variables in delay2 are renamed to avoid conflict with those in delay1.
    lps::detail::ultimate_delay combine_ultimate_delays(
                           const lps::detail::ultimate_delay& delay1,
                           const lps::detail::ultimate_delay& delay2)
    {
      // Make the bound variables of the second ultimate delay different from those in the first.
      variable_list renameable_variables=delay2.variables();
      maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
      alphaconvert(renameable_variables, sigma, delay1.variables(), data_expression_list());
      // Additionally map the time variable of the second ultimate delay to that of the first.
      sigma[delay2.time_var()]=delay1.time_var();
      data_expression new_constraint;
      optimized_and(new_constraint, delay1.constraint(), replace_variables_capture_avoiding_alt(delay2.constraint(),sigma));
      variable_list new_existential_variables = delay1.variables()+renameable_variables;

      // TODO: The new constraint can be simplified, as two conditions sharing the timed variable have been merged.
      return lps::detail::ultimate_delay(
                            delay1.time_var(),
                            new_existential_variables,
                            new_constraint);
    }

    void calculate_left_merge_action(
      const lps::detail::ultimate_delay& ultimate_delay_condition,
      const stochastic_action_summand_vector& action_summands1,
      const action_name_multiset_list& allowlist,  // This is a list of list of identifierstring.
      const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
      const bool is_block,
      stochastic_action_summand_vector& action_summands)
    {
      for (const stochastic_action_summand& summand1: action_summands1)
      {
        variable_list sumvars=ultimate_delay_condition.variables();
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(sumvars,sigma,summand1.summation_variables(),data_expression_list());

        variable_list sumvars1 = summand1.summation_variables() + sumvars;
        const action_list& multiaction1 = summand1.multi_action().actions();
        assert(std::is_sorted(multiaction1.begin(), multiaction1.end(), action_compare()));
        data_expression actiontime1 = summand1.multi_action().time();
        data_expression condition1 = summand1.condition();
        const assignment_list& nextstate1 = summand1.assignments();
        const stochastic_distribution& distribution1=summand1.distribution();
        bool has_time = summand1.has_time();

        if (multiaction1 != action_list({ terminationAction }))
        {
          if (is_allow && !allow_(allowlist,multiaction1,terminationAction))
          {
            continue;
          }
          if (is_block && encap(allowlist,multiaction1))
          {
            continue;
          }

          if (!options.ignore_time)
          {
            if (!has_time)
            {
              if (ultimate_delay_condition.constraint()!=sort_bool::true_())
              {
                actiontime1=ultimate_delay_condition.time_var();
                sumvars1.push_front(ultimate_delay_condition.time_var());
                condition1=lazy::and_(condition1,
                                      replace_variables_capture_avoiding_alt(ultimate_delay_condition.constraint(), sigma));
                has_time=true;
              }
            }
            else
            {
              /* Summand1 has time. Substitute the time expression for
                 timevar in ultimate_delay_condition, and extend the condition */
              const std::set<variable> variables_in_actiontime1=find_free_variables(actiontime1);
              sigma[ultimate_delay_condition.time_var()]=actiontime1;
              const data_expression intermediateultimatedelaycondition=
                         replace_variables_capture_avoiding_alt(ultimate_delay_condition.constraint(),sigma);
              optimized_and(condition1, condition1, intermediateultimatedelaycondition);
            }

            condition1=RewriteTerm(condition1);
          }

          if (condition1!=sort_bool::false_())
          {
            action_summands.emplace_back(sumvars1,
                condition1,
                has_time ? multi_action(multiaction1, actiontime1) : multi_action(multiaction1),
                nextstate1,
                distribution1);
          }
        }
      }
    }

    void calculate_left_merge_deadlock(
      const lps::detail::ultimate_delay& ultimate_delay_condition,
      const deadlock_summand_vector& deadlock_summands1,
      const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
      const bool is_block,
      const stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands)
    {
      bool inline_allow = is_allow || is_block;

      if (!inline_allow)
      {
        for (const deadlock_summand& summand1: deadlock_summands1)
        {
          variable_list sumvars=ultimate_delay_condition.variables();
          maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
          alphaconvert(sumvars,sigma,summand1.summation_variables(),data_expression_list());

          variable_list sumvars1=summand1.summation_variables() + sumvars;
          data_expression actiontime1=summand1.deadlock().time();
          data_expression condition1=summand1.condition();
          bool has_time=summand1.deadlock().has_time();

          if (!options.ignore_time)
          {
            if (!has_time)
            {
              if (ultimate_delay_condition.constraint()!=sort_bool::true_())
              {
                actiontime1=ultimate_delay_condition.time_var();
                sumvars1.push_front(ultimate_delay_condition.time_var());
                optimized_and(condition1, condition1,
                                         replace_variables_capture_avoiding_alt(ultimate_delay_condition.constraint(),
                                                                                  sigma));
                has_time=true;
              }
            }
            else
            {
              /* Summand1 has time. Substitute the time expression for
                 timevar in ultimate_delay_condition, and extend the condition */
              const std::set<variable> variables_in_actiontime1=find_free_variables(actiontime1);
              sigma[ultimate_delay_condition.time_var()]=actiontime1;
              const data_expression intermediateultimatedelaycondition=
                         replace_variables_capture_avoiding_alt(ultimate_delay_condition.constraint(),sigma);
              optimized_and(condition1, condition1, intermediateultimatedelaycondition);
            }

            condition1=RewriteTerm(condition1);
          }

          if (condition1!=sort_bool::false_() && !options.ignore_time)
          {
            insert_timed_delta_summand(action_summands,
                                       deadlock_summands,
                                       deadlock_summand(sumvars1,condition1, has_time?deadlock(actiontime1):deadlock()),
                                       options.ignore_time);
          }
        }
      }
    }

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
    static std::string log_leftmerge_application(const lps_statistics_t& lps_statistics_before,
        const lps_statistics_t& lps_statistics_result_before,
        const lps_statistics_t& lps_statistics_result_after,
        const bool is_allow,
        const bool is_block,
        const std::size_t num_allow_block_actions,
        size_t indent = 0)
    {
      std::string indent_str(indent, ' ');
      std::ostringstream os;

      os << indent_str << "- operator: leftmerge" << std::endl;

      indent += 2;
      indent_str = std::string(indent, ' ');
      os << indent_str << "is_allow: " << std::boolalpha << is_allow << std::endl;
      os << indent_str << "is_block: " << std::boolalpha << is_block << std::endl;
      if (is_allow)
      {
        os << indent_str << "number of allow expressions: " << num_allow_block_actions << std::endl;
      }
      if (is_block)
      {
        os << indent_str << "number of block expressions: " << num_allow_block_actions << std::endl;
      }

      os << indent_str << "LPS before:" << std::endl
         << print(lps_statistics_before, indent + 2)
         << indent_str << "result before:" << std::endl
         << print(lps_statistics_result_before, indent + 2)
         << indent_str << "result after:" << std::endl
         << print(lps_statistics_result_after, indent + 2);

      return os.str();
    }
#endif

    void calculate_left_merge(
      const stochastic_action_summand_vector& action_summands1,
      const deadlock_summand_vector& deadlock_summands1,
      const lps::detail::ultimate_delay& ultimate_delay_condition2,
      const action_name_multiset_list& allowlist,  // This is a list of list of identifierstring.
      const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
      const bool is_block,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands)
    {
#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_before = get_statistics(action_summands1, deadlock_summands1);
      lps_statistics_t lps_statistics_result_before = get_statistics(action_summands, deadlock_summands);
#endif

      calculate_left_merge_deadlock(ultimate_delay_condition2, deadlock_summands1,
                                    is_allow, is_block, action_summands, deadlock_summands);
      calculate_left_merge_action(ultimate_delay_condition2, action_summands1,
                                    allowlist, is_allow, is_block, action_summands);

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_result_after = get_statistics(action_summands, deadlock_summands);

      std::cout << log_leftmerge_application(lps_statistics_before,
          lps_statistics_result_before,
          lps_statistics_result_after,
          is_allow,
          is_block,
          (is_block ? allowlist.front().size() : allowlist.size()), 4);
#endif
    }


    void calculate_communication_merge_action_summands(
          const stochastic_action_summand_vector& action_summands1,
          const stochastic_action_summand_vector& action_summands2,
          const action_name_multiset_list& allowlist,   // This is a list of list of identifierstring.
          const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
          const bool is_block,
          stochastic_action_summand_vector& action_summands)
    {
      // First combine the action summands.
      for (const stochastic_action_summand& summand1: action_summands1)
      {
        const variable_list& sumvars1=summand1.summation_variables();
        const action_list multiaction1=summand1.multi_action().actions();
        const data_expression& actiontime1=summand1.multi_action().time();
        const data_expression& condition1=summand1.condition();
        const assignment_list& nextstate1=summand1.assignments();
        const stochastic_distribution& distribution1=summand1.distribution();

        for (const stochastic_action_summand& summand2: action_summands2)
        {
          const variable_list& sumvars2=summand2.summation_variables();
          const action_list multiaction2=summand2.multi_action().actions();
          const data_expression& actiontime2=summand2.multi_action().time();
          const data_expression& condition2=summand2.condition();
          const assignment_list& nextstate2=summand2.assignments();
          const stochastic_distribution& distribution2=summand2.distribution();

          if ((multiaction1 == action_list({ terminationAction })) == (multiaction2 == action_list({ terminationAction })))
          {
            action_list multiaction3;
            if ((multiaction1 == action_list({ terminationAction })) && (multiaction2 == action_list({ terminationAction })))
            {
              multiaction3.push_front(terminationAction);
            }
            else
            {
              multiaction3=linMergeMultiActionList(multiaction1,multiaction2);
            }

            if (is_allow && !allow_(allowlist,multiaction3,terminationAction))
            {
              continue;
            }
            if (is_block && encap(allowlist,multiaction3))
            {
              continue;
            }

            const variable_list allsums=sumvars1+sumvars2;
            data_expression condition3= lazy::and_(condition1,condition2);
            data_expression action_time3;
            bool has_time3=summand1.has_time()||summand2.has_time();

            if (!summand1.has_time())
            {
              if (summand2.has_time())
              {
                /* summand 2 has time*/
                action_time3=actiontime2;
              }
            }
            else
            {
              /* summand 1 has time */
              if (!summand2.has_time())
              {
                action_time3=actiontime1;
              }
              else
              {
                /* both summand 1 and 2 have time */
                action_time3=actiontime1;
                condition3=lazy::and_(
                             condition3,
                             equal_to(actiontime1,actiontime2));
              }
            }

            const assignment_list nextstate3=nextstate1+nextstate2;
            const stochastic_distribution distribution3(
                                              distribution1.variables()+distribution2.variables(),
                                              real_times_optimized(distribution1.distribution(),distribution2.distribution()));

            condition3=RewriteTerm(condition3);
            if (condition3!=sort_bool::false_())
            {
              assert(std::is_sorted(multiaction3.begin(), multiaction3.end(), action_compare()));
              action_summands.emplace_back(allsums,
                  condition3,
                  has_time3 ? multi_action(multiaction3, action_time3) : multi_action(multiaction3),
                  nextstate3,
                  distribution3);
            }
          }
        }
      }
    }

    void calculate_communication_merge_action_deadlock_summands(
          const stochastic_action_summand_vector& action_summands1,
          const deadlock_summand_vector& deadlock_summands1,
          const stochastic_action_summand_vector& action_summands,
          deadlock_summand_vector& deadlock_summands)
    {
      for (const stochastic_action_summand& summand1: action_summands1)
      {
        const variable_list& sumvars1=summand1.summation_variables();
        const data_expression& actiontime1=summand1.multi_action().time();
        const data_expression& condition1=summand1.condition();

        for (const deadlock_summand& summand2: deadlock_summands1)
        {
          const variable_list& sumvars2=summand2.summation_variables();
          const data_expression& actiontime2=summand2.deadlock().time();
          const data_expression& condition2=summand2.condition();

          const variable_list allsums=sumvars1+sumvars2;
          const data_expression condition3= lazy::and_(condition1,condition2);
          data_expression action_time3;
          bool has_time3=summand1.has_time()||summand2.has_time();

          if (!summand1.has_time())
          {
            if (summand2.has_time())
            {
              /* summand 2 has time*/
              action_time3=actiontime2;
            }
          }
          else
          {
            /* summand 1 has time */
            if (!summand2.has_time())
            {
              action_time3=actiontime1;
            }
            else
            {
              /* both summand 1 and 2 have time */
              action_time3=RewriteTerm(sort_real::minimum(actiontime1,actiontime2));
            }
          }

          if (condition3!=sort_bool::false_() && !options.ignore_time)
          {
            insert_timed_delta_summand(action_summands,
                                       deadlock_summands,
                                       deadlock_summand(allsums,
                                                        condition3,
                                                        has_time3?deadlock(action_time3):deadlock()),
                                       options.ignore_time);
          }

        }
      }
    }

    void calculate_communication_merge_deadlock_summands(
          const deadlock_summand_vector& deadlock_summands1,
          const deadlock_summand_vector& deadlock_summands2,
          const stochastic_action_summand_vector& action_summands,
          deadlock_summand_vector& deadlock_summands)
    {
      for (const deadlock_summand& summand1: deadlock_summands1)
      {
        const variable_list& sumvars1=summand1.summation_variables();
        const data_expression& actiontime1=summand1.deadlock().time();
        const data_expression& condition1=summand1.condition();

        for (const deadlock_summand& summand2: deadlock_summands2)
        {
          const variable_list& sumvars2=summand2.summation_variables();
          const data_expression& actiontime2=summand2.deadlock().time();
          const data_expression& condition2=summand2.condition();

          const variable_list allsums=sumvars1+sumvars2;
          const data_expression condition3= lazy::and_(condition1,condition2);
          data_expression action_time3;
          bool has_time3=summand1.has_time()||summand2.has_time();

          if (!summand1.has_time())
          {
            if (summand2.has_time())
            {
              /* summand 2 has time*/
              action_time3=actiontime2;
            }
          }
          else
          {
            /* summand 1 has time */
            if (!summand2.has_time())
            {
              action_time3=actiontime1;
            }
            else
            {
              /* both summand 1 and 2 have time */
              action_time3=RewriteTerm(sort_real::minimum(actiontime1,actiontime2));
            }
          }

          if (condition3!=sort_bool::false_() && !options.ignore_time)
          {
            insert_timed_delta_summand(action_summands,
                                       deadlock_summands,
                                       deadlock_summand(allsums,
                                                        condition3,
                                                        has_time3?deadlock(action_time3):deadlock()),
                                       options.ignore_time);
          }

        }
      }
    }

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
    static std::string log_communicationmerge_application(const lps_statistics_t& lps1_statistics_before,
        const lps_statistics_t& lps2_statistics_before,
        const lps_statistics_t& lps_result_statistics_before,
        const lps_statistics_t& lps_statistics_after,
        const bool is_allow,
        const bool is_block,
        const std::size_t num_allow_block_actions,
        size_t indent = 0)
    {
      std::string indent_str(indent, ' ');
      std::ostringstream os;

      os << indent_str << "- operator: communicationmerge" << std::endl;

      indent += 2;
      indent_str = std::string(indent, ' ');
      os << indent_str << "is_allow: " << std::boolalpha << is_allow << std::endl;
      os << indent_str << "is_block: " << std::boolalpha << is_block << std::endl;
      if (is_allow)
      {
        os << indent_str << "number of allow expressions: " << num_allow_block_actions << std::endl;
      }
      if (is_block)
      {
        os << indent_str << "number of block expressions: " << num_allow_block_actions << std::endl;
      }

      os << indent_str << "first LPS before:" << std::endl << print(lps1_statistics_before, indent + 2)
         << indent_str << "second LPS before:" << std::endl << print(lps2_statistics_before, indent + 2)
         << indent_str << "result LPS before:" << std::endl << print(lps_result_statistics_before, indent + 2)
         << indent_str << "result LPS after:" << std::endl << print(lps_statistics_after, indent + 2);

      return os.str();
    }
#endif

    void calculate_communication_merge(
          const stochastic_action_summand_vector& action_summands1,
          const deadlock_summand_vector& deadlock_summands1,
          const stochastic_action_summand_vector& action_summands2,
          const deadlock_summand_vector& deadlock_summands2,
          const action_name_multiset_list& allowlist,   // This is a list of list of identifierstring.
          const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
          const bool is_block,
          stochastic_action_summand_vector& action_summands,
          deadlock_summand_vector& deadlock_summands)
    {
#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps1_statistics_before = get_statistics(action_summands1, deadlock_summands1);
      lps_statistics_t lps2_statistics_before = get_statistics(action_summands2, deadlock_summands2);
      lps_statistics_t lps_result_statistics_before = get_statistics(action_summands, deadlock_summands);
#endif

      calculate_communication_merge_action_summands(action_summands1, action_summands2, allowlist, is_allow, is_block, action_summands);
      calculate_communication_merge_action_deadlock_summands(action_summands1, deadlock_summands2, action_summands, deadlock_summands);
      calculate_communication_merge_action_deadlock_summands(action_summands2, deadlock_summands1, action_summands, deadlock_summands);
      calculate_communication_merge_deadlock_summands(deadlock_summands1, deadlock_summands2, action_summands, deadlock_summands);

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_after = get_statistics(action_summands, deadlock_summands);

      std::cout << log_communicationmerge_application(lps1_statistics_before, lps2_statistics_before, lps_result_statistics_before,
          lps_statistics_after, is_allow, is_block, (is_block ? allowlist.front().size() : allowlist.size()), 4);
#endif
    }


    void combine_summand_lists(
      const stochastic_action_summand_vector& action_summands1,
      const deadlock_summand_vector& deadlock_summands1,
      const lps::detail::ultimate_delay& ultimate_delay_condition1,
      const stochastic_action_summand_vector& action_summands2,
      const deadlock_summand_vector& deadlock_summands2,
      const lps::detail::ultimate_delay& ultimate_delay_condition2,
      const variable_list& par1,
      const variable_list& par3,
      const action_name_multiset_list& allowlist1,  // This is a list of list of identifierstring.
      const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
      const bool is_block,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands)
    {
      assert(action_summands.size()==0);
      assert(deadlock_summands.size()==0);

      variable_list allpars = par1 + par3;

      if (is_allow || is_block)
      {
        // Inline allow is only supported for ignore_time,
        // for in other cases generation of delta summands cannot be inlined in any simple way.
        assert(!options.nodeltaelimination && options.ignore_time);
        deadlock_summands.emplace_back(variable_list(), sort_bool::true_(), deadlock());
      }

      /* first we enumerate the summands of t1 */

      action_name_multiset_list allowlist((is_allow)?sort_multi_action_labels(allowlist1):allowlist1);
      calculate_left_merge(action_summands1, deadlock_summands1,
                           ultimate_delay_condition2, allowlist, is_allow, is_block,
                           action_summands, deadlock_summands);

      /* second we enumerate the summands of sumlist2 */
      calculate_left_merge(action_summands2, deadlock_summands2,
                           ultimate_delay_condition1, allowlist, is_allow, is_block,
                           action_summands, deadlock_summands);

      /* thirdly we enumerate all multi actions*/

      calculate_communication_merge(action_summands1, deadlock_summands1, action_summands2, deadlock_summands2,
                                    allowlist, is_allow, is_block, action_summands, deadlock_summands);
    }


#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
    static
    std::string log_parallelcomposition_application_start(
                                    const lps_statistics_t& lps1_statistics_before,
                                    const lps_statistics_t& lps2_statistics_before,
                                    const bool is_allow,
                                    const bool is_block,
                                    const std::size_t num_allow_block_actions,
                                    size_t indent = 0)
    {
      std::string indent_str(indent, ' ');
      std::ostringstream os;

      os << indent_str << "- operator: parallel" << std::endl;

      indent += 2;
      indent_str = std::string(indent, ' ');
      os << indent_str << "is_allow: " << std::boolalpha << is_allow << std::endl;
      os << indent_str << "is_block: " << std::boolalpha << is_block << std::endl;
      if (is_allow)
      {
        os << indent_str << "number of allow expressions: " << num_allow_block_actions << std::endl;
      }
      if (is_block)
      {
        os << indent_str << "number of block expressions: " << num_allow_block_actions << std::endl;
      }

      os << indent_str << "first LPS before:" << std::endl << print(lps1_statistics_before, indent+2)
         << indent_str << "second LPS before:" << std::endl << print(lps2_statistics_before, indent+2);

      os << indent_str << "subcalls:" << std::endl;

      return os.str();
    }

  static
  std::string log_parallelcomposition_application_end(
                                      const lps_statistics_t& lps_statistics_after,
                                      size_t indent = 0)
    {
      indent += 2;
      std::string indent_str(indent, ' ');
      std::ostringstream os;

      os << indent_str << "after:" << std::endl << print(lps_statistics_after, indent+2);

      return os.str();
    }
#endif

    void parallelcomposition(
      const stochastic_action_summand_vector& action_summands1,
      const deadlock_summand_vector& deadlock_summands1,
      const variable_list& pars1,
      const data_expression_list& init1,
      const stochastic_distribution& initial_stochastic_distribution1,
      const lps::detail::ultimate_delay& ultimate_delay_condition1,
      const stochastic_action_summand_vector& action_summands2,
      const deadlock_summand_vector& deadlock_summands2,
      const variable_list& pars2,
      const data_expression_list& init2,
      const stochastic_distribution& initial_stochastic_distribution2,
      const lps::detail::ultimate_delay& ultimate_delay_condition2,
      const action_name_multiset_list& allowlist1,  // This is a list of list of identifierstring.
      const bool is_allow,                          // If is_allow or is_block is set, perform inline allow/block filtering.
      const bool is_block,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      variable_list& pars_result,
      data_expression_list& init_result,
      stochastic_distribution& initial_stochastic_distribution,
      lps::detail::ultimate_delay& ultimate_delay_condition)
    {
      mCRL2log(mcrl2::log::verbose) <<
            (is_allow ? "- calculating the parallel composition modulo the allow operator: " :
             is_block ? "- calculating the parallel composition modulo the block operator: " :
                        "- calculating the parallel composition: ") <<
            action_summands1.size() <<
            " action summands + " << deadlock_summands1.size() <<
            " deadlock summands || " << action_summands2.size() <<
            " action summands + " << deadlock_summands2.size() << " deadlock summands = ";

      // At this point the parameters of pars1 and pars2 are unique, except for
      // those that are constant in both processes.

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps1_statistics_before = get_statistics(action_summands1, deadlock_summands1);
      lps_statistics_t lps2_statistics_before = get_statistics(action_summands2, deadlock_summands2);

      std::cout << log_parallelcomposition_application_start(lps1_statistics_before, lps2_statistics_before, is_allow, is_block, (is_block?allowlist1.front().size():allowlist1.size()));
#endif

      variable_list pars3;
      for (const variable& v: pars2)
      {
        if (std::find(pars1.begin(),pars1.end(),v)==pars1.end())
        {
          // *i does not occur in pars1.
          pars3.push_front(v);
        }
      }

      pars3=reverse(pars3);
      assert(action_summands.size()==0);
      assert(deadlock_summands.size()==0);
      combine_summand_lists(action_summands1,deadlock_summands1,ultimate_delay_condition1,
                            action_summands2,deadlock_summands2,ultimate_delay_condition2,
                            pars1,pars3,allowlist1,is_allow,is_block,action_summands,deadlock_summands);

      mCRL2log(mcrl2::log::verbose) << action_summands.size() << " actions and " << deadlock_summands.size() << " delta summands.\n";
      pars_result=pars1+pars3;
      init_result=init1 + init2;
      initial_stochastic_distribution=stochastic_distribution(
                                          initial_stochastic_distribution1.variables()+initial_stochastic_distribution2.variables(),
                                          sort_real::times(initial_stochastic_distribution1.distribution(),
                                                                 initial_stochastic_distribution2.distribution()));
      if (!options.ignore_time)
      {
        ultimate_delay_condition=combine_ultimate_delays(ultimate_delay_condition1, ultimate_delay_condition2);
      }

#ifdef MCRL2_LOG_LPS_LINEARISE_STATISTICS
      lps_statistics_t lps_statistics_after = get_statistics(action_summands, deadlock_summands);

      std::cout << log_parallelcomposition_application_end(lps_statistics_after);
#endif
    }

    /**************** GENERaTE LPEmCRL **********************************/


    /// \brief Linearise a process indicated by procIdDecl.
    /// \details Returns actions_summands, deadlock_summands, the process parameters
    ///              and the initial assignment list.

    void generateLPEmCRLterm(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const process_expression& t,
      const bool regular,
      const bool rename_variables,
      variable_list& pars,
      data_expression_list& init,
      stochastic_distribution& initial_stochastic_distribution,
      lps::detail::ultimate_delay& ultimate_delay_condition)
    {
      if (is_process_instance_assignment(t))
      {
        generateLPEmCRL(action_summands,deadlock_summands,process_instance_assignment(t).identifier(),
                        regular,pars,init,initial_stochastic_distribution,ultimate_delay_condition);

        objectdatatype& object=objectIndex(process_instance_assignment(t).identifier());

        // Now apply the assignment in this process to the obtained initial process and the distribution.
        maintain_variables_in_rhs<mutable_map_substitution<> > sigma;
        for (const assignment& a: process_instance_assignment(t).assignments())
        {
          sigma[a.lhs()]=a.rhs();
        }

        init=replace_variables_capture_avoiding_alt(init,sigma);
        initial_stochastic_distribution =
                   stochastic_distribution(
                         initial_stochastic_distribution.variables(),
                         replace_variables_capture_avoiding_alt(initial_stochastic_distribution.distribution(), sigma));

        // Make the bound variables and parameters in this process unique.

        if ((object.processstatus==GNF)||
            (object.processstatus==pCRL)||
            (object.processstatus==GNFalpha))
        {
          make_parameters_and_sum_variables_unique(action_summands,deadlock_summands,pars,ultimate_delay_condition,std::string(object.objectname));
        }
        else
        {
          if (rename_variables)
          {
            make_parameters_and_sum_variables_unique(action_summands,deadlock_summands,pars,ultimate_delay_condition);
          }
        }

        if (regular && !options.do_not_apply_constelm)
        {
          // We apply constant elimination on the obtained linear process.
          // In order to do so, we have to create a complete process specification first, as
          // this is what the interface of constelm requires.
          // Note that this is only useful, in regular mode. This does not make sense if
          // stacks are being used.

          stochastic_linear_process lps(pars,deadlock_summands,action_summands);
          stochastic_process_initializer initializer(init,initial_stochastic_distribution);

          stochastic_specification temporary_spec(data,acts,global_variables,lps,initializer);
          constelm_algorithm < rewriter, stochastic_specification > alg(temporary_spec,rewr);

          // Remove constants from the specification, where global variables are
          // also instantiated if they exist.
          data::mutable_map_substitution<> sigma = alg.compute_constant_parameters(true); // Instantiate global variables; Take conditions into account.
          alg.remove_parameters(sigma);
          alg.remove_trivial_summands();

          if (!options.ignore_time)
          {
            ultimate_delay_condition.constraint()=data::replace_variables(ultimate_delay_condition.constraint(),sigma);
          }

          // Reconstruct the variables from the temporary specification
          init=temporary_spec.initial_process().expressions();
          pars=temporary_spec.process().process_parameters();
          assert(init.size()==pars.size());
          // Add all free variables in object.parameters that are not already in the parameter list
          // and are not global variables to pars. This can occur when a parameter of the process is replaced
          // by a constant, which by itself is a parameter.

          std::set <variable> variable_list = lps::find_free_variables(temporary_spec.process().action_summands());
          const std::set <variable> variable_list1 = lps::find_free_variables(temporary_spec.process().deadlock_summands());
          variable_list.insert(variable_list1.begin(),variable_list1.end());
          for (const variable& v: variable_list)
          {
            if (std::find(pars.begin(),pars.end(),v)==pars.end() && // The free variable is not in pars,
                global_variables.find(v)==global_variables.end() // it is neither a global variable
                // (lps::search_free_variable(temporary_spec.process().action_summands(),v) || lps::search_free_variable(temporary_spec.process().deadlock_summands(),v))
               )          // and it occurs in the summands.
            {
              pars.push_front(v);
            }
          }

          action_summands=temporary_spec.process().action_summands();
          deadlock_summands=temporary_spec.process().deadlock_summands();
        }
        // Now constelm has been applied.
        return;
      } // End process assignment.

      if (is_merge(t))
      {
        variable_list pars1;
        variable_list pars2;
        data_expression_list init1;
        data_expression_list init2;
        stochastic_distribution initial_stochastic_distribution1;
        stochastic_distribution initial_stochastic_distribution2;
        stochastic_action_summand_vector action_summands1;
        stochastic_action_summand_vector action_summands2;
        deadlock_summand_vector deadlock_summands1;
        deadlock_summand_vector deadlock_summands2;
        lps::detail::ultimate_delay ultimate_delay_condition1;
        lps::detail::ultimate_delay ultimate_delay_condition2;
        generateLPEmCRLterm(action_summands1,deadlock_summands1,process::merge(t).left(),
                              regular,rename_variables,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1);
        generateLPEmCRLterm(action_summands2,deadlock_summands2,process::merge(t).right(),
                              regular,true,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2);
        parallelcomposition(action_summands1,deadlock_summands1,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1,
                              action_summands2,deadlock_summands2,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2,
                              action_name_multiset_list(),false,false,
                              action_summands,deadlock_summands,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        return;
      }

      if (is_hide(t))
      {
        generateLPEmCRLterm(action_summands,deadlock_summands,hide(t).operand(),
                              regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        hidecomposition(hide(t).hide_set(),action_summands);
        return;
      }

      if (is_allow(t))
      {
        process_expression par = allow(t).operand();
        if (!options.nodeltaelimination && options.ignore_time && is_merge(par))
        {
          // Perform parallel composition with inline allow.
          variable_list pars1;
          variable_list pars2;
          data_expression_list init1;
          data_expression_list init2;
          stochastic_distribution initial_stochastic_distribution1;
          stochastic_distribution initial_stochastic_distribution2;
          stochastic_action_summand_vector action_summands1;
          stochastic_action_summand_vector action_summands2;
          deadlock_summand_vector deadlock_summands1;
          deadlock_summand_vector deadlock_summands2;
          lps::detail::ultimate_delay ultimate_delay_condition1;
          lps::detail::ultimate_delay ultimate_delay_condition2;
          generateLPEmCRLterm(action_summands1,deadlock_summands1,process::merge(par).left(),
                                regular,rename_variables,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1);
          generateLPEmCRLterm(action_summands2,deadlock_summands2,process::merge(par).right(),
                                regular,true,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2);
          parallelcomposition(action_summands1,deadlock_summands1,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1,
                                action_summands2,deadlock_summands2,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2,
                                allow(t).allow_set(),true,false,
                                action_summands,deadlock_summands,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
          return;
        }
        else if (!options.nodeltaelimination && options.ignore_time && is_comm(par))
        {
          generateLPEmCRLterm(action_summands,deadlock_summands,comm(par).operand(),
                                regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
          communicationcomposition(comm(par).comm_set(),allow(t).allow_set(),true,false,action_summands,deadlock_summands,terminationAction,options.nosumelm,options.nodeltaelimination,options.ignore_time,[this](const data::data_expression& e) { return RewriteTerm(e); });
          return;
        }

        generateLPEmCRLterm(action_summands,deadlock_summands,par,regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        allowblockcomposition(allow(t).allow_set(),true,action_summands,deadlock_summands,terminationAction,options.ignore_time,options.nodeltaelimination);
        return;
      }

      if (is_block(t))
      {
        process_expression par = block(t).operand();
        if (!options.nodeltaelimination && options.ignore_time && is_merge(par))
        {
          // Perform parallel composition with inline block.
          variable_list pars1;
          variable_list pars2;
          data_expression_list init1;
          data_expression_list init2;
          stochastic_distribution initial_stochastic_distribution1;
          stochastic_distribution initial_stochastic_distribution2;
          stochastic_action_summand_vector action_summands1;
          stochastic_action_summand_vector action_summands2;
          deadlock_summand_vector deadlock_summands1;
          deadlock_summand_vector deadlock_summands2;
          lps::detail::ultimate_delay ultimate_delay_condition1;
          lps::detail::ultimate_delay ultimate_delay_condition2;
          generateLPEmCRLterm(action_summands1,deadlock_summands1,process::merge(par).left(),
                                regular,rename_variables,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1);
          generateLPEmCRLterm(action_summands2,deadlock_summands2,process::merge(par).right(),
                                regular,true,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2);
          // Encode the actions of the block list in one multi action.
          parallelcomposition(action_summands1,deadlock_summands1,pars1,init1,initial_stochastic_distribution1,ultimate_delay_condition1,
                                action_summands2,deadlock_summands2,pars2,init2,initial_stochastic_distribution2,ultimate_delay_condition2,
                                action_name_multiset_list({action_name_multiset(block(t).block_set())}),false,true,
                                action_summands,deadlock_summands,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
          return;
        }
        else if (!options.nodeltaelimination && options.ignore_time && is_comm(par))
        {
          generateLPEmCRLterm(action_summands,deadlock_summands,comm(par).operand(),
                                regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
          // Encode the actions of the block list in one multi action.
          communicationcomposition(comm(par).comm_set(),action_name_multiset_list( { action_name_multiset(block(t).block_set())} ),
                                                     false,true,action_summands,deadlock_summands,terminationAction,options.nosumelm,options.nodeltaelimination,options.ignore_time,[this](const data::data_expression& e) { return RewriteTerm(e); });
          return;
        }

        generateLPEmCRLterm(action_summands,deadlock_summands,par,regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        // Encode the actions of the block list in one multi action.
        allowblockcomposition(action_name_multiset_list({action_name_multiset(block(t).block_set())}),false,action_summands,deadlock_summands,terminationAction,options.ignore_time,options.nodeltaelimination);
        return;
      }

      if (is_rename(t))
      {
        generateLPEmCRLterm(action_summands,deadlock_summands,process::rename(t).operand(),
                              regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        lps::rename(process::rename(t).rename_set(),action_summands);
        return;
      }

      if (is_comm(t))
      {
        generateLPEmCRLterm(action_summands,deadlock_summands,comm(t).operand(),
                              regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        communicationcomposition(comm(t).comm_set(),action_name_multiset_list(),false,false,action_summands,deadlock_summands,terminationAction,options.nosumelm,options.nodeltaelimination,options.ignore_time,[this](const data::data_expression& e) { return RewriteTerm(e); });
        return;
      }

      if (is_stochastic_operator(t))
      {
        /* YYYYYYYYYY */
        const stochastic_operator& sto=atermpp::down_cast<stochastic_operator>(t);
        generateLPEmCRLterm(action_summands,deadlock_summands,sto.operand(),
                              regular,rename_variables,pars,init,initial_stochastic_distribution,ultimate_delay_condition);
        variable_list stochvars=sto.variables();
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconvert(stochvars,sigma,pars + initial_stochastic_distribution.variables(), data_expression_list());
        initial_stochastic_distribution=stochastic_distribution(
                                          stochvars+initial_stochastic_distribution.variables(),
                                          sort_real::times(replace_variables_capture_avoiding_alt(
                                                                                    sto.distribution(), sigma),
                                                                 initial_stochastic_distribution.distribution()));
        /* Reset the bound variables in the initial_stochastic_distribution, to avoid erroneous renaming in the body of the process */
        for(const variable& v: initial_stochastic_distribution.variables())
        {
          sigma[v]=v;
        }
        init=replace_variables_capture_avoiding_alt(init,sigma);
        return;
      }

      throw mcrl2::runtime_error("Internal error. Expect an mCRL term " + process::pp(t) +".");
    }

    /**************** GENERaTE LPEmCRL **********************************/

    /* The result are a list of action summands, deadlock summand, the parameters of this
       linear process and its initial values. A initial stochastic distribution that must
       precede the initial linear process and the ultimate delay condition of this
       linear process that can be used or be ignored.

    */

    void generateLPEmCRL(
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      const process_identifier& procIdDecl,
      const bool regular,
      variable_list& pars,
      data_expression_list& init,
      stochastic_distribution& initial_stochastic_distribution,
      lps::detail::ultimate_delay& ultimate_delay_condition)
    {
      /* If regular=1, then a regular version of the pCRL processes
         must be generated */

      objectdatatype& object=objectIndex(procIdDecl);

      if ((object.processstatus==GNF)||
          (object.processstatus==pCRL)||
          (object.processstatus==GNFalpha)||
          (object.processstatus==multiAction))
      {
        generateLPEpCRL(action_summands,deadlock_summands,procIdDecl,
                               object.containstime,regular,pars,init,initial_stochastic_distribution);
        if (options.ignore_time)
        {
          ultimate_delay_condition=lps::detail::ultimate_delay();
        }
        else
        {
          ultimate_delay_condition=getUltimateDelayCondition(action_summands,deadlock_summands,pars);

          try
          {
            // The ultimate_delay_condition can be complex. Try to simplify it with a fourier_motzkin reduction.
            data_expression simplified_ultimate_delay_condition;
            variable_list reduced_sumvars;

            fourier_motzkin(ultimate_delay_condition.constraint(),
                            ultimate_delay_condition.variables(),
                            simplified_ultimate_delay_condition,
                            reduced_sumvars,
                            rewr);
            std::swap(ultimate_delay_condition.constraint(), simplified_ultimate_delay_condition);
            std::swap(ultimate_delay_condition.variables(), reduced_sumvars);
          }
          catch (mcrl2::runtime_error& e)
          {
            // Applying Fourier Motzkin failed. Continue working with the old ultimate delay condition.
            mCRL2log(mcrl2::log::debug) << "Simplifying a condition using Fourier-Motzkin reduction failed (I). \n" << e.what() << std::endl;
          }
        }

        assert(init.size() == pars.size());
        return;
      }
      /* process is a mCRLdone */
      if ((object.processstatus==mCRLdone)||
          (object.processstatus==mCRLlin)||
          (object.processstatus==mCRL))
      {
        object.processstatus=mCRLlin;
        generateLPEmCRLterm(action_summands, deadlock_summands, object.processbody,
                                   regular, false, pars, init, initial_stochastic_distribution, ultimate_delay_condition);
        assert(init.size() == pars.size());
        return;
      }

      throw mcrl2::runtime_error("laststatus: " + std::to_string(object.processstatus));
    }

    /**************** alphaconversion ********************************/

    process_expression alphaconversionterm(
      const process_expression& t,
      const variable_list& parameters,
      maintain_variables_in_rhs<mutable_map_substitution<> > sigma)
    {
      if (is_choice(t))
      {
        return choice(
                 alphaconversionterm(choice(t).left(),parameters,sigma),
                 alphaconversionterm(choice(t).right(),parameters,sigma));
      }

      if (is_seq(t))
      {
        return seq(
                 alphaconversionterm(seq(t).left(),parameters,sigma),
                 alphaconversionterm(seq(t).right(),parameters,sigma));
      }

      if (is_sync(t))
      {
        return process::sync(
                 alphaconversionterm(process::sync(t).left(),parameters,sigma),
                 alphaconversionterm(process::sync(t).right(),parameters,sigma));
      }

      if (is_bounded_init(t))
      {
        return bounded_init(
                 alphaconversionterm(bounded_init(t).left(),parameters,sigma),
                 alphaconversionterm(bounded_init(t).right(),parameters,sigma));
      }

      if (is_merge(t))
      {
        return process::merge(
                 alphaconversionterm(process::merge(t).left(),parameters,sigma),
                 alphaconversionterm(process::merge(t).right(),parameters,sigma));
      }

      if (is_left_merge(t))
      {
        return left_merge(
                 alphaconversionterm(left_merge(t).left(),parameters,sigma),
                 alphaconversionterm(left_merge(t).right(),parameters,sigma));
      }

      if (is_at(t))
      {
        return at(alphaconversionterm(at(t).operand(),parameters,sigma),
                  replace_variables_capture_avoiding_alt(at(t).time_stamp(),sigma));
      }

      if (is_if_then(t))
      {
        return if_then(
                 replace_variables_capture_avoiding_alt(if_then(t).condition(), sigma),
                 alphaconversionterm(if_then(t).then_case(),parameters,sigma));
      }

      if (is_sum(t))
      {
        variable_list sumvars=sum(t).variables();
        maintain_variables_in_rhs<mutable_map_substitution<> > local_sigma=sigma;

        alphaconvert(sumvars,local_sigma,variable_list(),variable_list_to_data_expression_list(parameters));
        return sum(sumvars,alphaconversionterm(sum(t).operand(), sumvars+parameters, local_sigma));
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(t);
        variable_list sumvars=sto.variables();
        maintain_variables_in_rhs<mutable_map_substitution<> > local_sigma=sigma;

        alphaconvert(sumvars,local_sigma,variable_list(),variable_list_to_data_expression_list(parameters));
        return stochastic_operator(
                       sumvars,
                       replace_variables_capture_avoiding_alt(sto.distribution(), sigma),
                       alphaconversionterm(sto.operand(), sumvars+parameters, local_sigma));
      }

      if (is_process_instance_assignment(t))
      {
        const process_identifier procId=process_instance_assignment(t).identifier();
        objectdatatype& object=objectIndex(procId);

        const variable_list instance_parameters=object.parameters;
        alphaconversion(procId,instance_parameters);

        const process_instance_assignment result(procId,
                                                 substitute_assignmentlist(process_instance_assignment(t).assignments(),
                                                                           instance_parameters,false, true,sigma));
        assert(check_valid_process_instance_assignment(result.identifier(),result.assignments()));
        return result;
      }

      if (is_action(t))
      {
        return action(action(t).label(), replace_variables_capture_avoiding_alt(action(t).arguments(), sigma));
      }

      if (is_delta(t)||
          is_tau(t))
      {
        return t;
      }

      if (is_hide(t))
      {
        return alphaconversionterm(hide(t).operand(),parameters,sigma);
      }

      if (is_rename(t))
      {
        return alphaconversionterm(process::rename(t).operand(),parameters,sigma);
      }

      if (is_comm(t))
      {
        return alphaconversionterm(comm(t).operand(),parameters,sigma);
      }

      if (is_allow(t))
      {
        return alphaconversionterm(allow(t).operand(),parameters,sigma);
      }

      if (is_block(t))
      {
        return alphaconversionterm(block(t).operand(),parameters,sigma);
      }

      throw mcrl2::runtime_error("unexpected process format in alphaconversionterm " + process::pp(t) +".");
      return process_expression();
    }

    void alphaconversion(const process_identifier& procId, const variable_list& parameters)
    {
      objectdatatype& object=objectIndex(procId);

      if ((object.processstatus==GNF)||
          (object.processstatus==multiAction))
      {
        object.processstatus=GNFalpha;
        // tempvar below is needed as objectdata may be reallocated
        // during a call to alphaconversionterm.
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        const process_expression tempvar=alphaconversionterm(object.processbody,parameters,sigma);
        object.processbody=tempvar;
      }
      else if (object.processstatus==mCRLdone)
      {
        maintain_variables_in_rhs< mutable_map_substitution<> > sigma;
        alphaconversionterm(object.processbody,parameters,sigma);

      }
      else if (object.processstatus==GNFalpha)
      {
        return;
      }
      else
      {
        throw mcrl2::runtime_error("unknown type " + std::to_string(object.processstatus) +
                                                " in alphaconversion of " + process::pp(procId) +".");
      }
      return;
    }

    /***** determinewhetherprocessescontaintime; **********/

    bool containstimebody(
      const process_expression& t,
      bool* stable,
      std::set < process_identifier >& visited,
      bool allowrecursion,
      bool& contains_if_then)
    {
      if (is_merge(t))
      {
        /* the construction below is needed to guarantee that
           both subterms are recursively investigated */
        bool r1=containstimebody(process::merge(t).left(),stable,visited,allowrecursion,contains_if_then);
        bool r2=containstimebody(process::merge(t).right(),stable,visited,allowrecursion,contains_if_then);
        return r1||r2;
      }

      if (is_process_instance_assignment(t))
      {
        if (allowrecursion)
        {
          return (containstime_rec(process_instance_assignment(t).identifier(),stable,visited,contains_if_then));
        }
        return objectIndex(process_instance_assignment(t).identifier()).containstime;
      }

      if (is_hide(t))
      {
        return containstimebody(hide(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_rename(t))
      {
        return containstimebody(process::rename(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_allow(t))
      {
        return containstimebody(allow(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_block(t))
      {
        return containstimebody(block(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_comm(t))
      {
        return containstimebody(comm(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_choice(t))
      {
        bool r1=containstimebody(choice(t).left(),stable,visited,allowrecursion,contains_if_then);
        bool r2=containstimebody(choice(t).right(),stable,visited,allowrecursion,contains_if_then);
        return r1||r2;
      }

      if (is_seq(t))
      {
        bool r1=containstimebody(seq(t).left(),stable,visited,allowrecursion,contains_if_then);
        bool r2=containstimebody(seq(t).right(),stable,visited,allowrecursion,contains_if_then);
        return r1||r2;
      }

      if (is_if_then(t))
      {
        // If delta is added, c->p is translated into c->p<>delta,
        // otherwise into c->p<>delta@0. In this last case the process
        // contains time.
        contains_if_then=true;
        if (options.ignore_time)
        {
          return containstimebody(if_then(t).then_case(),stable,visited,allowrecursion,contains_if_then);
        }
        else
        {
          return true;
        }
      }

      if (is_if_then_else(t))
      {
        bool r1=containstimebody(if_then_else(t).then_case(),stable,visited,allowrecursion,contains_if_then);
        bool r2=containstimebody(if_then_else(t).else_case(),stable,visited,allowrecursion,contains_if_then);
        return r1||r2;
      }

      if (is_sum(t))
      {
        return containstimebody(sum(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_stochastic_operator(t))
      {
        return containstimebody(stochastic_operator(t).operand(),stable,visited,allowrecursion,contains_if_then);
      }

      if (is_action(t)||
          is_delta(t)||
          is_tau(t))
      {
        return false;
      }

      if (is_at(t))
      {
        return true;
      }

      if (is_sync(t))
      {
        bool r1=containstimebody(process::sync(t).left(),stable,visited,allowrecursion,contains_if_then);
        bool r2=containstimebody(process::sync(t).right(),stable,visited,allowrecursion,contains_if_then);
        return r1||r2;
      }

      throw mcrl2::runtime_error("unexpected process format in containstime " + process::pp(t) +".");
      return false;
    }

    bool containstime_rec(
      const process_identifier& procId,
      bool* stable,
      std::set < process_identifier >& visited,
      bool& contains_if_then)
    {
      objectdatatype& object=objectIndex(procId);

      if (visited.count(procId)==0)
      {
        visited.insert(procId);
        const bool ct=containstimebody(object.processbody,stable,visited,true,contains_if_then);
        static bool show_only_once=true;
        if (ct && options.ignore_time && show_only_once)
        {
          mCRL2log(mcrl2::log::warning) << "process " << procId.name() <<
              " contains time, which is now not preserved. \n"  <<
              "Use --timed or -T, or untick `add deadlocks' for a correct timed linearisation...\n";
          show_only_once=false;
        }
        if (object.containstime!=ct)
        {
          object.containstime=ct;
          if (stable!=nullptr)
          {
            *stable=false;
          }
        }
      }
      return (object.containstime);
    }

    bool containstimebody(const process_expression& t)
    {
      std::set < process_identifier > visited;
      bool stable;
      bool contains_if_then;
      return containstimebody(t,&stable,visited,false,contains_if_then);
    }

    bool determinewhetherprocessescontaintime(const process_identifier& procId)

    {
      /* This function sets for all reachable processes in the array objectdata
         whether they contain time in the field containstime. In verbose mode
         it prints the process variables that contain time. Furtermore, it returns
         whether there are processes that contain an if-then that will be translated
         to an if-then-else with an delta@0 in the else branch, introducing time */
      bool stable=false;
      bool contains_if_then=false;

      while (!stable)
      {
        std::set < process_identifier > visited;
        stable=true;
        containstime_rec(procId,&stable,visited,contains_if_then);
      }
      return contains_if_then;
    }
    /***** transform_initial_distribution_term **********/

    process_expression transform_initial_distribution_term(
                            const process_expression& t,
                            const std::map < process_identifier, process_pid_pair >& processes_with_initial_distribution)
    {
      if (is_process_instance_assignment(t))
      {
        const process_instance_assignment u=atermpp::down_cast<process_instance_assignment>(t);
        const process_expression new_process=processes_with_initial_distribution.at(u.identifier()).process_body();
        if (is_stochastic_operator(new_process))
        {
          // Add the initial stochastic_distribution.
          const process_identifier& new_identifier=processes_with_initial_distribution.at(u.identifier()).process_id();
          objectdatatype& object=objectIndex(new_identifier);
          const variable_list new_parameters=object.parameters;
          const stochastic_operator& sto=down_cast<stochastic_operator>(new_process);
          assignment_list new_assignments;

          const variable_list relevant_stochastic_variables=parameters_that_occur_in_body(sto.variables(),object.processbody);
          assert(relevant_stochastic_variables.size()<=new_parameters.size());
          // The variables in sto.variables() may clash with local variables u and therefore need to be rename.
          variable_list renamed_sto_variables=sto.variables();
          mutable_indexed_substitution<> local_sigma1;
          alphaconvertprocess(renamed_sto_variables, local_sigma1, t);


          variable_list::const_iterator i=new_parameters.begin();
          for(const variable& v: relevant_stochastic_variables)
          {
            new_assignments=push_back(new_assignments,assignment(*i,local_sigma1(v)));
            ++i;
          }
          data_expression new_distribution = data::replace_variables_capture_avoiding(sto.distribution(),local_sigma1);

          // Some of the variables may occur only in the distribution, which is now moved out.
          // Therefore, the assignments must be filtered.
          new_assignments=filter_assignments(new_assignments + u.assignments(),object.parameters);

          // The variables bound in the distribution may conflict

          // Furthermore, the old assignment must be applied to the distribution, when it is moved
          // outside of the process body.
          mutable_indexed_substitution<> local_sigma;
          for(const assignment& a:u.assignments())
          {
            local_sigma[a.lhs()]=a.rhs();
          }
          return stochastic_operator(renamed_sto_variables,
                                     data::replace_variables_capture_avoiding(new_distribution, local_sigma),
                                     process_instance_assignment(new_identifier,new_assignments));
        }
        return t;

      }

      if (is_choice(t))
      {
        const process_expression r1_=transform_initial_distribution_term(choice(t).left(),processes_with_initial_distribution);
        const process_expression r2_=transform_initial_distribution_term(choice(t).right(),processes_with_initial_distribution);
        if (is_stochastic_operator(r1_))
        {
          const stochastic_operator& r1=down_cast<const stochastic_operator>(r1_);
          if (is_stochastic_operator(r2_))
          {
            /* both r1_ and r2_ are stochastic */
            const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
            const process_expression& new_body1=r1.operand();
            process_expression new_body2=r2.operand();
            variable_list stochvars=r2.variables();
            maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
            alphaconvert(stochvars,local_sigma, r1.variables(),data_expression_list());
            new_body2=substitute_pCRLproc(new_body2, local_sigma);
            const data_expression new_distribution=
                   process::replace_variables_capture_avoiding_with_an_identifier_generator(r2.distribution(),local_sigma,fresh_identifier_generator);

            return stochastic_operator(r1.variables() + stochvars,
                                       real_times_optimized(r1.distribution(),new_distribution),
                                       choice(new_body1,new_body2));
          }
          /* r1 is and r2_ is not a stochastic operator */
          return stochastic_operator(r1.variables(),r1.distribution(),choice(r1.operand(),r2_));
        }
        if (is_stochastic_operator(r2_))
        {
          /* r1_ is not and r2_ is a stochastic operator */
          const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
          return stochastic_operator(r2.variables(),r2.distribution(),choice(r1_,r2.operand()));
        }
        /* neither r1_ nor r2_ is stochastic */
        return choice(r1_,r2_);
      }

      if (is_seq(t))
      {
        const process_expression r1=transform_initial_distribution_term(seq(t).left(),processes_with_initial_distribution);
        const process_expression r2=transform_initial_distribution_term(seq(t).right(),processes_with_initial_distribution);
        if (is_stochastic_operator(r1))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r1);
          return stochastic_operator(r.variables(),r.distribution(),seq(r.operand(),r2));
        }
        return seq(r1,r2);
      }

      if (is_if_then(t))
      {
        const if_then& t_=atermpp::down_cast<if_then>(t);
        const process_expression r=transform_initial_distribution_term(t_.then_case(),processes_with_initial_distribution);
        if (is_stochastic_operator(r))
        {
          const stochastic_operator& r_=down_cast<const stochastic_operator>(r);
          return stochastic_operator(r_.variables(),
                                     if_(t_.condition(),r_.distribution(),if_(variables_are_equal_to_default_values(r_.variables()),sort_real::real_one(),sort_real::real_zero())),
                                     if_then(t_.condition(),r_.operand()));
        }
        return if_then(t_.condition(),r);
      }

      if (is_if_then_else(t))
      {
        const if_then_else& t_=atermpp::down_cast<if_then_else>(t);
        const process_expression r1_=transform_initial_distribution_term(t_.then_case(),processes_with_initial_distribution);
        const process_expression r2_=transform_initial_distribution_term(t_.else_case(),processes_with_initial_distribution);
        if (is_stochastic_operator(r1_))
        {
          const stochastic_operator& r1=down_cast<const stochastic_operator>(r1_);
          if (is_stochastic_operator(r2_))
          {
            /* both r1_ and r2_ are stochastic */

            const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
            const process_expression& new_body1=r1.operand();
            process_expression new_body2=r2.operand();
            variable_list stochvars=r2.variables();
            maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
            alphaconvert(stochvars,local_sigma, r1.variables(),data_expression_list());
            new_body2=substitute_pCRLproc(new_body2, local_sigma);
            const data_expression new_distribution=
                   process::replace_variables_capture_avoiding_with_an_identifier_generator(r2.distribution(),local_sigma,fresh_identifier_generator);

            return stochastic_operator(r1.variables() + stochvars,
                                       if_(t_.condition(),
                                               lazy::and_(variables_are_equal_to_default_values(stochvars),r1.distribution()),
                                               lazy::and_(variables_are_equal_to_default_values(r1.variables()),new_distribution)),
                                       if_then_else(t_.condition(),new_body1,new_body2));
          }
          /* r1 is and r2_ is not a stochastic operator */
          return stochastic_operator(r1.variables(),
                                     if_(t_.condition(),r1.distribution(),
                                               if_(variables_are_equal_to_default_values(r1.variables()),sort_real::real_one(),sort_real::real_zero())),
                                     if_then_else(if_then_else(t).condition(),r1.operand(),r2_));
        }
        if (is_stochastic_operator(r2_))
        {
          /* r1_ is not and r2_ is a stochastic operator */
          const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
          return stochastic_operator(r2.variables(),
                                     if_(t_.condition(),
                                         if_(variables_are_equal_to_default_values(r2.variables()),sort_real::real_one(),sort_real::real_zero()),
                                         r2.distribution()),
                                     if_then_else(t_.condition(),r1_,r2.operand()));
        }
        /* neither r1_ nor r2_ is stochastic */
        return if_then_else(t_.condition(),r1_,r2_);

      }

      if (is_sum(t))
      {
        const sum& s=down_cast<const sum>(t);
        const process_expression r_= transform_initial_distribution_term(s.operand(),processes_with_initial_distribution);
        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          std::set <variable> variables_in_distribution=find_all_variables(r.distribution());
          for(const variable& v: s.variables())
          {
            if (variables_in_distribution.count(v)>0)
            {
              throw mcrl2::runtime_error("Cannot commute a sum operator over a stochastic operator in " +
                                                   process::pp(t) + ".\n" +
                                         "The problematic variable is " + pp(v) + ":" + pp(v.sort()) + ".");
            }
          }
          return stochastic_operator(r.variables(),
                                     r.distribution(),
                                     sum(s.variables(),r.operand()));
        }
        return sum(s.variables(),r_);
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(t);
        const process_expression r_= transform_initial_distribution_term(sto.operand(),processes_with_initial_distribution);

        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(sto.variables()+r.variables(),
                                     real_times_optimized(sto.distribution(),r.distribution()),
                                     r.operand());
        }
        return stochastic_operator(sto.variables(),sto.distribution(),r_);
      }

      if (is_action(t))
      {
        return t;
      }

      if (is_delta(t))
      {
        return t;
      }

      if (is_tau(t))
      {
        return t;
      }

      if (is_at(t))
      {
        const process_expression r_=transform_initial_distribution_term(at(t).operand(),processes_with_initial_distribution);
        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(r.variables(),r.distribution(),at(r.operand(),at(t).time_stamp()));
        }
        return t;
      }

      if (is_sync(t))
      {
        return t;
      }

      throw mcrl2::runtime_error("unexpected process format in transform_initial_distribution_term " + process::pp(t) +".");
    }

    /***** obtain_initial_distribution **********/

    process_expression obtain_initial_distribution_term(const process_expression& t)
    {
      /* This function obtains the initial distribution of a pCRL term that is in Greibach normal form.  */
      if (is_process_instance_assignment(t))
      {
        assert(0); // This is not possible.
        return t;
      }

      if (is_choice(t))
      {
        const process_expression r1_=obtain_initial_distribution_term(choice(t).left());
        const process_expression r2_=obtain_initial_distribution_term(choice(t).right());
        if (is_stochastic_operator(r1_))
        {
          const stochastic_operator& r1=down_cast<const stochastic_operator>(r1_);
          if (is_stochastic_operator(r2_))
          {
            /* both r1_ and r2_ are stochastic */
            const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
            const process_expression& new_body1=r1.operand();
            process_expression new_body2=r2.operand();
            variable_list stochvars=r2.variables();
            maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
            alphaconvert(stochvars,local_sigma, r1.variables(),data_expression_list());
            new_body2=substitute_pCRLproc(new_body2, local_sigma);
            const data_expression new_distribution=
                   process::replace_variables_capture_avoiding_with_an_identifier_generator(r2.distribution(),local_sigma,fresh_identifier_generator);

            return stochastic_operator(r1.variables() + stochvars,
                                       real_times_optimized(r1.distribution(),new_distribution),
                                       choice(new_body1,new_body2));
          }
          /* r1 is and r2_ is not a stochastic operator */
          return stochastic_operator(r1.variables(),r1.distribution(),choice(r1.operand(),r2_));
        }
        if (is_stochastic_operator(r2_))
        {
          /* r1_ is not and r2_ is a stochastic operator */
          const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
          return stochastic_operator(r2.variables(),r2.distribution(),choice(r1_,r2.operand()));
        }
        /* neither r1_ nor r2_ is stochastic */
        return t;
      }

      if (is_seq(t))
      {
        const process_expression r_=obtain_initial_distribution_term(down_cast<seq>(t).left());
        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(r.variables(),r.distribution(),seq(r.operand(),seq(t).right()));
        }
        return t;
      }

      if (is_if_then(t))
      {
        const process_expression r_=obtain_initial_distribution_term(down_cast<if_then>(t).then_case());
        if (is_stochastic_operator(r_))
        {
          const if_then& t_if_then=atermpp::down_cast<if_then>(t);
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(r.variables(),
                                     if_(t_if_then.condition(),r.distribution(),if_(variables_are_equal_to_default_values(r.variables()),sort_real::real_one(),sort_real::real_zero())),
                                     if_then(t_if_then.condition(),r.operand()));
        }
        return t;
      }

      if (is_if_then_else(t))
      {
        const process_expression r1_=obtain_initial_distribution_term(down_cast<if_then_else>(t).then_case());
        const process_expression r2_=obtain_initial_distribution_term(down_cast<if_then_else>(t).else_case());
        if (is_stochastic_operator(r1_))
        {
          const stochastic_operator& r1=down_cast<const stochastic_operator>(r1_);
          if (is_stochastic_operator(r2_))
          {
            /* both r1_ and r2_ are stochastic */

            const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
            const process_expression& new_body1=r1.operand();
            process_expression new_body2=r2.operand();
            variable_list stochvars=r2.variables();
            maintain_variables_in_rhs< mutable_map_substitution<> > local_sigma;
            alphaconvert(stochvars,local_sigma, r1.variables(),data_expression_list());
            new_body2=substitute_pCRLproc(new_body2, local_sigma);
            const data_expression new_distribution=
                   process::replace_variables_capture_avoiding_with_an_identifier_generator(r2.distribution(),local_sigma,fresh_identifier_generator);

            return stochastic_operator(r1.variables() + stochvars,
                                       if_(if_then_else(t).condition(),
                                               lazy::and_(variables_are_equal_to_default_values(stochvars),r1.distribution()),
                                               lazy::and_(variables_are_equal_to_default_values(r1.variables()),new_distribution)),
                                       if_then_else(if_then_else(t).condition(),new_body1,new_body2));
          }
          /* r1 is and r2_ is not a stochastic operator */
          return stochastic_operator(r1.variables(),
                                     if_(if_then_else(t).condition(),r1.distribution(),
                                               if_(variables_are_equal_to_default_values(r1.variables()),sort_real::real_one(),sort_real::real_zero())),
                                     if_then_else(if_then_else(t).condition(),r1.operand(),r2_));
        }
        if (is_stochastic_operator(r2_))
        {
          /* r1_ is not and r2_ is a stochastic operator */
          const stochastic_operator& r2=down_cast<const stochastic_operator>(r2_);
          return stochastic_operator(r2.variables(),
                                     if_(if_then_else(t).condition(),
                                         if_(variables_are_equal_to_default_values(r2.variables()),sort_real::real_one(),sort_real::real_zero()),
                                         r2.distribution()),
                                     if_then_else(if_then_else(t).condition(),r1_,r2.operand()));
        }
        /* neither r1_ nor r2_ is stochastic */
        return t;

      }

      if (is_sum(t))
      {
        const sum& s=down_cast<const sum>(t);
        const process_expression r_= obtain_initial_distribution_term(s.operand());
        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          std::set <variable> variables_in_distribution=find_all_variables(r.distribution());
          for(const variable& v: s.variables())
          {
            if (variables_in_distribution.count(v)>0)
            {
              throw mcrl2::runtime_error("Cannot commute a sum operator over a stochastic operator in " +
                                                   process::pp(t) + ".\n" +
                                         "The problematic variable is " + pp(v) + ":" + pp(v.sort()) + ".");
            }
          }
          return stochastic_operator(r.variables(),
                                     r.distribution(),
                                     sum(s.variables(),r.operand()));
        }
        return t;
      }

      if (is_stochastic_operator(t))
      {
        const stochastic_operator& sto=down_cast<const stochastic_operator>(t);
        const process_expression r_= obtain_initial_distribution_term(sto.operand());

        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(sto.variables()+r.variables(),
                                     real_times_optimized(sto.distribution(),r.distribution()),
                                     r.operand());
        }
        return stochastic_operator(sto.variables(),sto.distribution(),r_);
      }

      if (is_action(t))
      {
        return t;
      }

      if (is_delta(t))
      {
        return t;
      }

      if (is_tau(t))
      {
        return t;
      }

      if (is_at(t))
      {
        const process_expression r_=obtain_initial_distribution_term(at(t).operand());
        if (is_stochastic_operator(r_))
        {
          const stochastic_operator& r=down_cast<const stochastic_operator>(r_);
          return stochastic_operator(r.variables(),r.distribution(),at(r.operand(),at(t).time_stamp()));
        }
        return t;
      }

      if (is_sync(t))
      {
        return t;
      }

      throw mcrl2::runtime_error("unexpected process format in canterminate " + process::pp(t) +".");
    }


    process_expression obtain_initial_distribution(const process_identifier& procId)
    {
      objectdatatype& object=objectIndex(procId);
      const process_expression initial_distribution_=obtain_initial_distribution_term(object.processbody);
      if (!is_stochastic_operator(initial_distribution_))
      {
        return process_instance_assignment(procId,assignment_list());
      }
      const stochastic_operator& initial_distribution=down_cast<const stochastic_operator>(initial_distribution_);
      if (!is_process_instance_assignment(initial_distribution.operand()))
      {
        const process_identifier new_procId=
                 newprocess(object.parameters + initial_distribution.variables(),
                            initial_distribution.operand(),
                            pCRL, canterminatebody(initial_distribution.operand()), containstimebody(initial_distribution.operand()));

        return stochastic_operator(initial_distribution.variables(),
                                   initial_distribution.distribution(),
                                   process_instance_assignment(new_procId,assignment_list())); // TODO add correct assignment here.
      }
      return initial_distribution;
    }

    /***** determinewhetherprocessescanterminate(init); **********/

    bool canterminatebody(
      const process_expression& t,
      bool& stable,
      std::set < process_identifier >& visited,
      const bool allowrecursion)
    {
      if (is_merge(t))
      {
        /* the construction below is needed to guarantee that
           both subterms are recursively investigated */
        const bool r1=canterminatebody(process::merge(t).left(),stable,visited,allowrecursion);
        const bool r2=canterminatebody(process::merge(t).right(),stable,visited,allowrecursion);
        return r1&&r2;
      }

      if (is_process_instance_assignment(t))
      {
        const process_instance_assignment u(t);
        if (allowrecursion)
        {
          return (canterminate_rec(u.identifier(),stable,visited));
        }
        return objectIndex(u.identifier()).canterminate;
      }

      if (is_hide(t))
      {
        return (canterminatebody(hide(t).operand(),stable,visited,allowrecursion));
      }

      if (is_rename(t))
      {
        return (canterminatebody(process::rename(t).operand(),stable,visited,allowrecursion));
      }

      if (is_allow(t))
      {
        return (canterminatebody(allow(t).operand(),stable,visited,allowrecursion));
      }

      if (is_block(t))
      {
        return (canterminatebody(block(t).operand(),stable,visited,allowrecursion));
      }

      if (is_comm(t))
      {
        return (canterminatebody(comm(t).operand(),stable,visited,allowrecursion));
      }

      if (is_choice(t))
      {
        const bool r1=canterminatebody(choice(t).left(),stable,visited,allowrecursion);
        const bool r2=canterminatebody(choice(t).right(),stable,visited,allowrecursion);
        return r1||r2;
      }

      if (is_seq(t))
      {
        const bool r1=canterminatebody(seq(t).left(),stable,visited,allowrecursion);
        const bool r2=canterminatebody(seq(t).right(),stable,visited,allowrecursion);
        return r1&&r2;
      }

      if (is_if_then(t))
      {
        return canterminatebody(if_then(t).then_case(),stable,visited,allowrecursion);
      }

      if (is_if_then_else(t))
      {
        const bool r1=canterminatebody(if_then_else(t).then_case(),stable,visited,allowrecursion);
        const bool r2=canterminatebody(if_then_else(t).else_case(),stable,visited,allowrecursion);
        return r1||r2;
      }

      if (is_sum(t))
      {
        return (canterminatebody(sum(t).operand(),stable,visited,allowrecursion));
      }

      if (is_stochastic_operator(t))
      {
        return (canterminatebody(stochastic_operator(t).operand(),stable,visited,allowrecursion));
      }

      if (is_action(t))
      {
        return true;
      }

      if (is_delta(t))
      {
        return false;
      }

      if (is_tau(t))
      {
        return true;
      }

      if (is_at(t))
      {
        return canterminatebody(at(t).operand(),stable,visited,allowrecursion);
      }

      if (is_sync(t))
      {
        const bool r1=canterminatebody(process::sync(t).left(),stable,visited,allowrecursion);
        const bool r2=canterminatebody(process::sync(t).right(),stable,visited,allowrecursion);
        return r1&&r2;
      }

      throw mcrl2::runtime_error("unexpected process format in canterminate " + process::pp(t) +".");
      return false;
    }

    bool canterminate_rec(
      const process_identifier& procId,
      bool& stable,
      std::set < process_identifier >& visited)
    {
      objectdatatype& object=objectIndex(procId);

      if (visited.count(procId)==0)
      {
        visited.insert(procId);
        const bool ct=canterminatebody(object.processbody,stable,visited,true);
        if (object.canterminate!=ct)
        {
          object.canterminate=ct;
          stable=false;
        }
      }
      return (object.canterminate);
    }

    bool canterminatebody(const process_expression& t)
    {
      std::set < process_identifier > visited;
      bool stable=false;
      return canterminatebody(t,stable,visited,false);
    }

    void determinewhetherprocessescanterminate(const process_identifier& procId)
    {
      bool stable=false;
      while (!stable)
      {
        std::set < process_identifier > visited;
        stable=true;
        canterminate_rec(procId,stable,visited);
      }
    }

    /*****  distinguishmCRLandpCRLprocsAndAddTerminatedAction  ******/

    process_identifier split_process(const process_identifier& procId,
                                     std::map < process_identifier,process_identifier >& visited_id,
                                     std::map < process_expression,process_expression >& visited_proc)
    {
      if (visited_id.count(procId)>0)
      {
        return visited_id[procId];
      }

      objectdatatype& object=objectIndex(procId);

      if ((object.processstatus!=mCRL) &&
          (object.canterminate==0))
      {
        /* no new process needs to be constructed */
        return procId;
      }

      if (object.processstatus==mCRL)
      {
        visited_id[procId]=procId;
        object.processbody = split_body(object.processbody,
                                               visited_id,visited_proc,
                                               object.parameters);
        return procId;
      }

      const process_identifier newProcId(
                   fresh_identifier_generator(procId.name()),
                   object.parameters);
      visited_id[procId]=newProcId;

      if (object.canterminate)
      {
        assert(check_valid_process_instance_assignment(terminatedProcId,assignment_list()));
        insert_process_declaration(
          newProcId,
          object.parameters,
          seq(object.processbody, process_instance_assignment(terminatedProcId,assignment_list())),
          pCRL,canterminatebody(object.processbody),
          containstimebody(object.processbody));
        return newProcId;
      }
      return procId;
    }

/* expand_process_instance_assignment takes a process instance assignment X(...) and
 * replaces it by Y if X=Y. This process is repeated until the right hand side of
 * an equation is not a process instance anymore */

    process_instance_assignment expand_process_instance_assignment(
                                 const process_instance_assignment& t,
                                 std::set<process_identifier>& visited_processes)
    {
      if (visited_processes.count(t.identifier())>0)
      {
        throw mcrl2::runtime_error("Process " + pp(t.identifier()) + " is unguardedly defined in itself");
      };

      visited_processes.insert(t.identifier());
      objectdatatype& object=objectIndex(t.identifier());
      if (is_process_instance_assignment(object.processbody))
      {
        const process_instance_assignment q=expand_process_instance_assignment(
                                                 down_cast<process_instance_assignment>(object.processbody),
                                                 visited_processes);

        maintain_variables_in_rhs<mutable_map_substitution<> > sigma;
        process_instance_assignment assign(t);
        for (const assignment& a: assign.assignments())
        {
          sigma[a.lhs()]=a.rhs();
        }

        return down_cast<process_instance_assignment>(substitute_pCRLproc(q,sigma));
      }
      return t;
    }

    process_instance_assignment expand_process_instance_assignment(const process_instance_assignment& t)
    {
      std::set<process_identifier> visited_processes;
      return expand_process_instance_assignment(t,visited_processes);
    }

/* Transform process_arguments
 *   This function replaces process_instances by process_instance_assignments.
 *   All assignments in a process_instance_assignment are ordered in the same
 *   sequence as the parameters belonging to that assignment.
 *   All assignments in a process_instance_assignment of the form x=x where
 *   x is not a bound variable are removed.
 *   Furthermore, process occurrences X where X is defined as X=X1, X1=X2...Xn-1=Xn
 *   are replaced by Xn, with the necessary substitutions for data parameters.
 *   The reason for this is that only Xn will become a process in the linear
 *   process, with its own state. Otherwise there is a risk that all process
 *   variables X, X1, X2, etc. have separate states, and worse, for all of them
 *   a copy of the summands will be added.
*/


    /* This function replaces all process instances by a process instance assignment,
       furthermore, if a process consists of only a process instantionation, i.e., X=Y, it is removed.
       Nested stochastic operators are merged into one. */
    void transform_process_arguments(
            const process_identifier& procId,
            std::set<process_identifier>& visited_processes)
    {
      if (visited_processes.count(procId)==0)
      {
        visited_processes.insert(procId);
        objectdatatype& object=objectIndex(procId);
        const std::set<variable> bound_variables;
        object.processbody=transform_process_arguments_body(
                                         object.processbody,
                                         bound_variables,
                                         visited_processes);
      }
    }

    /* This function replaces all process instances by a process instance assignment and
       merges nested stochastic operators into one.  */
    void transform_process_arguments(const process_identifier& procId)
    {
      std::set<process_identifier> visited_processes;
      transform_process_arguments(procId,visited_processes);
    }

    /* This function replaces all process instances by a process instance assignment
       and merges nested stochastic operators into one. */
    process_expression transform_process_arguments_body(
      const process_expression& t,
      const std::set<variable>& bound_variables,
      std::set<process_identifier>& visited_processes)
    {
      if (is_process_instance(t))
      {
        transform_process_arguments(process_instance(t).identifier(),visited_processes);
        process_instance_assignment t1=transform_process_instance_to_process_instance_assignment(atermpp::down_cast<process_instance>(t),bound_variables);
        return expand_process_instance_assignment(t1);
      }
      if (is_process_instance_assignment(t))
      {
        transform_process_arguments(process_instance_assignment(t).identifier(),visited_processes);
        const process_instance_assignment u(t);
        objectdatatype& object=objectIndex(u.identifier());
        assert(check_valid_process_instance_assignment(u.identifier(),
                 sort_assignments(u.assignments(),object.parameters)));
        process_instance_assignment t1(u.identifier(),
                                       sort_assignments(u.assignments(),object.parameters));
        return expand_process_instance_assignment(t1);
      }
      if (is_hide(t))
      {
        return hide(hide(t).hide_set(),
                    transform_process_arguments_body(hide(t).operand(),bound_variables,visited_processes));
      }
      if (is_rename(t))
      {
        return process::rename(
                 process::rename(t).rename_set(),
                 transform_process_arguments_body(process::rename(t).operand(),bound_variables,visited_processes));
      }
      if (is_allow(t))
      {
        return allow(allow(t).allow_set(),
                     transform_process_arguments_body(allow(t).operand(),bound_variables,visited_processes));
      }
      if (is_block(t))
      {
        return block(block(t).block_set(),
                     transform_process_arguments_body(block(t).operand(),bound_variables,visited_processes));
      }
      if (is_comm(t))
      {
        return comm(comm(t).comm_set(),
                    transform_process_arguments_body(comm(t).operand(),bound_variables,visited_processes));
      }
      if (is_merge(t))
      {
        return merge(
                 transform_process_arguments_body(merge(t).left(),bound_variables,visited_processes),
                 transform_process_arguments_body(merge(t).right(),bound_variables,visited_processes));
      }
      if (is_choice(t))
      {
        return choice(
                 transform_process_arguments_body(choice(t).left(),bound_variables,visited_processes),
                 transform_process_arguments_body(choice(t).right(),bound_variables,visited_processes));
      }
      if (is_seq(t))
      {
        return seq(
                 transform_process_arguments_body(seq(t).left(),bound_variables,visited_processes),
                 transform_process_arguments_body(seq(t).right(),bound_variables,visited_processes));
      }
      if (is_if_then_else(t))
      {
        return if_then_else(
                 if_then_else(t).condition(),
                 transform_process_arguments_body(if_then_else(t).then_case(),bound_variables,visited_processes),
                 transform_process_arguments_body(if_then_else(t).else_case(),bound_variables,visited_processes));
      }
      if (is_if_then(t))
      {
        return if_then(
                 if_then(t).condition(),
                 transform_process_arguments_body(if_then(t).then_case(),bound_variables,visited_processes));
      }
      if (is_sum(t))
      {
        std::set<variable> bound_variables1=bound_variables;
        const variable_list sum_vars=sum(t).variables();
        bound_variables1.insert(sum_vars.begin(),sum_vars.end());
        return sum(
                 sum_vars,
                 transform_process_arguments_body(sum(t).operand(),bound_variables1,visited_processes));
      }
      if (is_action(t))
      {
        return t;
      }
      if (is_delta(t))
      {
        return t;
      }
      if (is_tau(t))
      {
        return t;
      }
      if (is_at(t))
      {
        return at(
                 transform_process_arguments_body(at(t).operand(),bound_variables,visited_processes),
                 at(t).time_stamp());
      }
      if (is_sync(t))
      {
        return process::sync(
                 transform_process_arguments_body(process::sync(t).left(),bound_variables,visited_processes),
                 transform_process_arguments_body(process::sync(t).right(),bound_variables,visited_processes));
      }
      if (is_stochastic_operator(t))
      {
        const stochastic_operator& tso=down_cast<const stochastic_operator>(t);
        // Join nested stochastic operators into one.
        variable_list vars = tso.variables();
        data_expression dist = tso.distribution();
        process_expression body = tso.operand();
        while (is_stochastic_operator(body))
        {
          const stochastic_operator& tso1 = down_cast<const stochastic_operator>(body);
          vars = vars + tso1.variables();
          dist = sort_real::times(dist, tso1.distribution());
          body = tso1.operand();
        }
        return stochastic_operator(
                 vars,
                 dist,
                 transform_process_arguments_body(body,bound_variables,visited_processes));
      }
      throw mcrl2::runtime_error("unexpected process format in transform_process_arguments_body " + process::pp(t) +".");
    }

/* -----------------------------   split body  --------------------------- */

    process_expression split_body(
      const process_expression& t,
      std::map < process_identifier,process_identifier >& visited_id,
      std::map < process_expression,process_expression>& visited_proc,
      const variable_list& parameters)
    {
      /* Replace pCRL process terms that occur in the scope of mCRL processes
         by a process identifier. E.g. (a+b)||c is replaced by X||c and
         a new process equation X=a+b is added. Furthermore, if the replaced
         process can terminate a termination action is put behind it.
         In the example X=(a+b).terminate.delta@0.

         Besides this each ProcessAssignment is transformed into a Process. */

      process_expression result;

      if (visited_proc.count(t)>0)
      {
        return visited_proc[t];
      }

      if (is_merge(t))
      {
        result=process::merge(
                 split_body(process::merge(t).left(),visited_id,visited_proc,parameters),
                 split_body(process::merge(t).right(),visited_id,visited_proc,parameters));
      }
      else if (is_process_instance_assignment(t))
      {
        const process_instance_assignment u(t);
        objectdatatype& object=objectIndex(u.identifier());
        assert(check_valid_process_instance_assignment(split_process(u.identifier(),visited_id,visited_proc),
                 sort_assignments(u.assignments(),object.parameters)));
        result=process_instance_assignment(
                 split_process(u.identifier(),visited_id,visited_proc),
                 sort_assignments(u.assignments(),object.parameters));
      }
      else if (is_hide(t))
      {
        result=hide(hide(t).hide_set(),
                    split_body(hide(t).operand(),visited_id,visited_proc,parameters));
      }
      else if (is_rename(t))
      {
        result=process::rename(
                 process::rename(t).rename_set(),
                 split_body(process::rename(t).operand(),visited_id,visited_proc,parameters));
      }
      else if (is_allow(t))
      {
        result=allow(allow(t).allow_set(),
                     split_body(allow(t).operand(),visited_id,visited_proc,parameters));
      }
      else if (is_block(t))
      {
        result=block(block(t).block_set(),
                     split_body(block(t).operand(),visited_id,visited_proc,parameters));
      }
      else if (is_comm(t))
      {
        result=comm(comm(t).comm_set(),
                    split_body(comm(t).operand(),visited_id,visited_proc,parameters));
      }
      else if (is_stochastic_operator(t))
      {
        const stochastic_operator& t1 = atermpp::down_cast<stochastic_operator>(t);
        result=stochastic_operator(t1.variables(),
                                   t1.distribution(),
                                   split_body(t1.operand(),visited_id,visited_proc,parameters));
      }
      else if (is_choice(t)||
               is_seq(t)||
               is_if_then_else(t)||
               is_if_then(t)||
               is_sum(t)||
               is_action(t)||
               is_delta(t)||
               is_tau(t)||
               is_at(t)||
               is_sync(t))
      {
        if (canterminatebody(t))
        {
          assert(check_valid_process_instance_assignment(terminatedProcId,assignment_list()));
          const process_identifier p=newprocess(parameters,
                                                seq(t,process_instance_assignment(terminatedProcId,assignment_list())),
                                                pCRL,
                                                false,
                                                containstimebody(t));
          assert(check_valid_process_instance_assignment(p,assignment_list()));
          result=process_instance_assignment(p,assignment_list());
          visited_proc[t]=result;
        }
        else
        {
          const process_identifier p=newprocess(parameters,t,pCRL,false,containstimebody(t));
          assert(check_valid_process_instance_assignment(p,assignment_list()));
          result=process_instance_assignment(p,assignment_list());
          visited_proc[t]=result;
        }
      }
      else
      {
        throw mcrl2::runtime_error("unexpected process format in split process " + process::pp(t) +".");
      }

      return result;
    }

    process_identifier splitmCRLandpCRLprocsAndAddTerminatedAction(
                             const process_identifier& procId)
    {
      std::map < process_identifier,process_identifier> visited_id;
      std::map < process_expression,process_expression> visited_proc;
      return split_process(procId,visited_id,visited_proc);
    }

    /**************** AddTerminationActionIfNecessary ****************/

    void AddTerminationActionIfNecessary(const stochastic_action_summand_vector& summands)
    {
      for (const stochastic_action_summand& smd: summands)
      {
        const action_list multiaction=smd.multi_action().actions();
        if (multiaction == action_list({ terminationAction }))
        {
          acts.push_front(terminationAction.label());
          mCRL2log(mcrl2::log::warning) << "The action " << process::pp(terminationAction) <<
                           " followed by a deadlock is added to signal termination of the linear process. \n";
          return;
        }
      }
    }

  public:
    /********************** SieveProcDataVars ***********************/

    variable_list SieveProcDataVarsSummands(
      const std::set <variable>& vars,
      const stochastic_action_summand_vector& action_summands,
      const deadlock_summand_vector& deadlock_summands,
      const variable_list& parameters)
    {
      /* In this routine it is checked which free variables
         in vars occur in the summands. Those variables
         that occur in the summands are returned. The
         parameters are needed to check occurrences of vars
         in the assignment list */

      std::set < variable > vars_set(vars.begin(),vars.end());
      std::set < variable > vars_result_set;

      for (const deadlock_summand& smd: deadlock_summands)
      {
        if (smd.deadlock().has_time())
        {
          filter_vars_by_term(smd.deadlock().time(),vars_set,vars_result_set);
        }
        filter_vars_by_term(smd.condition(),vars_set,vars_result_set);
      }
      for (const stochastic_action_summand& smd: action_summands)
      {
        filter_vars_by_multiaction(smd.multi_action().actions(),vars_set,vars_result_set);
        filter_vars_by_assignmentlist(smd.assignments(),parameters,vars_set,vars_result_set);

        if (smd.multi_action().has_time())
        {
          filter_vars_by_term(smd.multi_action().time(),vars_set,vars_result_set);
        }
        filter_vars_by_term(smd.condition(),vars_set,vars_result_set);
        filter_vars_by_term(smd.distribution().distribution(),vars_set,vars_result_set);
      }
      variable_list result;
      for (std::set < variable >::reverse_iterator i=vars_result_set.rbegin();
           i!=vars_result_set.rend() ; ++i)
      {
        result.push_front(*i);
      }

      return result;
    }

  public:
    variable_list SieveProcDataVarsAssignments(
      const std::set <variable>& vars,
      const data_expression_list& initial_state_expressions)
    {
      const std::set < variable > vars_set(vars.begin(),vars.end());
      std::set < variable > vars_result_set;


      filter_vars_by_termlist(initial_state_expressions.begin(),
                              initial_state_expressions.end(),
                              vars_set,
                              vars_result_set);

      variable_list result;
      for (std::set < variable >::reverse_iterator i=vars_result_set.rbegin();
           i!=vars_result_set.rend() ; ++i)
      {
        result.push_front(*i);
      }

      return result;
    }

    /**************** transform **************************************/
  public:
    void transform(
      const process_identifier& init,
      stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      variable_list& parameters,
      data_expression_list& initial_state,
      stochastic_distribution& initial_stochastic_distribution)
    {
      /* Then select the BPA processes, and check that the others
         are proper parallel processes */
      transform_process_arguments(init);   // Also merges nested stochastic operators.
      determine_process_status(init,mCRL);
      determinewhetherprocessescanterminate(init);
      const process_identifier init_=splitmCRLandpCRLprocsAndAddTerminatedAction(init);
      determinewhetherprocessescontaintime(init_);

      std::vector <process_identifier> pcrlprocesslist;
      collectPcrlProcesses(init_,pcrlprocesslist);

      if (pcrlprocesslist.size()==0)
      {
        throw mcrl2::runtime_error("There are no processes to be linearised. This is most likely due to the use of unguarded recursion in process equations");
        // Note that this can occur with a specification
        // proc P(x:Int) = P(x); init P(1);
      }

      /* Second, transform into GNF with possibly variables as a head,
         but no actions in the tail */
      procstovarheadGNF(pcrlprocesslist);

      /* Third, transform to GNF by subsitution, such that the
         first variable in a sequence is always an actionvariable */
      procstorealGNF(init_,options.lin_method!=lmStack);

      lps::detail::ultimate_delay dummy_ultimate_delay_condition;
      generateLPEmCRL(action_summands,
                      deadlock_summands,
                      init_,
                      options.lin_method!=lmStack,
                      parameters,
                      initial_state,
                      initial_stochastic_distribution,
                      dummy_ultimate_delay_condition);
      allowblockcomposition(action_name_multiset_list({action_name_multiset()}),false,action_summands,deadlock_summands,terminationAction,options.ignore_time,options.nodeltaelimination); // This removes superfluous delta summands.
      if (options.final_cluster)
      {
        cluster_actions(action_summands,deadlock_summands,parameters);
      }

      AddTerminationActionIfNecessary(action_summands);
    }

}; // End of the class specification basictype

/**************** linearise **************************************/

mcrl2::lps::stochastic_specification mcrl2::lps::linearise(
  const mcrl2::process::process_specification& type_checked_spec,
  const mcrl2::lps::t_lin_options& lin_options)
{
  mCRL2log(mcrl2::log::verbose) << "linearising the process specification using the '" << lin_options.lin_method << "' method.\n";
  mcrl2::process::process_specification input_process=type_checked_spec;
  data_specification data_spec=input_process.data();

  if (lin_options.balance_summands) // Make a balanced tree of long expressions of the shape p1 + p2 + p3 + ... + p4.
                                    // By default the parser provides a skewed tree, and for very long sequences of summands this overflows the
                                    // stack.
  {
    balance_summands(input_process);
  }

  if (lin_options.apply_alphabet_axioms) // Apply alphabet reduction if requested.
  {
    alphabet_reduce(input_process, 1000ul);
  }

  std::set<data::sort_expression> s;
  process::find_sort_expressions(input_process.action_labels(), std::inserter(s, s.end()));
  process::find_sort_expressions(input_process.equations(), std::inserter(s, s.end()));
  process::find_sort_expressions(input_process.init(), std::inserter(s, s.end()));
  s.insert(sort_real::real_());
  data_spec.add_context_sorts(s);

  specification_basic_type spec(input_process.action_labels(),
                                input_process.equations(),
                                data::variable_list(input_process.global_variables().begin(),input_process.global_variables().end()),
                                data_spec,
                                input_process.global_variables(),
                                lin_options,
                                input_process);
  process_identifier init=spec.storeinit(input_process.init());

  //linearise spec
  variable_list parameters;
  data_expression_list initial_state;
  stochastic_action_summand_vector action_summands;
  deadlock_summand_vector deadlock_summands;
  stochastic_distribution initial_distribution(
                              variable_list(),
                              sort_real::real_one());
  spec.transform(init,action_summands,deadlock_summands,parameters,initial_state,initial_distribution);

  // compute global variables
  data::variable_list globals1 = spec.SieveProcDataVarsSummands(spec.global_variables,action_summands,deadlock_summands,parameters);
  data::variable_list globals2 = spec.SieveProcDataVarsAssignments(spec.global_variables,initial_state);
  std::set<data::variable> global_variables;
  global_variables.insert(globals1.begin(), globals1.end());
  global_variables.insert(globals2.begin(), globals2.end());
  stochastic_linear_process lps(parameters,
                                deadlock_summands,
                                action_summands);
  lps::stochastic_specification spec1(
      spec.data,
      spec.acts,
      global_variables,
      lps,
      stochastic_process_initializer(initial_state,initial_distribution));

  // add missing sorts to the data specification
  lps::complete_data_specification(spec1);
  return spec1;
}
