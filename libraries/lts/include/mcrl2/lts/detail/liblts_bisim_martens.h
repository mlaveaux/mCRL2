// Author(s): Jan Martens
//
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
/// \file lts/detail/liblts_bisim_m.h
///
/// \brief Partition refinement for Branching bisimularity reduction.
/// Inspired by Jules Jacobs, Thorsten Wi�mann,  "Fast coalgebraic bisimilarity minimization." - POPL2023.
///
/// \details 
#ifndef _LIBLTS_BISIM_MARTENS
#define _LIBLTS_BISIM_MARTENS
#include <fstream>
#include "mcrl2/modal_formula/state_formula.h"
#include "mcrl2/lts/lts_utilities.h"
#include "mcrl2/lts/detail/liblts_scc.h"
#include "mcrl2/lts/detail/liblts_merge.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/lts/lts_dot.h"
#include "mcrl2/utilities/execution_timer.h"
#include <boost/functional/hash.hpp>
#define UORDERED

namespace mcrl2
{
namespace lts
{
namespace detail
{
template < class LTS_TYPE>
class bisim_partitioner_martens
{
public:
  /** \brief Creates a bisimulation partitioner for an LTS.
    *  \details Based on the paper "Fast coalgebraic bisimilarity minimization." - Jacobs, Jules, and Thorsten Wi�mann. Proceedings of the ACM on Programming Languages 7.POPL (2023): 1514-1541.
    *  \warning Experimental.
    *  \param[in] l Reference to the LTS. */
  bisim_partitioner_martens(
    LTS_TYPE& l)
    : max_state_index(0),
    aut(l)
  {

    const std::vector<transition>& trans = aut.get_transitions();
    // Time this operation 
    auto startsort = std::chrono::high_resolution_clock::now();
    std::sort(aut.get_transitions().begin(), aut.get_transitions().end());
    auto stopsort = std::chrono::high_resolution_clock::now();
    auto sortduration = std::chrono::duration_cast<std::chrono::milliseconds>(stopsort - startsort);
    std::cout << "sort:" << sortduration.count() << std::endl;
    

    auto startinit = std::chrono::high_resolution_clock::now();
    //Initialize arrays for pred and suc, 
    pred = new std::vector<custom_transition_type>[aut.num_states()];
    suc = new std::vector<custom_transition_type>[aut.num_states()];
    trans_part = std::vector<trans_info>(aut.num_states());

    //TODO optimization: Derive these by partitioning the transition array
    //sil_pred = new std::set<state_type>[aut.num_states()];
    //sil_suc = new std::set<state_type>[aut.num_states()];

    blocks = new block_type[aut.num_states()];
    loc2state = new state_type[aut.num_states()];
    state2loc = new location_type[aut.num_states()];
    std::vector<std::size_t> state2in = std::vector<std::size_t>(aut.num_states(), 0);
    std::vector<std::size_t> state2out = std::vector<std::size_t>(aut.num_states(), 0);

    mCRL2log(mcrl2::log::debug) << "start moving transitions " << std::endl;
    for (auto r = trans.begin(); r != trans.end(); r++)
    {
      mCRL2log(mcrl2::log::debug) << r->to() << "- " << r->label() << " ->" << r->to() << std::endl;
      state2in[(*r).to()] += 1;
      state2out[(*r).from()] += 1;
    }

    //Count transitions per state
    for (state_type s = 0; s < aut.num_states(); ++s)
    {
      pred[s] = std::vector<custom_transition_type>(state2in[s]);
      suc[s] = std::vector<custom_transition_type>(state2out[s]);
      trans_part[s] = {0, state2in[s], 0, state2out[s], 0};

      //sil_pred[s] = std::set<state_type>();
      //sil_suc[s] = std::set<state_type>();
      state2loc[s] = s;
      loc2state[s] = s;
      blocks[s] = 0;
    }

    auto middleinit = std::chrono::high_resolution_clock::now();
    auto durationminit = std::chrono::duration_cast<std::chrono::milliseconds>(middleinit - startinit);
    std::cout << "m_init:" << durationminit.count() << std::endl;


    for (auto r = trans.begin(); r != trans.end(); r++)
    {
      state_type from = r->from();
      state_type to = r->to();
      if (!is_tau(r->label()))
      {
        pred[to][trans_part[to].mid_pred] = std::make_pair((*r).label(), (*r).from());
        trans_part[to].mid_pred += 1;
        suc[from][trans_part[from].mid_suc] = std::make_pair((*r).label(), (*r).to());
        trans_part[from].mid_suc += 1;
      }
      else
      {
        pred[to][trans_part[to].silent_pred - 1] = std::make_pair((*r).label(), (*r).from());
        trans_part[to].silent_pred -= 1;
        suc[from][trans_part[from].silent_suc - 1] = std::make_pair((*r).label(), (*r).to());
        trans_part[from].silent_suc -= 1;
      }
    }

    //Initialize block map
    worklist = std::deque<block_type>();
    block_map = std::vector<block>();

    block b0;
    b0.end = aut.num_states();
    b0.start = 0;
    b0.mid = aut.num_states();
    b0.frontier = aut.num_states();

    block_map.push_back(b0);
    //initialize worklist
    // Mark dirty bottom states, TODO: do this in the first initialization loop
    for (state_type s = 0; s < aut.num_states(); s++)
    {
      trans_part[s].silent_pred = pred[s].size();
      trans_part[s].silent_suc = suc[s].size();
      if (trans_part[s].mid_suc == suc[s].size())
      {
        // Bottom state:
        mark_dirty(s);
      }
    }
    auto endinit = std::chrono::high_resolution_clock::now();
    auto durationinit = std::chrono::duration_cast<std::chrono::milliseconds>(endinit - startinit);
    std::cout << "init:" << durationinit.count() << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    // Iterate refinement 
    refine();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << duration.count() << ":" << block_map.size() << std::endl;
    //cleanup
    delete[] pred;
    delete[] suc;
    delete[] loc2state;
    delete[] state2loc; 
    delete[] blocks;
  }

  /** \brief Destroys this partitioner. */
  ~bisim_partitioner_martens() = default;

  /** \brief Gives the number of bisimulation equivalence classes of the LTS.
   *  \return The number of bisimulation equivalence classes of the LTS.
   */
  std::size_t num_eq_classes() const
  {
    return block_map.size();
  }


  /** \brief Gives the bisimulation equivalence class number of a state.
   *  \param[in] s A state number.
   *  \return The number of the bisimulation equivalence class to which \e s belongs. */
  std::size_t get_eq_class(const std::size_t s) const
  {
    return blocks[s];
  }

  /** \brief Returns whether two states are in the same bisimulation equivalence class.
   *  \param[in] s A state number.
   *  \param[in] t A state number.
   *  \retval true if \e s and \e t are in the same bisimulation equivalence class;
   *  \retval false otherwise. */
  bool in_same_class(const std::size_t s, const std::size_t t) const
  {
    return get_eq_class(s) == get_eq_class(t);
  }

private:
  typedef std::size_t state_type;
  typedef std::size_t label_type;
  typedef std::size_t location_type;
  typedef std::size_t block_type;
  typedef std::pair<label_type, block_type> observation_type;
  typedef std::pair<label_type, state_type> custom_transition_type;

  //Typedef for signature
#ifdef UORDERED
  typedef std::vector<observation_type> signature_type;
#else
  typedef std::set<observation_type> signature_type;
#endif

  state_type max_state_index;
  LTS_TYPE& aut;

  // Array of vectors of predeccessors and successors
  // Form: mid_pred[s] = n; pred[s] = [pred_1, pred_2, \dots, pred_n, sil_pred_1, sil_pred_2, \dots, sil_pred_m]
  std::vector<custom_transition_type>* pred;
  // Form: mid_suc[s] = n; pred[s] = [pred_1, pred_2, \dots, pred_n, sil_pred_1, sil_pred_2, \dots, sil_pred_m]
  std::vector<custom_transition_type>* suc;

  //Struct block with start, mid, end being pointers
  // [start, mid) is clean, [mid, end) is dirty
  // [mid, bottom) is dirty and bottom.
  struct block
  {
    //   [   non-dirty |  dirty    |    dirty-frontier ]
    // start          mid         frontier
    location_type start;
    location_type mid;
    location_type frontier;
    location_type end;

    //   state_type end;
    //
    //std::set<state_type> states;
    //std::set<state_type> dirty_states;
    //std::set<state_type> frontier;
  };

  struct trans_info
  {
    // pred[s]  [   non-silent |  silent |  silent_marked ]
    //                        mid      silent
    // suc[s]   [   non-silent |  silent |  silent_marked ]
    //                        mid      silent
    location_type mid_pred;
    location_type silent_pred;
    location_type mid_suc;
    location_type silent_suc;
    std::size_t   num_dirty_suc;
  };

  block_type* blocks;
  std::vector<block> block_map; //TODO: change to vector
  //std::set<state_type> frontier;
  std::vector<trans_info> trans_part;
  //std::map<state_type, std::size_t> state2numdirtysuc; // NOOP!!!

  state_type* loc2state;
  location_type* state2loc;

  std::deque<block_type> worklist;

  bool is_tau(label_type l) 
  { 
    return aut.is_tau(aut.apply_hidden_label_map(l)); 
  }


  void swap(location_type loc, location_type loc2)
  {
    state_type s = loc2state[loc];
    state_type s2 = loc2state[loc2];
    loc2state[loc] = s2;
    loc2state[loc2] = s;
    state2loc[s] = loc2;
    state2loc[s2] = loc;
  }

  void move2dirty(state_type s)
  {
    location_type loc = state2loc[s];
    block* B = &block_map[blocks[s]];
    if (loc >= B->mid)
    {
      //already dirty.
      return;
    }
    if (loc < B->start || loc >= B->end) {
      mCRL2log(mcrl2::log::error) << "state not in block " << s << ":" << state2loc[s] << std::endl;
      mCRL2log(mcrl2::log::error) << "block start:" << B->start << " mid:" << B->mid << " frontier:" << B->frontier << " end:" << B->end << std::endl;
      
      return;
    }
    //Swap mid-1 and s
    B->mid -= 1;
    location_type new_loc = B->mid;
    swap(loc, new_loc);
  }
  
  void move2frontier(state_type s)
  {
    location_type loc = state2loc[s];
    block* B = &block_map[blocks[s]];
    if (loc >= B->frontier)
    {
      //already in frontier
      return;
    }
    if (loc < B->start || loc >= B->end)
    {
      mCRL2log(mcrl2::log::error) << "state not in block m2f" << std::endl;
      return;
    }

    //Swap frontier-1 and s
    B->frontier -= 1;
    location_type new_loc = B->frontier;
    swap(loc, new_loc);
  }

  void movefrontier2clean(state_type s)
  {
    location_type loc = state2loc[s];
    block* B = &block_map[blocks[s]];
    if (loc < B->frontier)
    {
      //already clean
      return;
    }

    if (loc < B->start || loc >= B->end)
    {
      mCRL2log(mcrl2::log::error) << "state not in block mf2c" << std::endl;
      return;
    }

    //Swap frontier and s
    B->frontier += 1;
    location_type new_loc = B->frontier - 1;
    swap(loc, new_loc);
  }

  void mark_dirty_backwards_closure(state_type s)
  { 
    for (auto sit = trans_part[s].mid_pred; sit < pred[s].size(); sit++)
    {
      state_type spre = pred[s][sit].second;
      if (blocks[spre] == blocks[s])
      {
        trans_part[spre].num_dirty_suc += 1; 
        block* B = &block_map[blocks[spre]];
        location_type locs = state2loc[spre];
        if (B->mid <= locs)
        {
          // The state was already dirty, so we remove it from the frontier
          movefrontier2clean(spre);
        }
        else
        {
          move2dirty(spre);
          mark_dirty_backwards_closure(spre);
        }
      }
      else
      {
        // TODO: Is this necessary? It seems to be necessary (we can go silently through the large part?).
        //mark_dirty(spre);
      }
    }
  }

  void mark_dirty(state_type s)
  {
    //Add to frontier if no silent marked state
    block* B = &block_map[blocks[s]];
    
    if (B->mid > state2loc[s])
    {
      move2dirty(s);
      if (B->end - B->mid == 1 && B->end - B->start > 1) {
        // The state is the first dirty state, and the block has more than 1 state, hence we add it to the worklist. 
        worklist.push_back(blocks[s]);
      }
      //The state was not dirty, so we add it to the frontier and compute reverse closure.
      move2frontier(s);
      mark_dirty_backwards_closure(s);
    }
  }

  //Check if sig1 has a observation not in sig2. both signatures are sorted.
  bool is_sig_subset(signature_type& sig1, signature_type& sig2, block_type t)
  {
    auto it1 = sig1.begin();
    auto it2 = sig2.begin();
    while (it1 != sig1.end())
    {
      if (it1->second == t)
      {
        it1++;
        continue;
      }
      if (it2 == sig2.end())
      {
        return false;
      }
      if (*it1 < *it2)
      {
        return false;
      }
      if (*it1 == *it2)
      {
        it1++;
        it2++;
      }
      else
      {
        it2++;
      }
    }
    return true;
  }

  bool is_bottom(state_type s)
  {
    return (trans_part[s].silent_suc == suc[s].size());
  }

  //Signature of a state
  void sig(const state_type& s,
           signature_type& retsignature,
           std::vector<signature_type>& sigs,
           std::vector<block_type>& state2num,
      block& B
  , std::vector<block_type>& retnums_to_add)
  {
    for (auto sit = 0; sit < trans_part[s].mid_suc; sit++)
    {
      custom_transition_type& t = suc[s][sit];
      retsignature.push_back(std::make_pair(t.first, blocks[t.second]));
    }
    
    for (auto sit = trans_part[s].mid_suc; sit < suc[s].size(); sit++)
    {
      // Loop through tau actions
      custom_transition_type& t = suc[s][sit];
      // Union with signature of t if signature is known
      // TODO Double check, what if t.second is in same block, but was not dirty??.
      if (blocks[t.second] != blocks[s])
      {
        // If not silent add observation
        retsignature.push_back(std::make_pair(t.first, blocks[t.second]));
      }
      else
      {
        // If silent and dirty, add signature of t
        if (state2loc[t.second] >= B.mid)
        {
          retnums_to_add.push_back(block_map.size() + state2num[state2loc[t.second] - B.mid]);
          retsignature.push_back(std::make_pair(t.first, block_map.size() + state2num[state2loc[t.second] - B.mid]));
        } 
        else
        {
          retsignature.push_back(std::make_pair(t.first, blocks[t.second]));
        }
        // else ? silent to clean ( is this important )??!? need some theory for this. 
      }
    }
  }

  //Renumber states based on signature
  std::vector<block_type> renumber(std::vector<block_type>& state2num, block_type max, block_type B)
  {
    //Count occurences of each number:
    std::vector<std::size_t> block2size(max+1);
    block2size[0] = block_map[B].mid - block_map[B].start;

    std::size_t maxsize = block_map[B].mid - block_map[B].start; 
    std::size_t maxblock = 0;
    std::size_t oldnumblocks = block_map.size();
    std::vector<state_type> dirty_states(&loc2state[block_map[B].mid], &loc2state[block_map[B].end]);

    for (std::size_t i=0; i<dirty_states.size(); i++)
    {
      std::pair<state_type, block_type> s = std::make_pair(dirty_states[i], state2num[i]);
      block2size[s.second] += 1;
      if (block2size[s.second] > maxsize)
      {
        maxsize = block2size[s.second];
        maxblock = s.second;
      }
    }
    bool swap = (block2size[0] < block2size[maxblock]);
    bool remove = (block2size[0] == 0);
    //mCRL2log(mcrl2::log::debug) << "maxblock:" << maxblock << " maxsize:" << maxsize << std::endl;
    //mCRL2log(mcrl2::log::debug) << "cleanblocksize:" << block_map[B].mid - block_map[B].start << std::endl;
    
    if (swap)
    {
      block2size[maxblock] = block2size[0];
      block2size[0] = maxsize;
    }
    
    std::vector<std::size_t> block2sizesum(max+1,0);
    block2sizesum[0] = block_map[B].start + block2size[0];

    block_map[B].end = block_map[B].start + block2size[0];
    block_map[B].mid = block_map[B].start + block2size[0];
    block_map[B].frontier = block_map[B].start + block2size[0];
    //mCRL2log(mcrl2::log::debug) << "old block:" << 0 << "->" << B << "start : " << block_map[B].start
    //                            << " end : " << block_map[B].end << " mid: " << block_map[B].mid
    //                            << " frontier : " << block_map[B].frontier << std::endl;

    for (std::size_t i = 1; i < max+1; ++i)
    {
      if (block2size[i] == 0)
      {
        block2sizesum[i] = block2sizesum[i-1];
        continue;
      }
      block2sizesum[i] = block2size[i] + block2sizesum[i-1];
      // Create newblock
      block newblock;
      newblock.start = block2sizesum[i-1];
      newblock.end = block2sizesum[i];
      newblock.mid = block2sizesum[i];
      newblock.frontier = block2sizesum[i];
      //mCRL2log(mcrl2::log::debug) << "new block:" << i << "->" << block_map.size()
      //                            << "start : " << newblock.start << " end : " << newblock.end 
      //                            << " mid: " << newblock.mid << " frontier : " << newblock.frontier << std::endl;
      block_map.push_back(newblock);
    }

    // Move clean states if necessary
    if (swap && !remove)
    {
      mCRL2log(mcrl2::log::debug) << "swap" << std::endl;
      std::vector<state_type> clean_states(&loc2state[block_map[B].start], &loc2state[block_map[B].mid]);
      for (state_type s : clean_states)
      {
        //new block 
        block_type target = oldnumblocks + maxblock - 1;
        blocks[s] = target;
        /* mCRL2log(mcrl2::log::debug) << "cstate:" << s << " block:" << blocks[s]
                                         << "loc:" << block2sizesum[maxblock] - 1 << std::endl;*/
        location_type new_loc = block2sizesum[maxblock] - 1;
        block2sizesum[maxblock] -= 1;
        loc2state[new_loc] = s;
        state2loc[s] = new_loc;
      }
    }

    for (std::size_t i = 0; i < dirty_states.size(); i++)
    {
      std::pair<state_type, block_type> ssigs = std::make_pair(dirty_states[i], state2num[i]);
      state_type s = ssigs.first;
      block_type target = label2realblock(ssigs.second, B, maxblock, oldnumblocks, remove);
      blocks[s] = target;
      block_type sizetarget = (ssigs.second == maxblock) ? 0 : ssigs.second;
      /* mCRL2log(mcrl2::log::debug) << "state:" << s << " block:" << blocks[s] 
                                  << "loc:" << block2sizesum[sizetarget] - 1
                                  << std::endl;*/

      location_type new_loc = block2sizesum[sizetarget] - 1;
      block2sizesum[sizetarget] -= 1;
      loc2state[new_loc] = s;
      state2loc[s] = new_loc;
    }

    std::vector<block_type> newblocks;
    return newblocks;
  }

  //Split block based on signature
  void split(block_type bid)
  {
    block* B = &block_map[bid];
    if (B->frontier == B->end || B->end - B->start == 1)
    {
      //Block is clean (should not happen)
      mCRL2log(mcrl2::log::debug) << "block clean!?." << B->start << " , mid: " << B->mid << " frontier:"
        << B->frontier << " end: " << B->end << std::endl;
      return;
    }

#ifdef UORDERED
    std::unordered_map<signature_type, block_type, boost::hash<signature_type>> sig2num; // , signaturehash
    signature_type signature;
#else  
    std::map<signature_type, block_type> sig2block;
    signature_type signature;
#endif


    std::vector<signature_type> num2sig;
    std::vector<block_type> state2num(B->end - B->mid);
    block_type j = 1; //Always start at 1, since 0 is reserved for clean states.
    num2sig.push_back(signature_type()); // 0 is empty signature reserverd for old block
    for(location_type locs = B->end; locs > B->mid; locs--)
    {
      if (locs < B->frontier)
      {
        mCRL2log(mcrl2::log::error) << "state was not yet ready to be processed" << std::endl;
      }
      state_type s = loc2state[locs-1];
      /* mCRL2log(mcrl2::log::debug) << "loc:" << locs - 1 << " state: " << s
                                  << " block: " << blocks[s] << std::endl;*/

      signature.clear();
      signature.reserve(suc[s].size()); //Only silent actions, bench if this is faster?
      std::vector<state_type> retnums_to_add;
      retnums_to_add.reserve(suc[s].size());

      sig(s, signature, num2sig, state2num, (*B), retnums_to_add);

      std::sort(signature.begin(), signature.end());
      auto last = std::unique(signature.begin(), signature.end());
      signature.erase(last, signature.end());
      
      for (auto& t: retnums_to_add)
      {
        // Inductive signature trick.
        if (is_sig_subset(signature, num2sig[t-block_map.size()], t))
        {
          signature = num2sig[t-block_map.size()]; //Minus T ...
          break;
        }
      }

      auto ret = sig2num.insert(std::make_pair(signature, j));
      if (ret.second)
      {
        // New signature
        sig2num[signature] = j;
        num2sig.push_back(signature);
        j += 1; 
      }

      state2num[locs - B->mid - 1] = ret.first->second; 
      
      for (auto sit = trans_part[s].mid_pred; sit < pred[s].size(); sit++)
      {
        state_type spre = pred[s][sit].second;
        // Mark silent backwards states
        if (blocks[spre] == blocks[s])
        {
          trans_part[spre].num_dirty_suc -= 1;
          if (trans_part[spre].num_dirty_suc == 0)
          {
            // If all dirty successor states are marked, we move the state to frontier.
            move2frontier(spre);
          }
        }
      }
    }

    std::vector<block_type> newblocks = renumber(state2num, num2sig.size(), bid); 
    return;
  }

  block_type label2realblock(block_type label, block_type bid, block_type max_block,block_type oldnumblocks,  bool remove)
  {
    block_type block = (label == 0) ? max_block : (label==max_block) ? 0 : label;
    if (block == 0)
    {
      return bid;
    }
    return (remove && label > max_block) ? oldnumblocks + block - 2 : oldnumblocks + block - 1;
  }

  //Refine based on sigs
  void refine()
  {
    int iter = 0;
    int new_blocks = 0;
    int old_blocks = 0;
    while (!worklist.empty())
    {
      block_type b = worklist.back();
      worklist.pop_back();
      block_type old_blocks = block_map.size();
      split(b);
      // update dirty states.
      block_type new_blocks = block_map.size(); 
      if (new_blocks == old_blocks)
      {
        mCRL2log(mcrl2::log::debug) << "No split: old_blocks:" << old_blocks << " new_blocks:" << new_blocks << std::endl;
        mCRL2log(mcrl2::log::debug) << "block of size: " << block_map[b].end - block_map[b].start << std::endl;
      }

      for (block_type b = old_blocks; b < new_blocks; ++b)
      {
        std::vector<state_type> states(&loc2state[block_map[b].start], &loc2state[block_map[b].end]);

        for (auto s : states)
        {
          for (auto t : pred[s])
          {
            if (!is_tau(t.first) or blocks[t.second] != b)
            {
              mark_dirty(t.second);
            }
          }
        }
      }
      iter += 1;
    }
  }
};


template <class LTS_TYPE>
void bisimulation_reduce_martens(LTS_TYPE& l)
{
    if (1 >= l.num_states())
    {
        mCRL2log(log::warning) << "There is only 1 state in the LTS. It is not "
                "guaranteed that branching bisimulation minimisation runs in "
                "time O(m log n).\n";
    }

    // Line 1.2: Find tau-SCCs and contract each of them to a single state
    mcrl2::utilities::execution_timer timer;
    mCRL2log(log::verbose) << "Start SCC\n";
    timer.start("preprocess");
    scc_reduce(l, false);
    timer.finish("preprocess");

    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
    timer.start("reduction");
    bisim_partitioner_martens<LTS_TYPE> bisim_part(l);
    timer.finish("reduction");
    timer.report();
}


template <class LTS_TYPE>
bool destructive_bisimulation_compare_martens(LTS_TYPE& l1, LTS_TYPE& l2)
{
    std::size_t init_l2(l2.initial_state() + l1.num_states());
    detail::merge(l1, std::move(l2));
    l2.clear(); // No use for l2 anymore.

    // Line 2.1: Find tau-SCCs and contract each of them to a single state
    detail::scc_partitioner<LTS_TYPE> scc_part(l1);
    scc_part.replace_transition_system(false);
    init_l2 = scc_part.get_eq_class(init_l2);
    
    bisim_partitioner_martens<LTS_TYPE> bisim_part(l1);

    return bisim_part.in_same_class(l1.initial_state(), init_l2);
}

}
}
}


#endif

