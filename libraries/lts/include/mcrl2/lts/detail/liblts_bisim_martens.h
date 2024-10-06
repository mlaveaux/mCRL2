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

//#define UORDERED
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
    auto start = std::chrono::high_resolution_clock::now();

    const std::vector<transition>& trans = aut.get_transitions();

    //Initialize arrays for pred and suc, blocks and state2loc and loc2state
    pred = new std::vector<custom_transition_type>[aut.num_states()];
    suc = new std::vector<custom_transition_type>[aut.num_states()];

    //TODO optimization: Derive these by partitioning the transition array
    sil_pred = new std::set<state_type>[aut.num_states()];
    sil_suc = new std::set<state_type>[aut.num_states()];
    pre_marked = new std::set<state_type>[aut.num_states()];
    marked = new std::set<state_type>[aut.num_states()];

    blocks = new block_type[aut.num_states()];
    state2loc = new state_type[aut.num_states()];
    loc2state = new state_type[aut.num_states()];
    std::vector<int> state2in = std::vector<int>(aut.num_states(), 0);
    std::vector<int> state2out = std::vector<int>(aut.num_states(), 0);
    mCRL2log(mcrl2::log::debug) << "start moving transitions " << std::endl;

    //Count transitions per state
    for (auto r = trans.begin(); r != trans.end(); r++) {
      state2in[(*r).to()] += 1;
      state2out[(*r).from()] += 1;
    }

    for (state_type s = 0; s < aut.num_states(); ++s)
    {
      pred[s] = std::vector<custom_transition_type>(state2in[s]);
      suc[s] = std::vector<custom_transition_type>(state2out[s]);
      sil_pred[s] = std::set<state_type>();
      sil_suc[s] = std::set<state_type>();
      marked[s] = std::set<state_type>();
      pre_marked[s] = std::set<state_type>();
    }

    mCRL2log(mcrl2::log::debug) << "now placing in correct place " << std::endl;
    for (auto r = trans.begin(); r != trans.end(); r++)
    {
      state2in[(*r).to()] -= 1;
      pred[(*r).to()][state2in[(*r).to()]] = std::make_pair((*r).label(), (*r).from());
      state2out[(*r).from()] -= 1;
      suc[(*r).from()][state2out[(*r).from()]] = std::make_pair((*r).label(), (*r).to());
      if (is_tau((*r).label())) {
        sil_pred[(*r).to()].insert((*r).from());
        sil_suc[(*r).from()].insert((*r).to());
      }
    }
    mCRL2log(mcrl2::log::debug) << "moved all transitions" << std::endl;
    //Initialize block map
    worklist = std::queue<block_type>();
    block_map = std::map<block_type, block>();

    //Initialize blocks
    block b0 = block{ 0,(unsigned int)aut.num_states() ,(unsigned int)aut.num_states() };
    for (state_type i = 0; i < (unsigned int)aut.num_states(); ++i)
    {
      blocks[i] = 0;
    }

    block_map[0] = b0;
    //initialize worklist
    worklist.push(0);

    //Initialize state2loc and loc2state
    for (std::size_t i = 0; i < (unsigned int) aut.num_states(); ++i)
    {
      state2loc[i] = i;
      loc2state[i] = i;
    }
    // Mark dirty bottom states
    for (std::size_t i = 0; i < (unsigned int) aut.num_states(); ++i)
    {
      if (sil_suc[i].empty()) {
        mark_dirty(i);
      }
    }

    // Iterate refinement 
    refine();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    mCRL2log(mcrl2::log::info) << "Done s:" << duration.count() << std::endl;
    //cleanup
    delete[] pred;
    delete[] suc;
    delete[] state2loc;
    delete[] loc2state;
    delete[] blocks;
    delete[] sil_pred;
    delete[] sil_suc;
    delete[] marked;
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
    mCRL2log(mcrl2::log::debug) << "in_same_class " << s << " " << t << " " << get_eq_class(s) << " " << get_eq_class(t) << std::endl;
    return get_eq_class(s) == get_eq_class(t);
  }

private:
  typedef std::size_t state_type;
  typedef std::size_t label_type;
  typedef std::size_t block_type;
  typedef std::pair<label_type, block_type> observation_type;
  typedef std::pair<label_type, state_type> custom_transition_type;

  //Typedef for signature
#ifdef UORDERED
  typedef std::set<observation_type> signature_type;
#else
  typedef std::set<observation_type> signature_type;
#endif

  bool is_tau(label_type l)
  {
    return aut.is_tau(aut.apply_hidden_label_map(l));
  }

  state_type max_state_index;
  LTS_TYPE& aut;
  // Array of vectors of predeccessors and successors
  std::vector<custom_transition_type>* pred;
  std::vector<custom_transition_type>* suc;

  //Struct block with start, mid, end being pointers
  // [start, mid) is clean, [mid, end) is dirty
  // [mid, bottom) is dirty and bottom.
  struct block
  {
    state_type start;
    state_type mid;
    state_type end;
    state_type bottom;
    //std::unordered_set<label_type> unstable_labels;
  };

  struct sigHash
  {
    std::size_t operator()(signature_type const& s) const
    {
      std::size_t res = 0;
      for (auto t : s)
      {
        res = res ^ t.first;
        res = res ^ t.second;
      }
      return res;
    }
  };

  block_type* blocks;
  std::map<block_type, block> block_map;
  state_type* loc2state;
  state_type* state2loc;
  std::set<state_type>* sil_pred;
  std::set<state_type>* sil_suc;
  std::set<state_type>* pre_marked;
  std::set<state_type>* marked;

  //Implement this later.
  //std::size_t* state2silentout;
  //std::size_t* state2silentin;

  std::queue<block_type> worklist;

  void mark_dirty_silent_closure(state_type s, state_type origin)
  {
    if (pre_marked[s].empty()) {
      for (auto t : sil_pred[s]) {
        if (blocks[s] == blocks[t]) {
          mark_dirty_silent_closure(t, s);
        }
      }
      block B = block_map[blocks[s]];
      if (B.mid <= state2loc[s]) {
        // Was dirty, but now not bottom:
        swap_to_not_dirty(s);
      }
    }
    pre_marked[s].insert(origin);
  }

  void mark_dirty(state_type s)
  {
    block_type Bid = blocks[s];
    block B = block_map[Bid];
    state_type loc = state2loc[s];
    //Add label to unstable labels
    /*if (B.unstable_labels.find(a) == B.unstable_labels.end()) {
        block_map[Bid].unstable_labels.insert(a);
    }*/

    if (loc < B.start || loc > B.end) {
      mCRL2log(mcrl2::log::info) << "KAPUTt " << loc << " " << B.start << " " << B.mid << " " << B.end << std::endl;
    }

    if (B.mid <= loc or B.start >= B.end - 1) {
      //Already dirty or only 1 state
      return;
    }

    if (B.mid == B.end) {
      //First dirty state
      worklist.push(Bid);
    }

    if (pre_marked[s].empty()) {
      swap_to_dirty(s);
      for ( auto spre : sil_pred[s]) {
        if (blocks[spre] == blocks[s]) {
          mark_dirty_silent_closure(spre, s);
        }
      }
    }

    //Swap last clean state with s

  }

  void swap_to_dirty(state_type s)
  {
    mCRL2log(mcrl2::log::debug) << "swap_to_dirty " << s << std::endl;
    block_type Bid = blocks[s];
    if (block_map.find(Bid) == block_map.end()) {
      mCRL2log(mcrl2::log::info) << "Block not found !?" << s << ":" << Bid << std::endl;
    }
    block B = block_map[Bid];
    state_type loc = state2loc[s];
    if (loc < B.start || loc > B.end) {
      mCRL2log(mcrl2::log::info) << "KAPUT1 " << loc << " " << B.start << " " << B.mid << " " << B.end << std::endl;
    }
    if (B.mid <= loc or B.start >= B.end - 1) {
      mCRL2log(mcrl2::log::info) << "KAPUT2 " << loc << " " << B.start << " " << B.mid << " " << B.end << std::endl;
      return;
    }
    //Swap last clean state with s
    state_type new_id = B.mid - 1;
    state_type tmp = loc2state[new_id];
    block_map[Bid].mid = new_id;

    state2loc[tmp] = loc;
    loc2state[loc] = tmp;
    state2loc[s] = new_id;
    loc2state[new_id] = s;
  }

  void swap_to_not_dirty(state_type s)
  {
    block_type Bid = blocks[s];
    block B = block_map[Bid];
    state_type loc = state2loc[s];
    if (loc < B.start || loc > B.end) {
      mCRL2log(mcrl2::log::info) << "KAPUTb " << loc << " " << B.start << " " << B.mid << " " << B.end << std::endl;
    }
    if (B.mid > loc or B.start >= B.end - 1) {
      mCRL2log(mcrl2::log::info) << "KAPUTa " << loc << " " << B.start << " " << B.mid << " " << B.end << std::endl;
      return;
    }
    //Swap first dirty state with s
    state_type new_id = B.mid;
    state_type tmp = loc2state[new_id];
    block_map[Bid].mid = new_id+1;

    state2loc[tmp] = loc;
    loc2state[loc] = tmp;
    state2loc[s] = new_id;
    loc2state[new_id] = s;
  }

  //Signature of a state
  void sig(const state_type& s, signature_type& retsignature, std::map<state_type, signature_type>& sigs, block_type curblock)
  {
    for (auto t : suc[s])
    {
      if (blocks[t.second] != curblock or !is_tau(t.first)) {
        retsignature.insert(std::make_pair(t.first, blocks[t.second]));
        auto sig = std::make_pair(t.first, blocks[t.second]);
        // Add observation of r to sig
        retsignature.insert(sig);
      }
      else {
        // Silent tau, if it is marked, we should include the signature.
        if (marked[s].find(t.second) != marked[s].end()) {
          marked[s].erase(t.second);
          retsignature.insert(sigs[t.second].begin(), sigs[t.second].end());
        }
      }
    }
  }

  //Split block based on signature
  int split(block_type Bid)
  {
    block B = block_map[Bid];
    if (B.mid == B.end)
    {
      //Block is clean (should not happen)
      mCRL2log(mcrl2::log::debug) << "Block clean but in worklist!? should not happen." << std::endl;
      return 0;
    }
#ifdef UORDERED
    std::unordered_map<signature_type, block_type, sigHash> sig2block;
#else  
    std::map<signature_type, block_type> sig2block;
#endif
    int j = 0;
    std::vector<block_type> state2block = std::vector<block_type>(aut.num_states(), 0);

    mCRL2log(mcrl2::log::debug) << "Computing signatures from here. "<< Bid << ":" << B.start <<  " : " << B.mid << " : " << B.end << std::endl;

    //Add signature of one clean state
    signature_type signature;
    //TODO: This might hit complexity, maybe improve this by juggling references to correct signatures?
    std::map<state_type, signature_type> sigs;

    std::set<state_type> checksum;
    //Add signatures of dirty states
    for (state_type i = B.end; i > B.mid; i--)
    {
      state_type s = loc2state[i - 1];
      if(checksum.find(s) != checksum.end()) {
        mCRL2log(mcrl2::log::debug) << "Duplicate state: " << s << std::endl;
      }
      checksum.insert(s);
      
      mCRL2log(mcrl2::log::debug) << "Making sig " << s << ": " << state2loc[s] << std::endl;

      if(!pre_marked[s].empty()) 
      {
        mCRL2log (mcrl2::log::debug) << "Bugg: " << s << std::endl;
      }
      signature.clear();
      mCRL2log(mcrl2::log::debug) << "sig clear" << s << std::endl;

      sig(s, signature, sigs, Bid);
      mCRL2log(mcrl2::log::debug) << "sig comped" << s << std::endl;

      sigs[s] = signature;
      auto ret = sig2block.insert(std::make_pair(signature, j));
      if (ret.second) {
        j += 1;
      }
      mCRL2log(mcrl2::log::debug) << "silent closure and bookkeeping " << s << std::endl;

      state2block[s] = (*ret.first).second;
      for (auto t : sil_pred[s]) {
        mCRL2log(mcrl2::log::debug) << "1";
        if (blocks[t] == Bid) {
          mCRL2log(mcrl2::log::debug) << "2";
          pre_marked[t].erase(s);
          mCRL2log(mcrl2::log::debug) << "3";
          marked[t].insert(s);
          mCRL2log(mcrl2::log::debug) << "removing premark: " << t << ":" << loc2state[s] << std::endl;
          if(pre_marked[t].empty()) {
            swap_to_dirty(t);
          }
          mCRL2log(mcrl2::log::debug) << "5";
        }
      }
      B = block_map[Bid];
    }

    for (state_type i = B.mid; i < B.end; i++)
    {
      if(!pre_marked[loc2state[i]].empty()) 
      {
        mCRL2log (mcrl2::log::debug) << "Not empty premarked: " << loc2state[i] << std::endl;
      }
      if (!marked[loc2state[i]].empty())
      {
        mCRL2log(mcrl2::log::debug) << "Not empty marked: " << loc2state[i] << std::endl;
      }
    }

    mCRL2log(mcrl2::log::debug) << "Sigs computed. : " << B.start << ":" << B.mid << ":" << B.end << std::endl;
    for(state_type s= B.start; s < B.end; s++)
    {
      mCRL2log(mcrl2::log::debug) << ":" << loc2state[s] << std::endl;
      for(state_type sp : pre_marked[loc2state[s]]) {
       mCRL2log(mcrl2::log::debug) << "premarked " << sp << std::endl;
      }
    }

    if (j == 1) {
      //Only one signature, no need to split
      mCRL2log(mcrl2::log::debug) << "no new signatures." << B.end - B.start << " " << B.mid - B.start << std::endl;
      block_map[Bid].mid = B.end;
      return 0;
    }
    //Count number of occurrences each signature
    mCRL2log(mcrl2::log::debug) << "Sigs computed.2 :" << std::endl;

    int* Sizes = new int[j + 1];
    for (int i = 0; i <= j; i++)
    {
      Sizes[i] = 0;
    }
    //Sizes[0] = B.mid - B.start;
    for (state_type s = B.mid; s < B.end; s++)
    {
      mCRL2log(mcrl2::log::debug) << "Counting signatures:" << B.mid << std::endl;
      mCRL2log(mcrl2::log::debug) << ":" << s << ":" << Sizes[state2block[loc2state[s]]] << std::endl;
      Sizes[state2block[loc2state[s]]] += 1;
    }
    mCRL2log(mcrl2::log::debug) << "Create new blocks. :" << std::endl;
      
    //Create new blocks
    // argmax Sizes
    state_type max = 0;
    int max_index = 0;
    Sizes[0] += B.mid - B.start;

    for (int i = 0; i < j; i++)
    {
      if (Sizes[i] > max) {
        max = Sizes[i];
        max_index = i;
      }
    }
    Sizes[0] -= B.mid - B.start;

    //Prefix sum Sizes
    for (int i = 0; i < j; i++)
    {
      Sizes[i + 1] += Sizes[i];
    }
    mCRL2log(mcrl2::log::debug) << "Going to mark dirties. :" << std::endl;
    for (int i =0; i <= j ; i++)
    {
      mCRL2log(mcrl2::log::debug) << Sizes[i] << std::endl;
    }

    int num_dirty = B.end - B.mid;
    state_type* dirty = new state_type[num_dirty];
    std::copy(loc2state + B.mid, loc2state + B.end, dirty);
    //Reorder states

    mCRL2log(mcrl2::log::debug) << "Rearranging states." << std::endl;

    for (state_type i = 0; i < num_dirty; i++) {
      state_type si = dirty[i];
      if (Sizes[state2block[si]] == 0) {
        //Impossible
        mCRL2log(mcrl2::log::info) << "Impossible" << si << " " << state2block[si] << std::endl;
        continue;
      }
      Sizes[state2block[si]] -= 1;
      int tmp = B.mid + Sizes[state2block[si]];
      loc2state[tmp] = si;
      state2loc[si] = tmp;
    }
    delete[] dirty;
    // create new blocks
    state_type old_start = B.mid;
    state_type old_end = B.end;

    for (int i = 0; i < j; i++)
    {
      state_type newstart = old_start + Sizes[i];
      state_type newend = old_start + Sizes[i + 1];
      if (i == 0) {
        newstart = B.start;
      }

      assert(newstart >= old_start);
      assert(newend <= old_end);
      if (i == max_index)
      {
        B.start = newstart;
        B.mid = newend;
        B.end = newend;
        block_map[Bid] = B;
      } else {
        block_type newBid = block_map.size();
        block newB = block{ newstart, newend, newend };
        block_map[newBid] = newB;
        for (state_type locs = newB.start; locs < newB.end; locs++)
        {
          blocks[loc2state[locs]] = newBid;
        }
      }
    }
    mCRL2log(mcrl2::log::debug) << "blocks? " <<  block_map.size() <<  std::endl;
    // Done with rearranging states
    // Output partition for debug purposes.
    mCRL2log(mcrl2::log::debug) << "Done rearranging states." << std::endl;
    for (state_type i = 0; i < (unsigned int) aut.num_states(); i++)
    {
      mCRL2log(mcrl2::log::debug) << "\t loc:" << i << ":" << state2loc[loc2state[i]] <<":" << loc2state[i] << " " << blocks[loc2state[i]] << std::endl;
    }
    delete[] Sizes;
    /* delete[] Sizes;
    delete[] dirty;
    delete[] state2block;*/
    return j;
}

  //Refine based on sigs
  void refine()
  {
    int iter = 0;
    int new_blocks = 0;
    int old_blocks = 0;
    mCRL2log(mcrl2::log::info) << "Start refinement" << std::endl;
    while (!worklist.empty())
    {
      iter += 1;
      block_type Bid = worklist.front();
      worklist.pop();
      mCRL2log(mcrl2::log::debug) << "Sup!." << Bid << std::endl;

      block B = block_map[Bid];

      //Count states
      if (B.mid == B.end and B.end > B.start + 1)
      {
        //Block is clean (should not happen)
        mCRL2log(mcrl2::log::info) << "Block clean but in worklist!? should not happen." << std::endl;
      }
      else {
        mCRL2log(mcrl2::log::debug) << "start iter" << std::endl;

        old_blocks = block_map.size();
        int ret = split(Bid);
        mCRL2log(mcrl2::log::debug) << "new blocks created = " << block_map.size() - old_blocks << std::endl;


        for (auto b0 : block_map)
        {
          block B = b0.second;
          mCRL2log(mcrl2::log::debug) << "block {" << b0.first << ", start=" << B.start << "mid = " << B.mid << ", end = " << B.end << "}" << std::endl;
        }

        if(old_blocks != block_map.size()) {
          // Mark every bottom state dirty.
          for (std::size_t i = 0; i < (unsigned int) aut.num_states(); ++i)
          {
            bool bottom = true;
            for (auto t : sil_suc[i]) {
              if (blocks[t] == blocks[i]) {
                bottom = false;
                break;
              }
            }
            if (bottom) {
              mCRL2log(mcrl2::log::debug) << i << "bottom: " << i << ": B:" << blocks[i] << std::endl;
              mark_dirty(i);
            }
          }
        }
        

        //Temporrary disable
        if (false) {
          for (int blockid = old_blocks; blockid < block_map.size(); blockid++)
          {
            mCRL2log(mcrl2::log::debug) << "marking dirty " << blockid << std::endl;
            block newB = block_map[blockid];

            std::vector<observation_type> overlap;
            for (state_type s = newB.start; s < newB.end; s++) {
              for (custom_transition_type t : pred[loc2state[s]]) {
                if (blocks[t.second] == blockid) {
                  if (!is_tau(t.first)) {
                    overlap.push_back(t);
                  }
                }
                else {
                  mark_dirty(t.second);
                  mCRL2log(mcrl2::log::debug) << "marking dirty " << t.second << std::endl;
                }
              }
            }
            for (auto s : overlap) {
              mCRL2log(mcrl2::log::debug) << "marking dirty " << s.second << std::endl;
              mark_dirty(s.second);
            }
          }
        }
      }
    

      mCRL2log(mcrl2::log::debug) << "Done marking dirty now blocks:" << std::endl;
      for (auto b0 : block_map)
      {
        block B = b0.second;
        mCRL2log(mcrl2::log::debug) << "block {" << b0.first << ", start=" << B.start << ", mid = " << B.mid << ", end = " << B.end << "}" << std::endl;
        for (state_type i=B.start; i < B.end; i++) {
          mCRL2log(mcrl2::log::debug) << "\t state " << loc2state[i] << " " << blocks[loc2state[i]] << " "<< suc[loc2state[i]].size() << std::endl;
        }
      }
    }
    mCRL2log(mcrl2::log::info) << "Done total blocks: \"" << block_map.size() << "\"" << std::endl;
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
    mCRL2log(log::verbose) << "Start SCC\n";
    scc_reduce(l, false);

    // Now apply the branching bisimulation reduction algorithm.  If there
    // are no taus, this will automatically yield strong bisimulation.
    bisim_partitioner_martens<LTS_TYPE> bisim_part(l);
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

