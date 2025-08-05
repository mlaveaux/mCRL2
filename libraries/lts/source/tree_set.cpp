// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file tree_set.cpp

#include "mcrl2/utilities/detail/memory_utility.h"
#include "mcrl2/utilities/logger.h"
#include "mcrl2/lts/detail/tree_set.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

#define EMPTY_SET (-1)
#define EMPTY_LIST (-1)
#define EMPTY_TAG (-1)
#define HASH_CLASS 16
#define TAGS_BLOCK 15000
#define BUCKETS_BLOCK 25000
// simple hash function; uses two large primes
#define hash(l,r,m) ((36425657*(l) + 77673689*(r)) & (m))


namespace mcrl2::lts
{

tree_set_store::tree_set_store()
{
  buckets = nullptr;
  buckets_size = 0;
  buckets_next = 0;

  tags = nullptr;
  tags_size = 0;
  tags_next = 0;

  hashmask = (1 << HASH_CLASS) - 1;
  hashtable = (ptrdiff_t*)malloc((hashmask+1)*sizeof(ptrdiff_t));
  if (hashtable == nullptr)
  {
    throw mcrl2::runtime_error("Out of memory.");
  }
  for (ptrdiff_t i=0; i<=hashmask; ++i)
  {
    hashtable[i] = EMPTY_LIST;
  }
}

tree_set_store::~tree_set_store()
{
  if (tags != nullptr)
  {
    free(tags);
    tags = nullptr;
  }
  if (buckets != nullptr)
  {
    free(buckets);
    buckets = nullptr;
  }
  free(hashtable);
  hashtable = nullptr;
}

void tree_set_store::check_tags()
{
  if (tags_next >= tags_size)
  {
    tags_size += TAGS_BLOCK;
    tags = (ptrdiff_t*)realloc(tags,tags_size*sizeof(ptrdiff_t));
    if (tags == nullptr)
    {
      throw mcrl2::runtime_error("Out of memory.");
    }
  }
}

void tree_set_store::check_buckets()
{
  if (buckets_next >= buckets_size)
  {
    buckets_size += BUCKETS_BLOCK;
    buckets = (bucket*)realloc(buckets,buckets_size*sizeof(bucket));
    if (buckets == nullptr)
    {
      throw mcrl2::runtime_error("Out of memory.");
    }
  }
  if (buckets_next*4 >= hashmask*3)
  {
    hashmask = hashmask + hashmask + 1;
    hashtable = (ptrdiff_t*)realloc(hashtable,(hashmask+1)*sizeof(ptrdiff_t));
    if (hashtable == nullptr)
    {
      throw mcrl2::runtime_error("Out of memory.");
    }
    ptrdiff_t i;
    ptrdiff_t hc;
    for (i=0; i<=hashmask; ++i)
    {
      hashtable[i] = EMPTY_LIST;
    }
    for (i=0; i<buckets_next; ++i)
    {
      hc = hash(buckets[i].child_l,buckets[i].child_r,hashmask);
      buckets[i].next = hashtable[hc];
      hashtable[hc] = i;
    }
  }
}

ptrdiff_t tree_set_store::build_set(ptrdiff_t child_l,ptrdiff_t child_r)
{
  check_buckets();
  ptrdiff_t hc = hash(child_l,child_r,hashmask);
  buckets[buckets_next].child_l = child_l;
  buckets[buckets_next].child_r = child_r;
  buckets[buckets_next].tag     = EMPTY_TAG;
  buckets[buckets_next].next    = hashtable[hc];
  hashtable[hc] = buckets_next;
  return buckets_next++;
}

ptrdiff_t tree_set_store::find_set(ptrdiff_t child_l,ptrdiff_t child_r)
{
  ptrdiff_t hc = hash(child_l,child_r,hashmask);
  for (ptrdiff_t i=hashtable[hc]; i!=EMPTY_LIST; i=buckets[i].next)
  {
    if (buckets[i].child_l==child_l && buckets[i].child_r==child_r)
    {
      return i;
    }
  }
  return build_set(child_l,child_r);
}

ptrdiff_t tree_set_store::create_set(std::vector<ptrdiff_t> &elems)
{
  if (elems.size() == 0)
  {
    return EMPTY_SET;
  }

  ptrdiff_t* nodes = MCRL2_SPECIFIC_STACK_ALLOCATOR(ptrdiff_t,elems.size());
  std::size_t node_size = 0;
  std::size_t i;
  std::size_t j;
  for (i=0; i < elems.size(); ++i)
  {
    nodes[i]=find_set(elems[i],EMPTY_SET);
  }
  node_size = i;
  while (node_size > 1)
  {
    j = 0;
    for (i=0; i<node_size; i+=2)
    {
      if (i == node_size-1)
      {
        nodes[j] = nodes[i];
      }
      else
      {
        nodes[j] = find_set(nodes[i],nodes[i+1]);
      }
      ++j;
    }
    node_size = j;
  }
  ptrdiff_t r = nodes[0];
  return r;
}

ptrdiff_t tree_set_store::get_next_tag()
{
  return tags_next;
}

ptrdiff_t tree_set_store::get_set(ptrdiff_t tag)
{
  return tags[tag];
}

ptrdiff_t tree_set_store::get_set_child_left(ptrdiff_t set)
{
  return buckets[set].child_l;
}

ptrdiff_t tree_set_store::get_set_child_right(ptrdiff_t set)
{
  return buckets[set].child_r;
}

ptrdiff_t tree_set_store::get_set_size(ptrdiff_t set)
{
  if (buckets[set].child_r == EMPTY_SET)
  {
    return 1;
  }
  return get_set_size(buckets[set].child_l) +
         get_set_size(buckets[set].child_r);
}

bool tree_set_store::is_set_empty(ptrdiff_t set)
{
  return (set == EMPTY_SET);
}

ptrdiff_t tree_set_store::set_set_tag(ptrdiff_t set)
{
  if (buckets[set].tag != EMPTY_TAG)
  {
    return buckets[set].tag;
  }
  check_tags();
  tags[tags_next] = set;
  buckets[set].tag = tags_next;
  return tags_next++;
}

}

