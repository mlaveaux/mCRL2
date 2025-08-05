// Author(s): A.J. (Hannes) pretorius
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ./cluster.cpp

#include "cluster.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

// -- constructors and destructors ----------------------------------

Cluster::Cluster()
{
  index      = NON_EXISTING;
  parent     = 0;
  attribute  = 0;
  attrValIdx = NON_EXISTING;
}


Cluster::Cluster(const std::vector< std::size_t > &crd)
{
  coord      = crd;
  index      = NON_EXISTING;
  parent     = 0;
  attribute  = 0;
  attrValIdx = NON_EXISTING;
}


Cluster::Cluster(const Cluster& clst)
{
  coord      = clst.coord;
  index      = clst.index;
  parent     = clst.parent;
  children   = clst.children;
  nodes      = clst.nodes;
  attribute  = clst.attribute;
  attrValIdx = clst.attrValIdx;
  inBundles  = clst.inBundles;
  outBundles = clst.outBundles;
}


Cluster::~Cluster()
{
  clearParent();
  clearChildren();
  clearNodes();
  clearAttribute();
  clearInBundles();
  clearOutBundles();
}


// -- get functions -------------------------------------------------


std::size_t Cluster::getCoord(const std::size_t& idx)
{
  std::size_t result = NON_EXISTING;
  if (idx != NON_EXISTING && static_cast <std::size_t>(idx) < coord.size())
  {
    result = coord[idx];
  }
  return result;
}


Cluster* Cluster::getChild(const std::size_t& idx)
{
  if (idx != NON_EXISTING && static_cast <std::size_t>(idx) < children.size())
  {
    return children[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving cluster child.");
  }
}


Node* Cluster::getNode(const std::size_t& idx)
{
  if (idx != NON_EXISTING && static_cast <std::size_t>(idx) < nodes.size())
  {
    return nodes[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving cluster node.");
  }
}


std::size_t Cluster::getSizeDescNodes()
{
  std::size_t sum = 0;
  getSizeDescNodes(this, sum);
  return sum;
}


Bundle* Cluster::getInBundle(const std::size_t& idx)
{
  if (idx != NON_EXISTING && static_cast <std::size_t>(idx) < inBundles.size())
  {
    return inBundles[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving cluster incoming bundle.");
  }
}


Bundle* Cluster::getOutBundle(const std::size_t& idx)
{
  if (idx != NON_EXISTING && static_cast <std::size_t>(idx) < outBundles.size())
  {
    return outBundles[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving cluster outgoing bundle.");
  }
}


// -- clear functions -----------------------------------------------


void Cluster::clearChildren()
{
  for (std::size_t i = 0; i < children.size(); ++i)
  {
    children[i] = 0;
  }
  children.clear();
}


void Cluster::clearNodes()
{
  for (std::size_t i = 0; i < nodes.size(); ++i)
  {
    nodes[i] = 0;
  }
  nodes.clear();
}


void Cluster::clearAttribute()
{
  attribute  = 0;
  attrValIdx = NON_EXISTING;
}


void Cluster::clearInBundles()
{
  for (std::size_t i = 0; i < inBundles.size(); ++i)
  {
    inBundles[i] = 0;
  }
  inBundles.clear();
}


void Cluster::clearOutBundles()
{
  for (std::size_t i = 0; i < outBundles.size(); ++i)
  {
    outBundles[i] = 0;
  }
  outBundles.clear();
}


// -- utility functions ---------------------------------------------


void Cluster::getSizeDescNodes(
  Cluster* curClst,
  std::size_t& sum)
{
  if (curClst->getSizeChildren() == 0)
  {
    sum += curClst->getSizeNodes();
  }
  else
  {
    for (std::size_t i = 0; i < curClst->getSizeChildren(); ++i)
    {
      getSizeDescNodes(curClst->getChild(i), sum);
    }
  }
}



// -- end -----------------------------------------------------------
