// Author(s): A.J. (Hannes) pretorius
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ./bundle.cpp

#include "bundle.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

// -- constructors and destructors ----------------------------------

Bundle::Bundle()
{
  index = NON_EXISTING;
  parent = 0;
  inCluster = 0;
  outCluster = 0;
}


Bundle::Bundle(const std::size_t& idx)
{
  index      = idx;
  parent     = 0;
  inCluster  = 0;
  outCluster = 0;
}


Bundle::Bundle(const Bundle& bdl)
{
  index      = bdl.index;
  parent     = bdl.parent;
  children   = bdl.children;
  inCluster  = bdl.inCluster;
  outCluster = bdl.outCluster;
  edges      = bdl.edges;
  labels     = bdl.labels;
}


Bundle::~Bundle()
{
  clearParent();
  clearChildren();
  clearInCluster();
  clearOutCluster();
  clearEdges();
}


// -- set functions -------------------------------------------------


void Bundle::addEdge(Edge* e)
{
  edges.push_back(e);

  labels.insert(std::pair< std::string, std::string >(e->getLabel(), "MAY"));
}


void Bundle::setEdges(const std::vector< Edge* > &e)
{
  clearEdges();
  edges = e;

  labels.clear();
  for (std::size_t i = 0; i < edges.size(); ++i)
  {
    labels.insert(std::pair< std::string, std::string >(edges[i]->getLabel(), "MAY"));
  }
}


// -- get functions -------------------------------------------------


Bundle* Bundle::getChild(const std::size_t& idx)
{
  if (idx >= children.size())
    return 0;
  return children[idx];
}


Edge* Bundle::getEdge(const std::size_t& idx)
{
  if (idx >= edges.size())
    throw mcrl2::runtime_error("Error retrieving bundle edge.");
  return edges[idx];
}


void Bundle::getLabels(std::vector< std::string > &lbls)
{
  lbls.clear();

  std::map< std::string, std::string >::iterator it;
  for (it = labels.begin(); it != labels.end(); ++it)
  {
    lbls.push_back(it->first);
  }

}


void Bundle::getLabels(
  std::vector< std::string > &lbls,
  std::vector< std::string > &status)
{
  lbls.clear();
  status.clear();

  std::map< std::string, std::string >::iterator it;
  for (it = labels.begin(); it != labels.end(); ++it)
  {
    lbls.push_back(it->first);
    status.push_back(it->second);
  }
}


void Bundle::getLabels(
  std::string& separator,
  std::string& lbls)
{
  lbls = "";

  std::size_t size  = labels.size();
  std::size_t count = 0;

  std::map< std::string, std::string >::iterator it;
  for (it = labels.begin(); it != labels.end(); ++it)
  {
    count += 1;
    lbls  += it->first;
    lbls  += "[";
    lbls  += it->second;
    lbls  += "]";

    if (count < size)
    {
      lbls += separator;
    }
  }
}


// -- clear functions -----------------------------------------------


void Bundle::clearChildren()
{
  for (std::size_t i = 0; i < children.size(); ++i)
  {
    children[i] = 0;
  }
  children.clear();
}


void Bundle::clearEdges()
{
  for (std::size_t i = 0; i < edges.size(); ++i)
  {
    edges[i] = 0;
  }
  edges.clear();

  labels.clear();
}


// -- end -----------------------------------------------------------
