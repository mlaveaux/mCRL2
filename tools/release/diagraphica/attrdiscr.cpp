// Author(s): A.J. (Hannes) pretorius
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ./attrdiscr.cpp

#include "attrdiscr.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

// -- constructors and destructor -----------------------------------

AttrDiscr::AttrDiscr(
  QString name,
  QString type,
  const std::size_t& idx,
  const std::vector< std::string >& vals)
  : Attribute(name, type, idx)
{
  initValues(vals);
}


AttrDiscr::AttrDiscr(const AttrDiscr& attr)
  :Attribute(attr)
{
  {
    for (std::size_t i = 0; i < attr.origValues.size(); ++i)
    {
      origValues.push_back(new Value(*(attr.origValues[i])));
    }
  }

  {
    for (std::size_t i = 0; i < attr.curValues.size(); ++i)
    {
      curValues.push_back(new Value(*(attr.curValues[i])));
    }
  }

  Value** mapping;
  {
    for (std::size_t i = 0; i < attr.curMap.size(); ++i)
    {
      mapping  = new Value*;
      *mapping = curValues[(*attr.curMap[i])->getIndex() ];
      curMap.push_back(mapping);
    }
  }
  mapping = 0;
}


AttrDiscr::~AttrDiscr()
{
  deleteOrigValues();
  deleteCurValues();
  deleteCurMap();
}


// -- set functions -------------------------------------------------


void AttrDiscr::clusterValues(
  const std::vector< int >& indices,
  const std::string& newValue)
{
  try
  {
    // variables
    std::vector< int > sorted;
    Value*        value;
    std::vector< std::vector< Value* >::iterator > toRemove;

    // sort indices in ascending order
    sorted = indices;
    sort(sorted.begin(), sorted.end());

    // ptr to new value
    value = new Value(NON_EXISTING, newValue);

    // update current map
    {
      for (std::size_t i = 0; i < curMap.size(); ++i)
      {
        for (std::size_t j = 0; j < sorted.size(); ++j)
        {
          if (*curMap[i] == curValues[sorted[j]])
          {
            *curMap[i] = value;
          }
        }
      }
    }

    // update current domain
    // get iterators to values to delete
    {
      for (std::size_t i = 0; i < sorted.size(); ++i)
      {
        toRemove.push_back(curValues.begin() + sorted[i]);
      }
    }

    // delete these values
    {
      for (std::size_t i = 0; i < toRemove.size(); ++i)
      {
        // delete object
        delete *(toRemove[i]);
      }
    }
    {
      for (std::size_t i = toRemove.size(); i > 0; --i)
      {
        // remove ptr from current domain
        curValues.erase(toRemove[i-1]);
      }
    }

    // update index
    value->setIndex(sorted[0]);
    // insert new value
    curValues.insert(curValues.begin() + sorted[0], value);

    // update value indices after new one
    {
      for (std::size_t i = sorted[0]; i < curValues.size(); ++i)
      {
        curValues[i]->setIndex(i);
      }
    }

    // -*-
    value = 0;
  }
  catch (...)
  {
    throw mcrl2::runtime_error("Error clustering attribute domain.");
  }

  emit changed();
}


void AttrDiscr::moveValue(
  const std::size_t& idxFr,
  const std::size_t& idxTo)
{
  try
  {
    Value* temp = curValues[idxFr];

    // 2 cases to consider
    if (idxFr < idxTo)
    {
      // move all values after idxFr 1 pos up
      for (std::size_t i = idxFr; i < idxTo; ++i)
      {
        curValues[i] = curValues[i+1];
        curValues[i]->setIndex(i);
      }
      // update idxTo
      curValues[idxTo] = temp;
      curValues[idxTo]->setIndex(idxTo);
    }
    else if (idxTo < idxFr)
    {
      // move all values before idxFr 1 pos down
      for (std::size_t i = idxFr; i > idxTo; --i)
      {
        curValues[i] = curValues[i-1];
        curValues[i]->setIndex(i);
      }
      // update idxTo
      curValues[idxTo] = temp;
      curValues[idxTo]->setIndex(idxTo);
    }

    temp = 0;
  }
  catch (...)
  {
    throw mcrl2::runtime_error("Error moving attribute domain value.");
  }

  emit changed();
}


void AttrDiscr::configValues(
  const std::vector< std::string >& curDomain,
  std::map< std::size_t , std::size_t  >& origToCurDomain)
{
  try
  {
    // clear current values
    {
      for (std::size_t i = 0; i < curValues.size(); ++i)
      {
        delete curValues[i];
        curValues[i] = 0;
      }
    }
    curValues.clear();

    // clear current map
    curMap.clear();

    // re-init current values
    {
      for (std::size_t i = 0; i < curDomain.size(); ++i)
      {
        curValues.push_back(new Value(i, curDomain[i]));
      }
    }

    // re-init mapping to current values
    Value** mapping;
    {
      for (std::size_t i = 0; i < origToCurDomain.size(); ++i)
      {
        mapping  = new Value*;
        *mapping = curValues[ origToCurDomain[i] ];
        curMap.push_back(mapping);
      }
    }
    mapping = 0;
  }
  catch (const mcrl2::runtime_error& e)
  {
    throw mcrl2::runtime_error(std::string("Error configuring attribute domain values.\n") + std::string(e.what()));
  }

  emit changed();
}


// -- get functions -------------------------------------------------


std::size_t AttrDiscr::getSizeOrigValues()
{
  return origValues.size();
}


Value* AttrDiscr::getOrigValue(std::size_t idx)
{
  if (idx != NON_EXISTING && idx < origValues.size())
  {
    return origValues[idx];
  }
  else
    throw mcrl2::runtime_error(
      "Error retrieving attribute value from original domain.");
}


std::size_t AttrDiscr::getSizeCurValues()
{
  return curValues.size();
}


Value* AttrDiscr::getCurValue(std::size_t idx)
{
  if (idx != NON_EXISTING && idx < curValues.size())
  {
    return curValues[idx];
  }
  else
    throw mcrl2::runtime_error(
      "Error retrieving attribute value from current domain.");
}


std::size_t AttrDiscr::getSizeMap()
{
  return curMap.size();
}


Value* AttrDiscr::mapToValue(double key)
{
  std::size_t idx = static_cast <std::size_t>(key);
  if (idx < curMap.size())
  {
    return *curMap[idx];
  }
  else
    throw mcrl2::runtime_error(
      "Error mapping index to current attribute value.");
}


// -- clear functions -----------------------------------------------


void AttrDiscr::clearClusters()
{
  resetCurValues();
  emit changed();
}


// -- private utility functions -------------------------------------


void AttrDiscr::initValues(const std::vector< std::string >& vals)
{
  Value*  value0  = 0;
  Value*  value1  = 0;
  Value** mapping = 0;

  // init orig & current domain, current map
  for (std::size_t i = 0; i < vals.size(); ++i)
  {
    // init new values
    value0 = new Value(i, vals[i]);
    value1 = new Value(i, vals[i]);

    // init domains
    origValues.push_back(value0);
    curValues.push_back(value1);

    // init mapping: pointer to pointer
    mapping  = new Value*;
    *mapping = value1;
    curMap.push_back(mapping);
  }

  value0  = 0;
  value1  = 0;
  mapping = 0;
}


void AttrDiscr::resetCurValues()
{
  Value* value    = 0;
  Value** mapping = 0;

  // clear current domain & mapping
  deleteCurValues();
  deleteCurMap();

  // reset current domain to original & update mapping
  for (std::size_t i = 0; i < origValues.size(); ++i)
  {
    // call copy constructor
    value = new Value(
      origValues[i]->getIndex(),
      origValues[i]->getValue());

    // init domain
    curValues.push_back(value);

    // init mapping: pointer to pointer
    mapping  = new Value*;
    *mapping = value;
    curMap.push_back(mapping);
  }

  value   = 0;
  mapping = 0;
}


void AttrDiscr::deleteOrigValues()
{
  for (std::size_t i = 0; i < origValues.size(); ++i)
  {
    delete origValues[i];
    origValues[i] = 0;
  }
  origValues.clear();
}


void AttrDiscr::deleteCurValues()
{
  for (std::size_t i = 0; i < curValues.size(); ++i)
  {
    delete curValues[i];
    curValues[i] = 0;
  }
  curValues.clear();
}


void AttrDiscr::deleteCurMap()
{
  for (std::size_t i = 0; i < curMap.size(); ++i)
  {
    delete curMap[i];
    curMap[i] = 0;
  }
  curMap.clear();
}


// -- end -----------------------------------------------------------
