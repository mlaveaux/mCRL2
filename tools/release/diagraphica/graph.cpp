// Author(s): A.J. (Hannes) pretorius
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ./graph.cpp

#include "graph.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

Graph::Graph()
{
  initRoot();
}


Graph::~Graph()
{
  deleteClusters();
  deleteBundles();
  deleteAttributes();
  deleteNodes();
  deleteEdges();
}


// -- set functions -------------------------------------------------


void Graph::setFileName(QString filename)
{
  m_filename = filename;
}

/*
void Graph::addAttribute(
    const std::string &name,
    const std::string &type,
    const int &idx,
    const std::vector< std::string > &vals )
{
    Attribute* attribute = new Attribute(
        name,
        type,
        idx,
        vals );
    attributes.push_back( attribute );
    attribute = 0;
}


void Graph::addAttribute(
    const std::string &name,
    const std::string &type,
    const int &idx,
    const double &lwrBnd,
    const double &uprBnd )
{
    Attribute* attribute = new Attribute(
        name,
        type,
        idx,
        lwrBnd,
        uprBnd );
    attributes.push_back( attribute );
    attribute = 0;
}
*/

void Graph::addAttrDiscr(
  QString name,
  QString type,
  const std::size_t& idx,
  const std::vector< std::string > &vals)
{
  AttrDiscr* attr = new AttrDiscr(
    name,
    type,
    idx,
    vals);
  attributes.push_back(attr);
  attr = 0;
}


void Graph::moveAttribute(
  const std::size_t& idxFr,
  const std::size_t& idxTo)
{
  if (idxFr < attributes.size() &&
      idxTo < attributes.size())
  {
    // move attributes
    Attribute* temp = attributes[idxFr];

    // 2 cases to consider
    if (idxFr < idxTo)
    {
      // move all values after idxFr 1 pos up
      for (std::size_t i = idxFr; i < idxTo; ++i)
      {
        attributes[i] = attributes[i+1];
        attributes[i]->setIndex(i);
      }
      // update idxTo
      attributes[idxTo] = temp;
      attributes[idxTo]->setIndex(idxTo);
    }
    else if (idxTo < idxFr)
    {
      // move all values before idxFr 1 pos down
      for (std::size_t i = idxFr; i > idxTo; --i)
      {
        attributes[i] = attributes[i-1];
        attributes[i]->setIndex(i);
      }
      // update idxTo
      attributes[idxTo] = temp;
      attributes[idxTo]->setIndex(idxTo);
    }

    temp = 0;

    // move node tuple values
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
      nodes[i]->moveTupleVal(idxFr, idxTo);
    }
  }
  else
  {
    throw mcrl2::runtime_error("Error moving attribute.");
  }

  attributes[idxTo]->emitMoved(idxTo);
}


void Graph::configAttributes(
  std::map< std::size_t, std::size_t > &idcsFrTo,
  std::map< std::size_t, std::vector< std::string > > &attrCurDomains,
  std::map< std::size_t, std::map< std::size_t , std::size_t  > > &attrOrigToCurDomains)
{
  std::size_t sizeAttrs = attributes.size();

  if (idcsFrTo.size()             == sizeAttrs &&
      attrCurDomains.size()       == sizeAttrs &&
      attrOrigToCurDomains.size() == sizeAttrs)
  {
    try
    {
      std::vector< Attribute* > attrsNew;

      // init new list of attributes
      {
        for (std::size_t i = 0; i < idcsFrTo.size(); ++i)
        {
          attrsNew.push_back(0);
        }
      }

      // update new list of attributes
      {
        for (std::size_t i = 0; i < idcsFrTo.size(); ++i)
        {
          attrsNew[ idcsFrTo[i] ] = attributes[i];
          attrsNew[ idcsFrTo[i] ]->setIndex(idcsFrTo[i]);
        }
      }

      // update attribute domains
      {
        for (std::size_t i = 0; i < idcsFrTo.size(); ++i)
          attrsNew[ idcsFrTo[i] ]->configValues(
            attrCurDomains[i],
            attrOrigToCurDomains[i]);
      }

      // set attributes to new list of attributes
      attributes.clear();
      attributes = attrsNew;
      attrsNew.clear();

      // update nodes
      {
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
          nodes[i]->moveTupleVals(idcsFrTo);
        }
      }
    }
    catch (const mcrl2::runtime_error& e)
    {
      throw mcrl2::runtime_error(std::string("Error moving attributes.\n") + std::string(e.what()));
    }

  }
  else
  {
    throw mcrl2::runtime_error("Error configuring attributes.");
  }

}


void Graph::duplAttributes(const std::vector< std::size_t > &idcs)
{
  std::size_t insIdx = 0;
  std::vector< Attribute* > newAttrs;

  // get insertion index
  {
    for (std::size_t i = 0; i < idcs.size(); ++i)
    {
      if (idcs[i] > insIdx)
      {
        insIdx = idcs[i];
      }
    }
  }
  insIdx++;

  // add attributes at insertion index
  {
    for (std::size_t i = 0; i < idcs.size(); ++i)
    {
      // add new attribute
      // -*-
      attributes.insert(
        attributes.begin() + insIdx + i,
        new AttrDiscr(*((AttrDiscr*)attributes[ idcs[i] ])));

      attributes[ insIdx + i ]->setIndex(insIdx + i);
      attributes[ insIdx + i ]->setName(
        "Copy_of_" + attributes[ idcs[i] ]->name());
      {
        for (std::size_t j = insIdx + i + 1; j < attributes.size(); ++j)
        {
          attributes[j]->setIndex(j);
        }
      }

      // update nodes
      {
        for (std::size_t j = 0; j < nodes.size(); ++j)
          nodes[j]->addTupleVal(
            insIdx + i,
            nodes[j]->getTupleVal(idcs[i]));
      }

      attributes[idcs[i]]->emitDuplicated();
    }
  }
}


void Graph::deleteAttribute(const std::size_t& idx)
{
  Cluster* clst;
  std::vector< std::size_t > idcsCurClust;
  std::vector< std::size_t > idcsNewClust;

  // get indcs of attributes used in current clustering
  clst = root;
  while (clst->getSizeChildren() > 0)
  {
    clst = clst->getChild(0);
    idcsCurClust.push_back(clst->getAttribute()->getIndex());
  }
  clst = 0;

  // get intersection of idcsToDelete & idcsCurClust
  {
    for (std::size_t i = 0; i < idcsCurClust.size(); ++i)
    {
      // if not found, insert into new clustering
      if (idx != idcsCurClust[i])
      {
        idcsNewClust.push_back(idcsCurClust[i]);
      }
    }
  }

  // init new clustering
  if (idcsNewClust.size() < idcsCurClust.size())
  {
    clustNodesOnAttr(idcsNewClust);
  }

  Attribute *attribute = attributes[idx];
  // update attributes & nodes
  attributes.erase(attributes.begin() + idx);
  {
    for (std::size_t i = idx; i < attributes.size(); ++i)
    {
      attributes[i]->setIndex(i);
    }
  }

  // update nodes
  {
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
      nodes[i]->delTupleVal(idx);
    }
  }

  attribute->emitDeleted();
  delete attribute;
  emit deletedAttribute();
}


void Graph::addNode(const std::vector< double > &tpl)
{
  Node* n = new Node(
    nodes.size(),
    tpl);
  addNode(n);
  n = 0;
}


void Graph::addEdge(
  const std::string& lbl,
  const std::size_t& inNodeIdx,
  const std::size_t& outNodeIdx)
{
  Edge* e = new Edge(
    edges.size(),
    lbl);

  e->setInNode(nodes[inNodeIdx]);
  e->setOutNode(nodes[outNodeIdx]);
  nodes[inNodeIdx]->addOutEdge(e);
  nodes[outNodeIdx]->addInEdge(e);

  Bundle* b = getBundle(0);
  b->addEdge(e);
  e->setBundle(b);
  b = 0;

  addEdge(e);
  e = 0;
}


void Graph::initGraph()
// This function initializes the graph and should be called after all
// nodes and edges have been added.
{
  //initRoot();
  updateLeaves();
  updateBundles();
}


// -- get functions -------------------------------------------------


QString Graph::filename()
{
  return m_filename;
}


std::size_t Graph::getSizeAttributes()
{
  return attributes.size();
}


Attribute* Graph::getAttribute(const std::size_t& idx)
{
  assert(idx < attributes.size());
  return attributes[idx];
}


Attribute* Graph::getAttribute(QString name)
{
  Attribute* result = 0;

  for (std::size_t i = 0; i < attributes.size() && result == 0; ++i)
  {
    if (attributes[i]->name() == name)
    {
      result = attributes[i];
    }
  }

  return result;
}


std::size_t Graph::getSizeNodes()
{
  return nodes.size();
}


Node* Graph::getNode(const std::size_t& idx)
{
  if (idx < nodes.size())
  {
    return nodes[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving node.");
  }
}


std::size_t Graph::getSizeEdges()
{
  return edges.size();
}


Edge* Graph::getEdge(const std::size_t& idx)
{
  if (idx < edges.size())
  {
    return edges[idx];
  }
  else
  {
    throw mcrl2::runtime_error("Error retrieving edge.");
  }
}


Cluster* Graph::getRoot()
{
  return root;
}


Cluster* Graph::getCluster(const std::vector< std::size_t > coord)
{
  Cluster* result = 0;

  if (coord.size() > 1)
  {
    Cluster* temp = root;

    for (std::size_t i = 1; i < coord.size(); ++i)
    {
      if (coord[i] != NON_EXISTING  && coord[i] < temp->getSizeChildren())
      {
        temp = temp->getChild(coord[i]);
      }
      else
      {
        temp = 0;
        break;
      }
    }

    result = temp;
    temp   = 0;
  }
  else if (coord.size() == 1 && coord[0] == 0)
  {
    result = root;
  }

  return result;
}


Cluster* Graph::getLeaf(const std::size_t& idx)
{
  Cluster* result = 0;
  if (idx < leaves.size())
  {
    result = leaves[idx];
  }
  return result;
}


std::size_t Graph::getSizeLeaves()
{
  return leaves.size();
}


Bundle* Graph::getBundle(const std::size_t& idx)
{
  Bundle* result = 0;
  if (idx < bundles.size())
  {
    result = bundles[idx];
  }
  return result;
}


std::size_t Graph::getSizeBundles()
{
  return bundles.size();
}


// -- calculation functions -----------------------------------------


void Graph::calcAttrDistr(
  const std::size_t& attrIdx,
  std::vector< std::size_t > &distr)
{
  // vars
  Attribute* attribute = 0;
  std::size_t sizeDomain       = 0;
  std::size_t sizeNodes        = 0;

  // init vars
  attribute  = getAttribute(attrIdx);
  sizeDomain = attribute->getSizeCurValues();
  sizeNodes  = nodes.size();

  // init results
  distr.clear();
  if (sizeDomain > 0)
  {
    for (std::size_t i = 0; i < sizeDomain; ++i)
    {
      distr.push_back(0);
    }
  }

  // calc results
  if (sizeDomain > 0)
  {
    for (std::size_t i = 0; i < sizeNodes; ++i)
    {
      distr[ attribute->mapToValue(nodes[i]->getTupleVal(attrIdx))->getIndex() ] += 1;
    }
  }

  // reset ptr
  attribute = 0;
}


void Graph::calcAttrCorrl(
  const std::size_t& attrIdx1,
  const std::size_t& attrIdx2,
  std::vector< std::vector< std::size_t > > &corrlMap,
  std::vector< std::vector< int > > &number)
{
  // vars
  Attribute* attr1 = 0;
  Attribute* attr2 = 0;
  Node*      node  = 0;
  std::size_t sizeDomain1  = 0;
  std::size_t sizeDomain2  = 0;
  std::size_t sizeNodes    = 0;
  std::size_t domIdx1      = 0;
  std::size_t domIdx2      = 0;
  std::vector< std::size_t > toErase;

  // init vars
  attr1       = getAttribute(attrIdx1);
  attr2       = getAttribute(attrIdx2);
  sizeDomain1 = attr1->getSizeCurValues();
  sizeDomain2 = attr2->getSizeCurValues();
  sizeNodes   = nodes.size();

  // init results;
  corrlMap.clear();
  number.clear();
  if (sizeDomain1 > 0 && sizeDomain2 > 0)
  {
    for (std::size_t i = 0; i < sizeDomain1; ++i)
    {
      std::vector< int > tempNumVec;
      for (std::size_t j = 0; j < sizeDomain2; ++j)
      {
        tempNumVec.push_back(0);
      }
      number.push_back(tempNumVec);

      std::vector< std::size_t > tempMapVec;
      corrlMap.push_back(tempMapVec);
    }
  }

  // calc prelim results
  if (sizeDomain1 > 0 && sizeDomain2 > 0)
  {
    for (std::size_t i = 0; i < sizeNodes; ++i)
    {
      node    = nodes[i];
      /*
      domIdx1 = attr1->mapToValue(
          (int)node->getTupleVal(
              attrIdx1 ) )->getIndex();
      domIdx2 = attr2->mapToValue(
          (int)node->getTupleVal(
              attrIdx2 ) )->getIndex();
      */
      domIdx1 = attr1->mapToValue(
                  node->getTupleVal(
                    attrIdx1))->getIndex();
      domIdx2 = attr2->mapToValue(
                  node->getTupleVal(
                    attrIdx2))->getIndex();

      number[domIdx1][domIdx2] += 1;
    }
  }

  // update correlation map
  if (sizeDomain1 > 0 && sizeDomain2 > 0)
  {
    for (std::size_t i = 0; i < sizeDomain1; ++i)
    {
      for (std::size_t j = 0; j < sizeDomain2; ++j)
      {
        if (number[i][j] > 0)
        {
          corrlMap[i].push_back(j);
        }
      }
    }
  }

  // remove zero entries from number
  if (sizeDomain1 > 0 && sizeDomain2 > 0)
  {
    for (std::size_t i = 0; i < sizeDomain1; ++i)
    {
      toErase.clear();
      {
        for (std::size_t j = 0; j < sizeDomain2; ++j)
        {
          if (number[i][j] < 1)
          {
            toErase.push_back(j);
          }
        }
      }

      {
        for (std::size_t j = 0; j < toErase.size(); ++j)
        {
          // brackets around (toErase[j] - j) REALLY important
          number[i].erase(number[i].begin() + (toErase[j] - j));
        }
      }
    }
  }

  // reset ptrs
  attr1 = 0;
  attr2 = 0;
  node  = 0;
}


void Graph::calcAttrCombn(
  const std::vector< std::size_t > &attrIndcs,
  std::vector< std::vector< std::size_t > > &combs,
  std::vector< std::size_t > &number)
{
  if (root != 0)
  {
    calcAttrCombn(root, attrIndcs, combs, number);
  }
}


void Graph::calcAttrCombn(
  Cluster* clust,
  const std::vector< std::size_t > &attrIndcs,
  std::vector< std::vector< std::size_t > > &combs,
  std::vector< std::size_t > &number)
{
  std::size_t summand     = 0;
  std::size_t key         = 0;
  std::size_t card        = 0;
  Node*      node = 0;
  Attribute* attr = 0;
  std::vector< Node* >           nodesInClst;
  std::map< std::size_t , std::vector< std::size_t > > keyToCombn;
  std::map< std::size_t , std::size_t >           keyToNumber;

  // do initialization
  getDescNodesInCluster(clust, nodesInClst);

  // calc results
  for (std::size_t i = 0; i < nodesInClst.size(); ++i)
  {
    key  = 0;
    node = nodesInClst[i];

    // calc key
    for (std::size_t j = 0; j < attrIndcs.size(); ++j)
    {
      attr = getAttribute(attrIndcs[j]);
      card = attr->getSizeCurValues();
      if (card > 0)
        /*
        summand = attr->mapToValue(
            (int)node->getTupleVal( attrIndcs[j] ) )->getIndex();
        */
        summand = attr->mapToValue(
                    node->getTupleVal(attrIndcs[j]))->getIndex();
      else
      {
        summand = 0;
      }

      for (std::size_t k = j+1; k < attrIndcs.size(); ++k)
      {
        attr = getAttribute(attrIndcs[k]);
        card = attr->getSizeCurValues();

        if (card > 0)
        {
          summand *= card;
        }
      }

      key += summand;
    } // for j

    // insert combn 1st time it occurs
    if (keyToCombn.find(key) == keyToCombn.end())
    {
      std::vector< std::size_t > v;
      for (std::size_t j = 0; j < attrIndcs.size(); ++j)
      {
        attr = getAttribute(attrIndcs[j]);
        card = attr->getSizeCurValues();

        if (card > 0)
          /*
          v.push_back( attr->mapToValue(
              (int)node->getTupleVal( attrIndcs[j] ) )->getIndex() );
          */
          v.push_back(attr->mapToValue(
                        node->getTupleVal(attrIndcs[j]))->getIndex());
      }

      keyToCombn.insert(std::pair< std::size_t, std::vector< std::size_t > >(key, v));
      keyToNumber.insert(std::pair< std::size_t, std::size_t >(key, 1));
    }
    // increment number if combn exists
    else
    {
      keyToNumber[key] += 1;
    }

  } // for i

  // update combs
  combs.clear();

  std::map< std::size_t , std::vector< std::size_t > >::iterator it;
  for (it = keyToCombn.begin(); it != keyToCombn.end(); ++it)
  {
    combs.push_back((*it).second);
    number.push_back(keyToNumber[(*it).first]);
  }

  // clear memory
  node = 0;
  attr = 0;
  nodesInClst.clear();
}


void Graph::calcAttrCombn(
  const std::vector< std::size_t > &attrIndcs,
  std::vector< std::vector< std::size_t > > &combs)
{
  if (root != 0)
  {
    calcAttrCombn(root, attrIndcs, combs);
  }
}



void Graph::calcAttrCombn(
  Cluster* clust,
  const std::vector< std::size_t > &attrIndcs,
  std::vector< std::vector< std::size_t > > &combs)
{
  std::size_t summand     = 0;
  std::size_t key         = 0;
  std::size_t card        = 0;
  Node*      node = 0;
  Attribute* attr = 0;
  std::vector< Node* >           nodesInClst;
  std::map< std::size_t , std::vector< std::size_t > > keyToCombn;

  // do initialization
  getDescNodesInCluster(clust, nodesInClst);

  // calc results
  for (std::size_t i = 0; i < nodesInClst.size(); ++i)
  {
    key  = 0;
    node = nodesInClst[i];

    // calc key
    for (std::size_t j = 0; j < attrIndcs.size(); ++j)
    {
      attr = getAttribute(attrIndcs[j]);
      card = attr->getSizeCurValues();
      if (card > 0)
        /*
        summand = attr->mapToValue(
            (int)node->getTupleVal( attrIndcs[j] ) )->getIndex();
        */
        summand = attr->mapToValue(
                    node->getTupleVal(attrIndcs[j]))->getIndex();
      else
      {
        summand = 0;
      }

      for (std::size_t k = j+1; k < attrIndcs.size(); ++k)
      {
        attr = getAttribute(attrIndcs[k]);
        card = attr->getSizeCurValues();

        if (card > 0)
        {
          summand *= card;
        }
      }

      key += summand;
    } // for j

    // insert combn 1st time it occurs
    if (keyToCombn.find(key) == keyToCombn.end())
    {
      std::vector< std::size_t > v;
      for (std::size_t j = 0; j < attrIndcs.size(); ++j)
      {
        attr = getAttribute(attrIndcs[j]);
        card = attr->getSizeCurValues();

        if (card > 0)
          /*
          v.push_back( attr->mapToValue(
              (int)node->getTupleVal( attrIndcs[j] ) )->getIndex() );
          */
          v.push_back(attr->mapToValue(
                        node->getTupleVal(attrIndcs[j]))->getIndex());
      }

      keyToCombn.insert(std::pair< std::size_t, std::vector< std::size_t > >(
                          key,
                          v));
    }
  } // for i

  // update combs
  combs.clear();

  std::map< std::size_t, std::vector< std::size_t > >::iterator it;
  for (it = keyToCombn.begin(); it != keyToCombn.end(); ++it)
  {
    combs.push_back((*it).second);
  }

  // clear memory
  node = 0;
  attr = 0;
  nodesInClst.clear();
}


void Graph::calcAttrCombn(
  Cluster* clust,
  const std::vector< std::size_t > &attrIndcs,
  std::vector< std::vector< Node* > > &combs)
{
  std::size_t summand     = 0;
  std::size_t key         = 0;
  std::size_t card        = 0;
  Node*      node = 0;
  Attribute* attr = 0;
  std::vector< Node* >           nodesInClst;
  std::map< std::size_t , std::vector< Node* > > keyToCombn;

  // do initialization
  getDescNodesInCluster(clust, nodesInClst);

  // calc results
  for (std::size_t i = 0; i < nodesInClst.size(); ++i)
  {
    key  = 0;
    node = nodesInClst[i];

    // calc key
    for (std::size_t j = 0; j < attrIndcs.size(); ++j)
    {
      attr = getAttribute(attrIndcs[j]);
      card = attr->getSizeCurValues();
      if (card > 0)
        /*
        summand = attr->mapToValue(
            (int)node->getTupleVal( attrIndcs[j] ) )->getIndex();
        */
        summand = attr->mapToValue(
                    node->getTupleVal(attrIndcs[j]))->getIndex();
      else
      {
        summand = 0;
      }

      for (std::size_t k = j+1; k < attrIndcs.size(); ++k)
      {
        attr = getAttribute(attrIndcs[k]);
        card = attr->getSizeCurValues();

        if (card > 0)
        {
          summand *= card;
        }
      }

      key += summand;
    } // for j

    std::map< std::size_t, std::vector< Node* > >::iterator pos;
    pos = keyToCombn.find(key);
    // insert combn 1st time it occurs
    if (pos == keyToCombn.end())
    {
      std::vector< Node* > v;
      v.push_back(node);
      keyToCombn.insert(std::pair< std::size_t , std::vector< Node* > >(key, v));
    }
    // insert state if it occurs again
    else
    {
      pos->second.push_back(node);
    }
  } // for i

  // update combs
  combs.clear();
  std::map< std::size_t, std::vector< Node* > >::iterator it;
  for (it = keyToCombn.begin(); it != keyToCombn.end(); ++it)
  {
    combs.push_back((*it).second);
  }

  // clear memory
  node = 0;
  attr = 0;
  nodesInClst.clear();
}


void Graph::calcAttrCombn(
  Cluster* clust,
  const std::vector< Attribute* > &attrs,
  std::vector< Cluster* > &combs)
{
  std::size_t summand     = 0;
  std::size_t key         = 0;
  std::size_t card        = 0;
  Node*      node = 0;
  Attribute* attr = 0;
  std::vector< Node* >      nodesInClst;
  std::map< std::size_t, Cluster* > keyToCombn;

  // do initialization
  getDescNodesInCluster(clust, nodesInClst);

  // calc results
  for (std::size_t i = 0; i < nodesInClst.size(); ++i)
  {
    key  = 0;
    node = nodesInClst[i];

    // calc key
    for (std::size_t j = 0; j < attrs.size(); ++j)
    {
      attr = attrs[j];
      card = attr->getSizeCurValues();
      if (card > 0)
        /*
        summand = attr->mapToValue(
            (int)node->getTupleVal( attr->getIndex() ) )->getIndex();
        */
        summand = attr->mapToValue(
                    node->getTupleVal(attr->getIndex()))->getIndex();
      else
      {
        summand = 0;
      }

      for (std::size_t k = j+1; k < attrs.size(); ++k)
      {
        attr = attrs[k];
        card = attr->getSizeCurValues();

        if (card > 0)
        {
          summand *= card;
        }
      }

      key += summand;
    } // for j

    std::map< std::size_t, Cluster* >::iterator pos;
    pos = keyToCombn.find(key);
    // insert combn 1st time it occurs
    if (pos == keyToCombn.end())
    {
      Cluster* clst = new Cluster();
      clst->addNode(node);
      keyToCombn.insert(std::pair< std::size_t, Cluster* >(key, clst));
      clst = 0;
    }
    // insert state if it occurs again
    else
    {
      pos->second->addNode(node);
    }
  } // for i

  // update combs
  combs.clear();
  std::map< std::size_t, Cluster* >::iterator it;
  for (it = keyToCombn.begin(); it != keyToCombn.end(); ++it)
  {
    combs.push_back(it->second);
  }

  // clear memory
  node = 0;
  attr = 0;
  nodesInClst.clear();
}


bool Graph::hasMultAttrCombns(
  Cluster* clust,
  const std::vector< int > &attrIndcs)
// This function returns true as soon as more than 1 state is found
// with different combinations of values for the attributes specified.
{
  bool result = false;
  std::size_t  numAttrs  = 0;
  std::size_t  cardAttr  = 0;
  std::size_t  sizeNodes = 0;
  std::size_t  mapFirst  = 0;
  std::size_t  map       = 0;
  std::size_t  summand   = 0;
  Node*      node      = 0;
  Attribute* attribute = 0;
  std::vector< Node* > clustNodes;

  getDescNodesInCluster(clust, clustNodes);
  numAttrs  = attrIndcs.size();
  sizeNodes = clustNodes.size();

  // calc results
  {
    for (std::size_t i = 0; i < sizeNodes && result == false; ++i)
    {
      map  = 0;
      node = clustNodes[i];

      // calc map
      for (std::size_t j = 0; j < numAttrs && result == false; ++j)
      {
        attribute = getAttribute(attrIndcs[j]);
        cardAttr  = attribute->getSizeCurValues();
        if (cardAttr > 0)
        {
          /*
          summand = attribute->mapToValue(
              (int)node->getTupleVal( attrIndcs[j] ) )->getIndex();
          */
          summand = attribute->mapToValue(
                      node->getTupleVal(attrIndcs[j]))->getIndex();
        }
        else
        {
          summand = 0;
        }

        for (std::size_t k = j+1; k < numAttrs; ++k)
        {
          attribute = getAttribute(attrIndcs[k]);
          cardAttr  = attribute->getSizeCurValues();

          if (cardAttr > 0)
          {
            summand *= cardAttr;
          }
        }

        map += summand;
      }

      if (i == 0)
      {
        mapFirst = map;
      }
      else if (mapFirst != map)
      {
        result = true;
      }
    }
  }

  node = 0;
  clustNodes.clear();

  return result;
}


// -- cluster functions ---------------------------------------------


void Graph::clustNodesOnAttr(const std::vector< std::size_t > &attrIdcs)
{
  std::size_t progress       = 0;
  std::vector< std::size_t > idcs = attrIdcs;

  disconnect(this, SLOT(recluster()));
  for (std::size_t i = 0; i < attrIdcs.size(); ++i)
  {
    connect(getAttribute(attrIdcs[i]), SIGNAL(changed()), this, SLOT(recluster()));
  }

  // cluster nodes
  clustNodesOnAttr(root, idcs, progress);

  // update leaves
  updateLeaves();

  // update bundles
  progress     = 0;
  updateBundles(progress);

  idcs.clear();
  emit clusteringChanged();
}


void Graph::clearSubClusters(const std::vector< std::size_t > &coord)
{
  Cluster* clst = 0;

  // get cluster
  clst = getCluster(coord);
  // clear subclusters
  clearSubClusters(clst);
  // update leaves
  updateLeaves();
  // update bundles
  updateBundles();

  clst = 0;
}


std::size_t Graph::sumNodesInCluster(const std::vector< std::size_t > &coord)
{
  std::size_t      sum  = 0;
  Cluster* clst = 0;

  // get cluster
  clst = getCluster(coord);
  // get sum
  if (clst != 0)
  {
    sumNodesInCluster(clst, sum);
  }

  clst = 0;
  return sum;
}


void Graph::sumNodesInCluster(
  Cluster* clust,
  std::size_t& total)
{
  for (std::size_t i = 0; i < clust->getSizeChildren(); ++i)
  {
    sumNodesInCluster(
      clust->getChild(i),
      total);
  }
  total += clust->getSizeNodes();
}


void Graph::getDescNodesInCluster(
  const std::vector< std::size_t > &coord,
  std::vector< Node* > &nodes)
{
  Cluster* clst = 0;

  // get cluster
  clst = getCluster(coord);

  // get descendant nodes
  nodes.clear();
  if (clst != 0)
  {
    getDescNodesInCluster(clst, nodes);
  }

  clst = 0;
}


void Graph::getDescNodesInCluster(
  Cluster* clust,
  std::vector< Node* > &nodes)
{
  // call recursively on all child clusters
  {
    for (std::size_t i = 0; i < clust->getSizeChildren(); ++i)
    {
      getDescNodesInCluster(
        clust->getChild(i),
        nodes);
    }
  }

  // add nodes
  {
    for (std::size_t i = 0; i < clust->getSizeNodes(); ++i)
    {
      nodes.push_back(clust->getNode(i));
    }
  }
}


std::size_t Graph::calcMaxNumCombns(const std::vector< std::size_t > &attrIdcs)
{
  std::size_t combinations     = 1;
  Attribute* attribute = 0;
  std::size_t cardinality      = 0;

  for (std::size_t i = 0; i < attrIdcs.size(); ++i)
  {
    attribute   = getAttribute(attrIdcs[i]);
    cardinality = attribute->getSizeCurValues();
    if (cardinality > 0)
    {
      combinations *= cardinality;
    }
  }

  attribute = 0;
  return combinations;
}


void Graph::recluster()
{
  std::vector<std::size_t> attributes;
  for(Cluster *cluster = getLeaf(0); cluster != getRoot(); cluster = cluster->getParent())
  {
    attributes.insert(attributes.begin(), cluster->getAttribute()->getIndex());
  }
  clustNodesOnAttr(attributes);
}


// -- private utility functions -------------------------------------


void Graph::deleteAttributes()
{
  for (std::size_t i = 0; i < attributes.size(); ++i)
  {
    delete attributes[i];
    attributes[i] = 0;
  }
  attributes.clear();
}


void Graph::addNode(Node* n)
{
  nodes.push_back(n);

  root->addNode(n);
  n->setCluster(root);
}


void Graph::deleteNodes()
{
  for (std::size_t i = 0; i < nodes.size(); ++i)
  {
    delete nodes[i];
    nodes[i] = 0;
  }
  nodes.clear();
}


void Graph::addEdge(Edge* e)
{
  edges.push_back(e);
}


void Graph::deleteEdges()
{
  for (std::size_t i = 0; i < edges.size(); ++i)
  {
    delete edges[i];
    edges[i] = 0;
  }
  edges.clear();
}


void Graph::initRoot()
{
  std::vector< std::size_t > rootCoord;
  rootCoord.push_back(0);

  root = new Cluster(rootCoord);

  Bundle* bndl = new Bundle(0);
  bundles.push_back(bndl);
  bndl = 0;
}


void Graph::clustNodesOnAttr(
  Cluster* clust,
  std::vector< std::size_t > attrIdcs,
  std::size_t& progress)
{
  if (attrIdcs.size() > 0)
  {
    // cluster
    clustClusterOnAttr(
      clust,
      attrIdcs[0]);

    // update progress
    if (attrIdcs.size() == 1)
    {
      progress += getAttribute(attrIdcs[0])->getSizeCurValues();
    }

    // remove first attribute
    attrIdcs.erase(attrIdcs.begin());

    // repeat recursively on children
    for (std::size_t i = 0; i < clust->getSizeChildren(); ++i)
    {
      clustNodesOnAttr(
        clust->getChild(i),
        attrIdcs,
        progress);
    }
  }
}


void Graph::clustClusterOnAttr(
  const std::vector< std::size_t > coord,
  const std::size_t& attrIdx)
{
  Cluster* clst = 0;

  // get cluster
  clst = getCluster(coord);
  // do clustering
  clustClusterOnAttr(
    clst,
    attrIdx);
  // update leaves
  updateLeaves();
  // update bundles
  updateBundles();

  clst = 0;
}


void Graph::clustClusterOnAttr(
  Cluster* clust,
  const std::size_t& attrIdx)
{
  Attribute*         attr;
  Node*              node;
  std::vector< Cluster* > clstTmp;
  std::vector< std::size_t >      clstCoordNew;

  // delete descendants & move up nodes
  clearSubClusters(clust);

  // get attribute
  attr = getAttribute(attrIdx);

  // init temporary clusters
  {
    for (std::size_t i = 0; i < attr->getSizeCurValues(); ++i)
    {
      clstTmp.push_back(new Cluster());
    }
  }

  // do clustering & move nodes to new cluster
  {
    for (std::size_t i = 0; i < clust->getSizeNodes(); ++i)
    {
      node = clust->getNode(i);

      /*
      clstIdxTmp = attr->mapToValue(
          (int)node->getTupleVal(attrIdx) )->getIndex();
      */
      std::size_t clstIdxTmp = attr->mapToValue(
                            node->getTupleVal(attrIdx))->getIndex();

      clstTmp[clstIdxTmp]->addNode(node);
      node->setCluster(clstTmp[clstIdxTmp]);
    }
  }
  clust->clearNodes();

  // update clust
  {
    for (std::size_t i = 0; i < clstTmp.size(); ++i)
    {
      if (clstTmp[i]->getSizeNodes() > 0)
      {
        clust->getCoord(clstCoordNew);
        clstCoordNew.push_back(clust->getSizeChildren());

        clstTmp[i]->setCoord(clstCoordNew);
        clust->addChild(clstTmp[i]);
        clstTmp[i]->setParent(clust);

        clstTmp[i]->setAttribute(attr);
        clstTmp[i]->setAttrValIdx(i);
      }
      else
      {
        delete clstTmp[i];
        clstTmp[i] = 0;
      }
    }
  }

  // free memory
  attr = 0;
  node = 0;
  clstTmp.clear();
}


void Graph::clearSubClusters(Cluster* clust)
{
  if (clust != 0)
  {
    for (std::size_t i = 0; i < clust->getSizeChildren(); ++i)
    {
      Cluster* child = clust->getChild(i);

      // clear child
      clearSubClusters(child);

      // move up nodes
      for (std::size_t j = 0; j < child->getSizeNodes(); ++j)
      {
        clust->addNode(child->getNode(j));
        child->getNode(j)->setCluster(clust);
      }

      delete child;
      child = 0;
    }

    clust->clearChildren();
  }
}


void Graph::updateLeaves()
{
  // clear existing leaves
  clearLeaves();

  // init new leaves
  /*
  std::vector< Cluster* > v;
  v.push_back( root );
  updateLeaves( v );
  */
  updateLeaves(root);
}


void Graph::updateLeaves(Cluster* clust)
{
  if (clust->getSizeChildren() == 0)
  {
    // set index
    clust->setIndex(leaves.size());

    // update leaves
    leaves.push_back(clust);
  }
  else
  {
    // reset index
    clust->setIndex(NON_EXISTING);

    // update clusts
    for (std::size_t i = 0; i < clust->getSizeChildren(); ++i)
    {
      updateLeaves(clust->getChild(i));
    }
  }
}

/*
void Graph::updateLeaves( std::vector< Cluster* > &clusts )
{
    if ( clusts.size() > 0 )
    {
        Cluster* c = clusts[0];

        if ( c->getSizeChildren() > 0 )
        {
            // reset index
            c->setIndex( -1 );

            // update clusts
            for ( int i = 0; i < c->getSizeChildren(); ++i )
                clusts.push_back( c->getChild( i ) );
        }
        else
        {
            // set index
            c->setIndex( leaves.size() );

            // update leaves
            leaves.push_back( c );

            std::vector< int > coord;
            c->getCoord( coord );
        }

        c = 0;

        clusts.erase( clusts.begin() );
        updateLeaves( clusts );
    }
}
*/

void Graph::clearLeaves()
{
  for (std::size_t i = 0; i < leaves.size(); ++i)
  {
    leaves[i] = 0;
  }
  leaves.clear();
}


void Graph::deleteClusters()
{
  // delete all descendants
  clearSubClusters(root);

  // delete root
  delete root;
  root = 0;

  // clean up leaves
  clearLeaves();
}


void Graph::updateBundles(std::size_t& progress)
{
  std::vector< std::vector< Bundle* > > temp;

  // clear bundles
  deleteBundles();

  // init temp data struct
  {
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
      std::vector< Bundle* > v;
      temp.push_back(v);

      for (std::size_t j = 0; j < leaves.size(); ++j)
      {
        temp[i].push_back(0);
      }
    }
  }

  // iterate over edges & update temp data struct
  {
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
      std::size_t idxFr = edges[i]->getInNode()->getCluster()->getIndex();
      std::size_t idxTo = edges[i]->getOutNode()->getCluster()->getIndex();

      if (temp[idxFr][idxTo] == 0)
      {
        temp[idxFr][idxTo] = new Bundle();

        // update bundle
        temp[idxFr][idxTo]->setInCluster(edges[i]->getInNode()->getCluster());
        temp[idxFr][idxTo]->setOutCluster(edges[i]->getOutNode()->getCluster());

        // update clusters
        edges[i]->getInNode()->getCluster()->addOutBundle(temp[idxFr][idxTo]);
        edges[i]->getOutNode()->getCluster()->addInBundle(temp[idxFr][idxTo]);
      }

      temp[idxFr][idxTo]->addEdge(edges[i]);
      edges[i]->setBundle(temp[idxFr][idxTo]);

      // update progress
      ++progress;
    }
  }

  // iterate over temp data struct & update bundles
  int idxBdl = 0;
  {
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
      for (std::size_t j = 0; j < leaves.size(); ++j)
      {
        if (temp[i][j] != 0)
        {
          temp[i][j]->setIndex(idxBdl);
          bundles.push_back(temp[i][j]);

          idxBdl++;
        }
      }
    }
  }
  temp.clear();

  // iterate over bundles & update may or must
  std::vector< std::string > labels;
  Cluster* inClust;
  Cluster* outClust;
  Node*    node;
  {
    for (std::size_t i = 0; i < bundles.size(); ++i)
    {
      labels.clear();
      bundles[i]->getLabels(labels);
      inClust  = bundles[i]->getInCluster();
      outClust = bundles[i]->getOutCluster();

      for (std::size_t j = 0; j < labels.size(); ++j)
      {
        bool must = true;
        for (std::size_t k = 0; k < inClust->getSizeNodes() && must; ++k)
        {
          node = inClust->getNode(k);
          bool hasLbl = false;

          for (std::size_t m = 0; m < node->getSizeOutEdges() && !hasLbl; ++m)
          {
            if (node->getOutEdge(m)->getLabel() == labels[j] &&
                node->getOutEdge(m)->getOutNode()->getCluster() == outClust)
            {
              hasLbl = true;
            }
          }

          if (!hasLbl)
          {
            must = false;
          }
        }

        if (must)
        {
          bundles[i]->updateLabel(labels[j], "MUST");
        }
        /*
        // by default labels are may
        else
            bundles[i]->updateLabel( labels[j], "MAY" );
        */
      }
    }
  }

  inClust  = 0;
  outClust = 0;
  node     = 0;
}


void Graph::updateBundles()
{
  std::vector< std::vector< Bundle* > > temp;

  // clear bundles
  deleteBundles();

  // init temp data struct
  {
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
      std::vector< Bundle* > v;
      temp.push_back(v);

      for (std::size_t j = 0; j < leaves.size(); ++j)
      {
        temp[i].push_back(0);
      }
    }
  }

  // iterate over edges & update temp data struct
  {
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
      std::size_t idxFr = edges[i]->getInNode()->getCluster()->getIndex();
      std::size_t idxTo = edges[i]->getOutNode()->getCluster()->getIndex();

      if (temp[idxFr][idxTo] == 0)
      {
        temp[idxFr][idxTo] = new Bundle();

        // update bundle
        temp[idxFr][idxTo]->setInCluster(edges[i]->getInNode()->getCluster());
        temp[idxFr][idxTo]->setOutCluster(edges[i]->getOutNode()->getCluster());

        // update clusters
        edges[i]->getInNode()->getCluster()->addOutBundle(temp[idxFr][idxTo]);
        edges[i]->getOutNode()->getCluster()->addInBundle(temp[idxFr][idxTo]);
      }

      temp[idxFr][idxTo]->addEdge(edges[i]);
      edges[i]->setBundle(temp[idxFr][idxTo]);
    }
  }

  // iterate over temp data struct & update bundles
  int idxBdl = 0;
  {
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
      for (std::size_t j = 0; j < leaves.size(); ++j)
      {
        if (temp[i][j] != 0)
        {
          temp[i][j]->setIndex(idxBdl);
          bundles.push_back(temp[i][j]);

          idxBdl++;
        }
      }
    }
  }

  temp.clear();
}


void Graph::deleteBundles()
{
  for (std::size_t i = 0; i < bundles.size(); ++i)
  {
    delete bundles[i];
    bundles[i] = 0;
  }
  bundles.clear();
}


// -- end -----------------------------------------------------------
