// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file conedb.cpp
/// \brief Source file for the ConeDB class.

#include "conedb.h"
#include "mathutils.h"

#define ohash(k1,k2,b,m) ((k1*11408669 + k2*97416181 + b*71053447) & m)
#define thash(k,tb,m) ((k*11408669 + tb*97416181) & m)
#define BOT_BIT 1
#define TOP_BIT 2

ConeDB::ConeDB()
{
  int hashclass = 6;
  ohashtable.assign((1 << hashclass),-1);
  thashtable.assign((1 << hashclass),-1);
}

ConeDB::~ConeDB()
{
}

void ConeDB::addTruncatedCone(float r,bool t,bool b,int c)
{
  int k = compute_key(r);
  unsigned char tb = combine_top_bot(t,b);
  int i = find_tbucket(k,tb);
  if (i != -1)
  {
    // c will be the new cone belonging to k and tb
    tbuckets[i].cone = c;
  }
  else
  {
    check_thashtable();
    i = thash(k,tb,(thashtable.size()-1));
    tcone_bucket cb = { k, c, thashtable[i], tb };
    thashtable[i] = static_cast<int>(tbuckets.size());
    tbuckets.push_back(cb);
  }
}

int ConeDB::findTruncatedCone(float r,bool t,bool b)
{
  int k = compute_key(r);
  unsigned char tb = combine_top_bot(t,b);
  int i = find_tbucket(k,tb);
  if (i == -1)
  {
    return -1;
  }
  else
  {
    return tbuckets[i].cone;
  }
}

void ConeDB::check_thashtable()
{
  if (4*tbuckets.size() >= 3*thashtable.size())
  {
    // hash table has become too full; double its capacity and rehash
    thashtable.assign(2*thashtable.size(),-1);
    unsigned int i,h;
    for (i = 0; i < tbuckets.size(); ++i)
    {
      h = thash(tbuckets[i].key,tbuckets[i].top_bot,(thashtable.size()-1));
      tbuckets[i].next = thashtable[h];
      thashtable[h] = i;
    }
  }
}

int ConeDB::find_tbucket(int k,unsigned char tb)
{
  int h = thash(k,tb,(thashtable.size()-1));
  for (h = thashtable[h]; h != -1; h = tbuckets[h].next)
  {
    if (tbuckets[h].key == k  &&  tbuckets[h].top_bot == tb)
    {
      return h;
    }
  }
  return -1;
}

void ConeDB::addObliqueCone(float a,float r,float s,int c)
{
  int k1 = compute_key(a);
  int k2 = compute_key(r);
  bool b = s > 0.0f;
  int i = find_obucket(k1,k2,b);
  if (i != -1)
  {
    // c will be the new cone belonging to k and tb
    obuckets[i].cone = c;
  }
  else
  {
    check_ohashtable();
    i = ohash(k1,k2,b,(ohashtable.size()-1));
    ocone_bucket cb = { k1, k2, b, c, ohashtable[i] };
    ohashtable[i] = static_cast<int>(obuckets.size());
    obuckets.push_back(cb);
  }
}

int ConeDB::findObliqueCone(float a,float r,float s)
{
  int i = find_obucket(compute_key(a),compute_key(r),s > 0.0f);
  if (i == -1)
  {
    return -1;
  }
  else
  {
    return obuckets[i].cone;
  }
}

void ConeDB::check_ohashtable()
{
  if (4*obuckets.size() >= 3*ohashtable.size())
  {
    // hash table has become too full; double its capacity and rehash
    ohashtable.assign(2*ohashtable.size(),-1);
    unsigned int i,h;
    for (i = 0; i < obuckets.size(); ++i)
    {
      h = ohash(obuckets[i].alpha,obuckets[i].radius,obuckets[i].sign,
                (ohashtable.size()-1));
      obuckets[i].next = ohashtable[h];
      ohashtable[h] = i;
    }
  }
}

int ConeDB::find_obucket(int k1,int k2,bool b)
{
  int h = ohash(k1,k2,b,(ohashtable.size()-1));
  for (h = ohashtable[h]; h != -1; h = obuckets[h].next)
  {
    if (obuckets[h].alpha == k1 && obuckets[h].radius == k2 && obuckets[h].sign == b)
    {
      return h;
    }
  }
  return -1;
}

int ConeDB::compute_key(float r)
{
  return MathUtils::round_to_int(static_cast<float>(r * 100.0));
}

unsigned char ConeDB::combine_top_bot(bool t,bool b)
{
  unsigned char result = 0;
  if (t)
  {
    result |= TOP_BIT;
  }
  if (b)
  {
    result |= BOT_BIT;
  }
  return result;
}
