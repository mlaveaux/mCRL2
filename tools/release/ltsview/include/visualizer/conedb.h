// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file conedb.h
/// \brief Header file for the ConeDB class.

#ifndef CONEDB_H
#define CONEDB_H
#include <vector>

class ConeDB
{
  public:
    ConeDB();
    ~ConeDB();
    void addObliqueCone(float a,float r,float s,int c);
    int findObliqueCone(float a,float r,float s);
    void addTruncatedCone(float r,bool t,bool b,int c);
    int findTruncatedCone(float r,bool t,bool b);
  private:
    struct ocone_bucket
    {
      int alpha;
      int radius;
      bool sign;
      int cone;
      int next;
    };
    struct tcone_bucket
    {
      int key;
      int cone;
      int next;
      unsigned char top_bot;
    };
    std::vector<int> ohashtable;
    std::vector<int> thashtable;
    std::vector<ocone_bucket> obuckets;
    std::vector<tcone_bucket> tbuckets;

    void check_ohashtable();
    void check_thashtable();
    int find_obucket(int k1,int k2,bool b);
    int find_tbucket(int k,unsigned char tb);
    int compute_key(float r);
    unsigned char combine_top_bot(bool t,bool b);
};

#endif
