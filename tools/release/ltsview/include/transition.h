// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file transition.h
/// \brief Header file for Transition class

#ifndef TRANSITION_H
#define TRANSITION_H

class State;

class Transition
{
  public:
    Transition(State* bs,State* es,int lbl);
    ~Transition();
    State*  getBeginState() const;
    State*  getEndState() const;
    int     getLabel() const;
    bool    isBackpointer() const;
    bool    isSelfLoop() const;

  private:
    State*  beginState;
    State*  endState;
    int     label;
};

#endif
