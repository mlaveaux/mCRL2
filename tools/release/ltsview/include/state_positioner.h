// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATE_POSITIONER_H
#define STATE_POSITIONER_H

class LTS;
class State;

class StatePositioner
{
  public:
    StatePositioner(LTS* l): lts(l) {}
    virtual ~StatePositioner() {}
    virtual void positionStates() = 0;

  protected:
    LTS* lts;
};

#endif
