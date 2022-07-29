// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SP_STATE_POSITIONER_H
#define SP_STATE_POSITIONER_H

#include "state_positioner.h"

class LTS;

class SinglePassStatePositioner: public StatePositioner
{
  public:
    SinglePassStatePositioner(LTS* l);
    ~SinglePassStatePositioner();
    void positionStates();
};

#endif
