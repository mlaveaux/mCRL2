// Copyright (c) 2009-2013 University of Twente
// Copyright (c) 2009-2013 Michael Weber <michaelw@cs.utwente.nl>
// Copyright (c) 2009-2013 Maks Verver <maksverver@geocities.com>
// Copyright (c) 2009-2013 Eindhoven University of Technology
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_PG_LIFTING_STRATEGY_H
#define MCRL2_PG_LIFTING_STRATEGY_H

#include "mcrl2/pg/ParityGame.h"
#include "mcrl2/pg/RefCounted.h"

#include <string>

class SmallProgressMeasures;

/*! \defgroup LiftingStrategies Lifting strategies

    Lifting strategies for the Small Progress Measures solving algorithm.
*/


/*! \ingroup LiftingStrategies

    Version 1 API for lifting strategies.

    Instances of this class encapsulate vertex lifting strategies to be used
    with the small progress measures parity game solver.

    These are expected to initialize and maintain the set of dirty vertices
    themselves, and inform the SPM implementation when the set becomes empty by
    returning NO_VERTEX from next().
*/
class LiftingStrategy
{
public:
  virtual ~LiftingStrategy() = default;

  /*! Record that the given vertex was lifted: */
  virtual void lifted(verti vertex) = 0;

  /*! Select the next vertex to lift. This method is called repeatedly by the
      SPM solver until it returns NO_VERTEX to indicate the solution is
      complete.

      \see lifted(verti vertex)
  */
  virtual verti next() = 0;
};


/*! \ingroup LiftingStrategies

    Version 2 API for lifting strategies.

    Instances of this class encapsulate vertex lifting strategies to be used
    with the small progress measures parity game solver.

    The SPM implementation manages the set of dirty vertices; it will call
    push() to add a vertex to the set, pop() to remove one (which one is at the
    discretion of the lifting strategy implementation) and bump() when an
    already-queued vertex's progress measure has been increased.
*/
class LiftingStrategy2
{
public:
  virtual ~LiftingStrategy2() = default;

  // TODO: document each separately like LiftingStrategy above
  virtual void push(verti vertex) = 0;
  virtual void bump(verti vertex) = 0;
  virtual verti pop() = 0;
};


/*! \ingroup LiftingStrategies
    Abstract base class for lifting strategy factories. */
class LiftingStrategyFactory : public RefCounted
{
public:
  ~LiftingStrategyFactory() override;

  /*! Returns pre-formatted plain-text documentation of the description
      strings accepted by create(), intended to be shown to the user.
      \see create(const std::string &description) */
  static const char* usage();

  /*! Creates a lifting strategy factory from a string description.
  \returns A factory object or NULL if the description could not be parsed.
  \see usage() for a description of available format strings. */
  static LiftingStrategyFactory* create(const std::string& description);

  virtual bool supports_version(int version) { return version == 1; }

  /*! Create a lifting strategy for the given game, to be used by the given
      Small Progress Measures solver. */
  virtual LiftingStrategy* create(const ParityGame& game, const SmallProgressMeasures& spm) = 0;

  virtual LiftingStrategy2* create2(const ParityGame& /*game*/, const SmallProgressMeasures& /*spm*/)
  {
    return nullptr;
    }
};


#endif /* ndef MCRL2_PG_LIFTING_STRATEGY_H */
