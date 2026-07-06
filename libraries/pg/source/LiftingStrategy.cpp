// Copyright (c) 2009-2013 University of Twente
// Copyright (c) 2009-2013 Michael Weber <michaelw@cs.utwente.nl>
// Copyright (c) 2009-2013 Maks Verver <maksverver@geocities.com>
// Copyright (c) 2009-2013 Eindhoven University of Technology
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/pg/PredecessorLiftingStrategy.h"
#include "mcrl2/pg/FocusListLiftingStrategy.h"
#include "mcrl2/pg/MaxMeasureLiftingStrategy.h"
#include "mcrl2/pg/OldMaxMeasureLiftingStrategy.h"
#include "mcrl2/pg/LinPredLiftingStrategy.h"

#include <cstdlib>

/// \brief Convert string to lowercase.
/// \param s An input string
/// \return The text in \a in in lowercase
static inline
std::string tolower(std::string s)
{
  // dirty trick in the function to disambiguate between two tolower functions!
  std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
  return s;
}

LiftingStrategyFactory::~LiftingStrategyFactory() = default;

const char *LiftingStrategyFactory::usage()
{
    return
"linear:alternate\n"
"   Use a linear lifting strategy (aka swiping).\n"
"   - alternate: if 1, switches direction between forward and backward \n"
"     whenever the end of the list is reached (default: 0)\n"
"\n"
"predecessor:stack\n"
"   Use a predecessor lifting strategy (aka worklist).\n"
"   - stack: if 1, removes elements from the end of the queue instead \n"
"            of the beginning (default: 0)\n"
"\n"
"focuslist:alternate:max_size:lift_ratio\n"
"   Use swiping + focus list lifting strategy.\n"
"   - alternate: see 'linear' (default: 0)\n"
"   - max_size: the maximum size of the focus list, either as an absolute size\n"
"     greater than 1, or as a ratio (between 0 and 1) of the total number of\n"
"     vertices (default: 0.1)\n"
"   - lift_ratio: the maximum number of lifting attempts performed on the\n"
"     focus list before switching back to swiping, as a ratio of the maximum\n"
"     focus list size (default: 10.0)\n"
"\n"
"maxmeasure:order\n"
"   Maximum measure propagation; a variant of the predecessor lifting strategy\n"
"   that prefers to lift vertices with higher progress measures.\n"
"   - order: tie-breaking lifting order: 0 (queue-like), 1 (stack-like)\n"
"            or 2 (heap order) (default: 2)\n"
"\n"
"maxstep:order\n"
"   Maximum step variant of maximum measure propagation.\n"
"   - order: see 'maxmeasure'\n"
"\n"
"minmeasure:order\n"
"   Minimum measure propagation.\n"
"   - order: see 'maxmeasure'\n"
"\n"
"oldmaxmeasure\n"
"   Old implementation of max. measure lifting strategy.\n"
"   Included for regression testing purposes only.\n"
"\n"
"linpred\n"
"   Obsolete (roughly equivalent to predecessor:0:0).\n";
}

LiftingStrategyFactory *
    LiftingStrategyFactory::create(const std::string &description)
{
  if (description.empty())
  {
    return nullptr;
  }

    // Split into parts, separated by semicolon characters
    std::vector<std::string> parts;
    std::string::size_type i;
    std::string::size_type j;
    for (i = 0; (j = description.find(':', i)) != std::string::npos; i = j + 1)
    {
      parts.emplace_back(description, i, j - i);
    }
    parts.emplace_back(description, i, j);

    std::string case_ = tolower(parts[0]);
    if ( case_ == "linear" || case_ == "lin" )
    {
        bool alternate = parts.size() > 1 && 0 != atoi(parts[1].c_str());
        return new LinearLiftingStrategyFactory(alternate);
    }
    else if ( case_ == "predecessor" || case_ == "pred" )
    {
        bool stack    = parts.size() > 1 && 0 != atoi(parts[1].c_str());
        return new PredecessorLiftingStrategyFactory(stack);
    }
    else if ( case_ == "focuslist" || case_ == "focus" )
    {
        bool alternate    = parts.size() > 1 && 0 != atoi(parts[1].c_str());
        double max_size   = (parts.size() > 2 ? atof(parts[2].c_str()) : 0);
        double lift_ratio = (parts.size() > 3 ? atof(parts[3].c_str()) : 0);
        return new FocusListLiftingStrategyFactory(
            alternate, max_size, lift_ratio );
    }
    else if (case_ == "maxmeasure")
    {
        int order = (parts.size() > 1 ? atoi(parts[1].c_str()) : 2);
        return new MaxMeasureLiftingStrategyFactory(
            (MaxMeasureLiftingStrategy2::Order)order,
            MaxMeasureLiftingStrategy2::MAX_VALUE);
    }
    else if (case_ == "maxstep")
    {
        int order = (parts.size() > 1 ? atoi(parts[1].c_str()) : 2);
        return new MaxMeasureLiftingStrategyFactory(
            (MaxMeasureLiftingStrategy2::Order)order, 
            MaxMeasureLiftingStrategy2::MAX_STEP );
    }
    else if (case_ == "minmeasure")
    {
        int order = (parts.size() > 1 ? atoi(parts[1].c_str()) : 2);
        return new MaxMeasureLiftingStrategyFactory(
            (MaxMeasureLiftingStrategy2::Order)order, 
            MaxMeasureLiftingStrategy2::MIN_VALUE );
    }
    else if (case_ == "oldmaxmeasure")
    {
        return new OldMaxMeasureLiftingStrategyFactory();
    }
    else if (case_ == "linpred")
    {
        return new LinPredLiftingStrategyFactory();
    }
    else
    {
        // No suitable strategy found
        return nullptr;
    }
}
