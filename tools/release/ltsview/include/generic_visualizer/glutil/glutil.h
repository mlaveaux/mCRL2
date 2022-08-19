// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef MCRL2_OPENGL_GLUTIL_H
#define MCRL2_OPENGL_GLUTIL_H

#include "glcalcarc.h"
#include "glluts.h"
#include "glmeshes.h"
#include "glprimitivefactories.h"
#include "glshapes.h"
#include "gltree.h"
#include "glvectortree.h"

#ifndef NDEBUG
#define msgAssert(expr, msg) if (!expr) { std::cerr << "Assert failed: " << msg << std::endl; assert(expr); }
#else
#define msgAssert(expr, msg) do {} while(0); 
#endif

#endif