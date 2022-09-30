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
#include "mcrl2/utilities/logger.h"

#include <QOpenGLFunctions>


#ifndef NDEBUG
#include "windows.h"
#define msgAssert(expr, msg) if (!expr) { mCRL2log(mcrl2::log::error) << "Assert failed: " << msg << std::endl; assert(expr); }
inline GLenum glCheckError_(const char *file, int line)
{
    GLenum errorCode;
    while ((errorCode = QOpenGLContext::currentContext()->functions()->glGetError()) != GL_NO_ERROR)
    {
        std::string error;
        switch (errorCode)
        {
            case GL_INVALID_ENUM:                  error = "INVALID_ENUM"; break;
            case GL_INVALID_VALUE:                 error = "INVALID_VALUE"; break;
            case GL_INVALID_OPERATION:             error = "INVALID_OPERATION"; break;
            case GL_STACK_OVERFLOW:                error = "STACK_OVERFLOW"; break;
            case GL_STACK_UNDERFLOW:               error = "STACK_UNDERFLOW"; break;
            case GL_OUT_OF_MEMORY:                 error = "OUT_OF_MEMORY"; break;
            case GL_INVALID_FRAMEBUFFER_OPERATION: error = "INVALID_FRAMEBUFFER_OPERATION"; break;
        }
        mCRL2log(mcrl2::log::error) << error << " | " << file << " (" << line << ")" << std::endl;
    }
    return errorCode;
}
#define glCheckError() glCheckError_(__FILE__, __LINE__) 
#define execAssertGL(expr, msg) { bool execAssertGl_result = expr; msgAssert(execAssertGl_result, msg);}
#else
#define noop do {} while(0)
#define msgAssert(expr, msg) noop 
#define glCheckError() noop
#define execAssertGL(expr, msg) noop
#endif

#endif