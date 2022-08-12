// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <vector>
#include <QMatrix4x4>

#include "glvectortree.h"
#include "glshapes.h"
#include "gltree.h"

#ifndef MCRL2_OPENGL_VISTREE_H
#define MCRL2_OPENGL_VISTREE_H

namespace VisTree
{
/// TODO: Smart pointer
struct VisTreeData{
  GlUtil::Shape* shape;
  QMatrix4x4 matrix;
};

using VisTreeNode = GlUtil::PointerTree::PointerTree<VisTreeData>;

template <typename InputTreeType, typename Functor, typename ChildIterator>
VisTree::VisTreeNode*
generateVisTree(InputTreeType* root, Functor& f,
                std::function<ChildIterator(InputTreeType*)>& getChildBegin,
                std::function<ChildIterator(InputTreeType*)>& getChildEnd)
{
  return GlUtil::fold_tree<VisTreeNode, InputTreeType, Functor,
                           ChildIterator>(nullptr, root, VisTreeNode::setParent, VisTreeNode::addChild,
                           f, getChildBegin, getChildEnd);
}
} // namespace VisTree

#endif
