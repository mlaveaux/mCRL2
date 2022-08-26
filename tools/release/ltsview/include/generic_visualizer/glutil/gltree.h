#ifndef MCRL2_OPENGL_TREE_H
#define MCRL2_OPENGL_TREE_H

#include <functional>

namespace GlUtil{
template <typename OutputTreeType, typename InputTreeType, typename Functor,
          typename ChildIterator>
OutputTreeType*
fold_tree(OutputTreeType* parent, InputTreeType* current,
          std::function<void(OutputTreeType*, OutputTreeType*)>& setParent,
          std::function<void(OutputTreeType*, OutputTreeType*)>& addChild,
          Functor& f,
          std::function<ChildIterator(InputTreeType*)>& getChildBegin,
          std::function<ChildIterator(InputTreeType*)>& getChildEnd)
{
    /// TODO: Assert
  if (!current)
    return nullptr;
  OutputTreeType* node = f(parent, current);
  setParent(node, parent);
  auto child_it = getChildBegin(current);
  auto child_it_end = getChildEnd(current);
  for (; child_it != child_it_end; ++child_it)
    addChild(node, fold_tree(node, *child_it, setParent, addChild, f,
                             getChildBegin, getChildEnd));
  return node;
}
}
#endif