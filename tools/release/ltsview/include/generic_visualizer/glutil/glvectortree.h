#ifndef MCRL2_OPENGL_GLVECTORTREE_H
#define MCRL2_OPENGL_GLVECTORTREE_H

#include <functional>

namespace GlUtil
{
namespace VectorTree
{

template <typename T> struct VectorTree
{
  T data;
  VectorTree<T>* parent;
  std::vector<VectorTree<T>*> children;

  static std::function<void(VectorTree<T>*, VectorTree<T>*)> setParent;

  static std::function<void(VectorTree<T>*, VectorTree<T>*)> addChild;

  typedef typename std::vector<VectorTree<T>*>::iterator childIterator;

  static std::function<childIterator(VectorTree<T>*)> getChildBegin;
  static std::function<childIterator(VectorTree<T>*)> getChildEnd;
};

template <typename T>
std::function<void(VectorTree<T>*, VectorTree<T>*)> VectorTree<T>::setParent =
    [](VectorTree<T>* node, VectorTree<T>* parent) { node->parent = parent; };

template <typename T>
std::function<void(VectorTree<T>*, VectorTree<T>*)> VectorTree<T>::addChild =
    [](VectorTree<T>* node, VectorTree<T>* child)
{ node->children.emplace_back(child); };

template <typename T>
std::function<typename VectorTree<T>::childIterator(VectorTree<T>*)>
    VectorTree<T>::getChildBegin = [](VectorTree<T>* node)
{ return node->children.begin(); };

template <typename T>
std::function<typename VectorTree<T>::childIterator(VectorTree<T>*)>
    VectorTree<T>::getChildEnd = [](VectorTree<T>* node)
{ return node->children.end(); };
} // namespace VectorTree
} // namespace GlUtil
#endif