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

namespace PointerTree
{
template <typename T> struct PointerTree; // forward declaration

template <typename T> struct PointerTreeForwardIterator
{
  PointerTree<T>* node;
  void operator++();
  PointerTree<T>* operator*();
};

template <typename T>
bool operator!=(const PointerTreeForwardIterator<T>& a, const PointerTreeForwardIterator<T>& b){
    return a.node != b.node;
}

template <typename T> struct PointerTree
{
  T data;
  PointerTree<T>* parent;
  PointerTree<T>* next;
  PointerTree<T>* firstChild;
  PointerTree<T>* lastChild;
  static std::function<void(PointerTree<T>*, PointerTree<T>*)> setParent;
  static std::function<void(PointerTree<T>*, PointerTree<T>*)> addChild;
  typedef typename PointerTreeForwardIterator<T> childIterator;
  static std::function<childIterator(PointerTree<T>*)> getChildBegin;
  static std::function<childIterator(PointerTree<T>*)> getChildEnd;
};

template <typename T> void PointerTreeForwardIterator<T>::operator++()
{
  assert(node != nullptr); // should never call ++ when node == nullptr
  node = node->next;
}

template <typename T> PointerTree<T>* PointerTreeForwardIterator<T>::operator*()
{
  assert(node != nullptr); // should never dereference when node == nullptr
  return node;
}


template <typename T>
std::function<void(PointerTree<T>*, PointerTree<T>*)>
    PointerTree<T>::setParent = [](PointerTree<T>* node, PointerTree<T>* parent)
{ node->parent = parent; };

template <typename T>
std::function<void(PointerTree<T>*, PointerTree<T>*)> PointerTree<T>::addChild =
    [](PointerTree<T>* node, PointerTree<T>* child)
{
  if (node->firstChild == nullptr)
  {
    node->firstChild = child;
    node->lastChild = child;
  }
  else
  {
    node->lastChild->next = child;
    node->lastChild = child;
  }
};

template <typename T>
std::function<typename PointerTree<T>::childIterator(PointerTree<T>*)>
    PointerTree<T>::getChildBegin = [](PointerTree<T>* node)
{ PointerTree<T>::childIterator it;
it.node = node->firstChild; return it; };

template <typename T>
std::function<typename PointerTree<T>::childIterator(PointerTree<T>*)>
    PointerTree<T>::getChildEnd = [](PointerTree<T>* node)
{ PointerTree<T>::childIterator it; 
it.node = nullptr; return it; };
}; // namespace PointerTree

} // namespace GlUtil
#endif