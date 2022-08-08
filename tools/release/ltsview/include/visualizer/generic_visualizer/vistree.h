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

#include "glprimitives.h"

#ifndef MCRL2_OPENGL_VISTREE_H
#define MCRL2_OPENGL_VISTREE_H



namespace VisTree{
    struct VisTreeNode{
        Primitives::Shape shape;
        QMatrix4x4 matrix;
        VisTreeNode* parent;
        std::vector<VisTreeNode*> children;
    };

    template < typename OutputTreeType, typename InputTreeType, typename Functor, typename ChildIterator >
    OutputTreeType* fold_tree(
        OutputTreeType* parent, 
        InputTreeType* current,
        std::function<void(OutputTreeType *, OutputTreeType *)> setParent, 
        std::function<void(OutputTreeType *, OutputTreeType *)>  addChild, 
        Functor& f, 
        std::function<ChildIterator(InputTreeType *)>& getChildIterator
    ){
        OutputTreeType* node = f(parent, current);
        setParent(node, parent);
        for(auto child_it = getChildIterator(current); child_it.it != child_it.it_end; ++child_it) 
            addChild(node, fold_tree(node, *child_it, setParent, addChild, f, getChildIterator));
        return node;
    }

    template < typename InputTreeType, typename Functor, typename ChildIterator >
    VisTree::VisTreeNode* generateVisTree(InputTreeType* root, Functor& f, std::function<ChildIterator(InputTreeType *)>& getChildIterator) {
        return fold_tree<VisTree::VisTreeNode, InputTreeType, Functor, ChildIterator>(nullptr, 
                                                                                      root,
                                                                                      [](VisTree::VisTreeNode* node, VisTree::VisTreeNode* parent){ node->parent = parent; },
                                                                                      [](VisTree::VisTreeNode* node, VisTree::VisTreeNode* child){ node->children.emplace_back(child); },
                                                                                      f,
                                                                                      getChildIterator);
    }
}




#endif
