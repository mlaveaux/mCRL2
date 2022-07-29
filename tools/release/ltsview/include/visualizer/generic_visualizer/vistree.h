// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
#pragma once
#include <vector>
#include <QMatrix4x4>

#ifndef MCRL2_OPENGL_VISTREE_H
#define MCRL2_OPENGL_VISTREE_H

struct Shape;


namespace VisTree{
    
    enum ShapeType{
        SPHERE,
        HEMISPHERE,
        TUBE,
        CONE,
        TRUNCATED_CONE,
        OBLIQUE_CONE,
        DISC,
        RING,
    };

    struct Shape{
        ShapeType shape;
        QMatrix4x4 matrix;
    };

    struct VisTreeNode{
        Shape shape;
        VisTreeNode* parent;
        std::vector<VisTreeNode*> children;
    };


    namespace Shapes{
        struct Sphere : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::SPHERE;
            float radius;
        };
        struct HemiSphere : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::HEMISPHERE;
            float radius;
        };

        // without transform a tube is formed by circle in the xy0 and xy1 planes 
        struct Tube : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::TUBE;
            float radius;
            float height;
        };

        // without transform a cone is formed by a cirle in the xy0 plane and a point at 001 
        struct Cone : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::CONE;
            float radius;
            float height;
        };

        // without transform a disc is formed by a circle in the xy0 plane (filled)
        struct Disc : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::DISC;
            float radius;
        };

        // without transform a disc is formed by a circle in the xy0 plane (not filled)
        struct Ring : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::RING;
            float radius;
        };

        // Suppose we have an untransformed cone with a certain radius and height
        // We can truncate the cone by 'cutting' it open using two planes at z=a and z=b
        // This results in 4 options: 
        //  - a and b both miss -> regular cone
        //  - bottom is sliced off -> truncate
        //  - top is sliced off -> truncate
        //  - both are sliced off -> truncate
        // For every case where we truncate we can either fill or not fill the resulting cut
        // i.e. ring/disc
        struct TruncatedCone : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::TRUNCATED_CONE;
            float radius_bot, radius_top;
            float height;
            float bot, top; // if bot \in (0, height) a cut occurs
            bool fill_bot, fill_top;
        };


        // Think of a regular cone, but move the tip off-axis
        // In this case, untransformed:
        //   - tip = height*(cos(alpha), sin(alpha), 0) with alpha in [0, 90] degrees
        struct ObliqueCone : Shape{
            VisTree::ShapeType shape = VisTree::ShapeType::OBLIQUE_CONE;
            float radius;
            float height;
            float alpha; // in degrees
        };
    }


    template < typename InputTreeType, typename Functor, typename Iterator >
    void recurse(VisTreeNode* parent, InputTreeType* current, Functor& f, std::function<Iterator(InputTreeType *)>& getChildIterator){
        VisTreeNode* node = f(parent, current);
        node->parent = parent;
        parent->children.emplace_back(node);
        for(auto child_it = getChildIterator(current); child_it.it != child_it.it_end; ++child_it) 
            recurse(node, *child_it, f, getChildIterator);
    }

    template < typename InputTreeType, typename Functor, typename Iterator >
    VisTree::VisTreeNode* generateVisTree(InputTreeType* root, Functor& f, std::function<Iterator(InputTreeType *)>& getChildIterator){
        VisTreeNode* node = f(nullptr, root);
        node->parent = nullptr; // this is the root of the vistree
        for(auto child_it = getChildIterator(root); child_it.it != child_it.it_end; ++child_it) 
            VisTree::recurse(node, *child_it, f, getChildIterator);
        return node;
    }
}




#endif
