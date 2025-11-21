% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Bool data sort.

sort Bool <"bool_">;

cons true <"true_"> : Bool                                         external defined_by_rewrite_rules;
     false <"false_"> : Bool                                       external defined_by_rewrite_rules;

map ! <"not_"> : Bool <"arg"> -> Bool                                   external defined_by_rewrite_rules;
    && <"and_"> : Bool <"left"> # Bool <"right"> -> Bool                external defined_by_rewrite_rules;
    || <"or_"> : Bool <"left"> # Bool <"right"> -> Bool                 external defined_by_rewrite_rules;
    => <"implies"> : Bool <"left"> # Bool <"right"> -> Bool             external defined_by_rewrite_rules;

var b:Bool;
    c:Bool;
eqn !(true) = false;
    !(false) = true;
    !(!(b)) = b;
    &&(b,true) = b;
    &&(b,false) = false;
    &&(true,b) = b;
    &&(false,b) = false;
    ||(b,true) = true;
    ||(b,false) = b;
    ||(true,b) = true;
    ||(false,b) = b;
    =>(b,true) = true;
    =>(b,false) = !(b);
    =>(true,b) = b;
    =>(false,b) = true;
    ==(true,b) = b;
    ==(false,b) = !(b);
    ==(b,true) = b;
    ==(b,false) = !(b);
    <(false,b) = b;
    <(true,b) = false;
    <(b,false) = false;
    <(b,true) = !(b);
    <=(false,b) = true;
    <=(true,b) = b;
    <=(b,false) = !(b);
    <=(b,true) = true;
    less_total(b,c) = <(b,c);


