% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Pos data sort.

#include bool.spec

sort Pos <"pos">;

cons @c1 <"c1">:Pos                                                               internal defined_by_rewrite_rules;
     @cDub <"cdub">:Bool <"left"> #Pos <"right"> ->Pos                            internal defined_by_rewrite_rules;

map max <"maximum">:Pos <"left"> #Pos <"right">->Pos                              external defined_by_rewrite_rules;
    min <"minimum">:Pos <"left"> #Pos <"right">->Pos                              external defined_by_rewrite_rules;
    succ <"succ">:Pos <"arg">->Pos                                                external defined_by_rewrite_rules;
    @pospred <"pos_predecessor">:Pos <"arg">->Pos                                 internal defined_by_rewrite_rules;
    + <"plus">:Pos <"left"> #Pos <"right">->Pos                                   external defined_by_rewrite_rules;
    @addc <"add_with_carry">:Bool <"arg1"> #Pos <"arg2"> #Pos <"arg3">->Pos       internal defined_by_rewrite_rules;
    * <"times">:Pos <"left"> #Pos <"right">->Pos                                  external defined_by_rewrite_rules;
    @powerlog2 <"powerlog2_pos">:Pos <"arg"> -> Pos                               internal defined_by_rewrite_rules;

var b:Bool;
    c:Bool;
    p:Pos;
    q:Pos;
    p1:Pos;
    q1:Pos;
eqn ==(@c1, @cDub(b,p)) = false;
    ==(@cDub(b,p), @c1) = false;
% The tools run 10% faster with the rules below, compared to the rules in the book of 2014. 
    ==(@cDub(b,p), @cDub(b, q)) =  ==(p,q);
    ==(@cDub(false,p), @cDub(true, q)) = false;
    ==(@cDub(true,p), @cDub(false, q)) = false;
    ==(succ(p),@c1) = false;
    ==(@c1,succ(q)) = false;
    ==(succ(p),@cDub(c,q)) = ==(p,@pospred(@cDub(c,q)));
    ==(@cDub(b,p),succ(q)) = ==(@pospred(@cDub(b,p)),q);
    <(p, @c1) = false;
    <(@c1, @cDub(b,p)) = true;
% The equation below is required for the enumeration of lists
    <(@cDub(b,p), @cDub(c,q)) = if(=>(c,b), <(p,q), <=(p,q));
    <(succ(p),@cDub(c,q)) = <(p,@pospred(@cDub(c,q)));
    <(@cDub(b,p),succ(q)) = <=(@cDub(b,p),q);
    <(@c1,succ(q)) = true;
    <=(@c1, p) = true;
    <=(@cDub(b,p), @c1) = false;
% The equation below is required for the enumeration of lists
    <=(@cDub(b,p), @cDub(c,q)) = if(=>(b,c), <=(p,q), <(p,q));
    <=(succ(p),@cDub(c,q)) = <(p,@cDub(c,q));
    <=(@cDub(b,p),succ(q)) = <=(@pospred(@cDub(b,p)),q);
    <=(succ(p),@c1) = false;
    less_total(p,q) = <(p,q);
    max(p,q) = if(<=(p,q),q,p);
    min(p,q) = if(<=(p,q),p,q);
    succ(@c1) = @cDub(false, @c1);
    succ(@cDub(false,p)) = @cDub(true,p);
    succ(@cDub(true,p)) = @cDub(false,succ(p));
    @pospred(@c1) = @c1;
    @pospred(@cDub(false,@c1)) = @c1;
    @pospred(@cDub(false,@cDub(b,p))) = @cDub(true,@pospred(@cDub(b,p)));
    @pospred(@cDub(true,p)) = @cDub(false,p);
    +(p,q) = @addc(false,p,q);
    @addc(false,@c1,p) = succ(p);
    @addc(true,@c1,p) = succ(succ(p));
    @addc(false,p,@c1) = succ(p);
    @addc(true,p,@c1) = succ(succ(p));
    @addc(b,@cDub(c,p),@cDub(c,q)) = @cDub(b,@addc(c,p,q));
    @addc(b,@cDub(false,p),@cDub(true,q)) = @cDub(!(b), @addc(b,p,q));
    @addc(b,@cDub(true,p),@cDub(false,q)) = @cDub(!(b), @addc(b,p,q));
    *(@c1,p) = p;
    *(p,@c1) = p;
    *(@cDub(false,p),q) = @cDub(false,*(p,q));
    *(p,@cDub(false,q)) = @cDub(false,*(p,q));
    *(@cDub(true,p),@cDub(true,q)) = @cDub(true,@addc(false,p,@addc(false,q,@cDub(false,*(p,q)))));
    @powerlog2(@c1) = @c1;
    @powerlog2(@cDub(b,@c1)) = @c1;
    @powerlog2(@cDub(b,@cDub(c,p))) = @cDub(false,@powerlog2(p));
