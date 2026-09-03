% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Int data sort.

#include bool.spec
#include pos1.spec
#include nat1.spec
#supertypeof Pos
#supertypeof Nat

sort Int <"int_">;

cons @cInt <"cint"> : Nat <"arg"> -> Int                 internal defined_by_rewrite_rules;
     @cNeg <"cneg"> : Pos <"arg"> -> Int                 internal defined_by_rewrite_rules;

map Nat2Int <"nat2int"> : Nat <"arg"> -> Int             external defined_by_rewrite_rules;
    Int2Nat <"int2nat"> : Int <"arg"> -> Nat             external defined_by_rewrite_rules;
    Pos2Int <"pos2int"> : Pos <"arg"> -> Int             external defined_by_rewrite_rules;
    Int2Pos <"int2pos"> : Int <"arg"> -> Pos             external defined_by_rewrite_rules;
    max <"maximum">:Pos <"left"> #Int <"right">->Pos     external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Pos <"right">->Pos     external defined_by_rewrite_rules;
    max <"maximum">:Nat <"left"> #Int <"right">->Nat     external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Nat <"right">->Nat     external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Int <"right">->Int     external defined_by_rewrite_rules;
    min <"minimum">:Int <"left"> #Int <"right">->Int     external defined_by_rewrite_rules;
    abs <"abs">:Int <"arg">->Nat                         external defined_by_rewrite_rules;
    - <"negate">:Pos <"arg">->Int                        external defined_by_rewrite_rules;
    - <"negate">:Nat <"arg">->Int                        external defined_by_rewrite_rules;
    - <"negate">:Int <"arg">->Int                        external defined_by_rewrite_rules;
    succ <"succ">:Int <"arg">->Int                       external defined_by_rewrite_rules;
    pred <"pred">:Nat <"arg">->Int                       external defined_by_rewrite_rules;
    pred <"pred">:Int <"arg">->Int                       external defined_by_rewrite_rules;
    + <"plus">:Int <"left"> #Int <"right">->Int          external defined_by_rewrite_rules;
    - <"minus">:Pos <"left"> # Pos <"right">->Int        external defined_by_rewrite_rules;
    - <"minus">:Nat <"left"> # Nat <"right">->Int        external defined_by_rewrite_rules;
    - <"minus">:Int <"left"> # Int <"right">->Int        external defined_by_rewrite_rules;
    * <"times">:Int <"left"> # Int <"right">->Int        external defined_by_rewrite_rules;
    div <"div">: Int <"left"> # Pos <"right"> -> Int     external defined_by_rewrite_rules;
    mod <"mod">:Int <"left"> # Pos <"right"> -> Nat      external defined_by_rewrite_rules;
    exp <"exp">:Int <"left"> # Nat <"right"> -> Int      external defined_by_rewrite_rules;

var b:Bool;
    n:Nat;
    m:Nat;
    p:Pos;
    q:Pos;
    x:Int;
    y:Int;

eqn ==(@cInt(m),@cInt(n)) = ==(m,n);
    ==(@cInt(n), @cNeg(p)) = false;
    ==(@cNeg(p), @cInt(n)) = false;
    ==(@cNeg(p),@cNeg(q)) = ==(p,q);
    <(@cInt(m),@cInt(n)) = <(m,n);
    <(@cInt(n),@cNeg(p)) = false;
    <(@cNeg(p),@cInt(n)) = true;
    <(@cNeg(p),@cNeg(q)) = <(q,p);
    <=(@cInt(m),@cInt(n)) = <=(m,n);
    <=(@cInt(n),@cNeg(p)) = false;
    <=(@cNeg(p),@cInt(n)) = true;
    <=(@cNeg(p),@cNeg(q)) = <=(q,p);
    less_total(x,y) = <(x,y);
    Nat2Int(n) = @cInt(n);
    Int2Nat(@cInt(n)) = n;
    Pos2Int(p) = @cInt(@cNat(p));
    Int2Pos(@cInt(n)) = Nat2Pos(n);
    max(p,@cInt(n)) = max(p,n);
    max(p,@cNeg(q)) = p;
    max(@cInt(n),p) = max(n,p);
    max(@cNeg(q),p) = p;
    max(m,@cInt(n)) = if(<=(m,n),n,m);
    max(n,@cNeg(p)) = n;
    max(@cInt(m),n) = if(<=(m,n),n,m);
    max(@cNeg(p),n) = n;
    max(x,y) = if(<=(x,y),y,x);
    min(x,y) = if(<=(x,y),x,y);
    abs(@cInt(n)) = n;
    abs(@cNeg(p)) = @cNat(p);
    -(p) = @cNeg(p);
    -(@c0) = @cInt(@c0);
    -(@cNat(p)) = @cNeg(p);
    -(@cInt(n)) = -(n);
    -(@cNeg(p)) = @cInt(@cNat(p));
    succ(@cInt(n)) = @cInt(@cNat(succ(n)));
    succ(@cNeg(p)) = -(pred(p));
    pred(@c0) = @cNeg(@c1);
    pred(@cNat(p)) = @cInt(pred(p));
    pred(@cInt(n)) = pred(n);
    pred(@cNeg(p)) = @cNeg(succ(p));
    +(@cInt(m),@cInt(n)) = @cInt(+(m,n));
    +(@cInt(n),@cNeg(p)) = -(n,@cNat(p));
    +(@cNeg(p),@cInt(n)) = -(n,@cNat(p));
    +(@cNeg(p),@cNeg(q)) = @cNeg(@addc(false,p,q));
    <=(q,p) -> -(p,q) = @cInt(@gtesubtb(false,p,q));
    <(p,q) -> -(p,q) = -(@gtesubtb(false,q,p));
    <=(n,m) -> -(m,n) = @cInt(@monus(m,n));
    <(m,n) -> -(m,n) = -(@monus(n,m));
    -(x,y) = +(x, -(y));
    *(@cInt(m),@cInt(n)) = @cInt(*(m,n));
    *(@cInt(n),@cNeg(p)) = -(*(@cNat(p),n));
    *(@cNeg(p),@cInt(n)) = -(*(@cNat(p),n));
    *(@cNeg(p),@cNeg(q)) = @cInt(@cNat(*(p,q)));
    div(@cInt(n),p) = @cInt(div(n,p));
    div(@cNeg(p),q) = @cNeg(succ(div(pred(p),q)));
    mod(@cInt(n),p) = mod(n,p);
    mod(@cNeg(p),q) = Int2Nat(-(q,succ(mod(pred(p),q))));
    exp(@cInt(m),n) = @cInt(exp(m,n));
    @even(n) -> exp(@cNeg(p),n) = @cInt(@cNat(exp(p,n)));
    !(@even(n)) -> exp(@cNeg(p),n) = @cNeg(exp(p,n));
