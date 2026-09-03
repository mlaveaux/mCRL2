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
#include pos64.spec
#include nat64.spec
#supertypeof Pos
#supertypeof Nat

sort Int <"int_">;

cons @cInt <"cint"> : Nat <"arg"> -> Int                      internal defined_by_rewrite_rules;
     @cNeg <"cneg"> : Pos <"arg"> -> Int                      internal defined_by_rewrite_rules;

map Nat2Int <"nat2int"> : Nat <"arg"> -> Int                  external defined_by_rewrite_rules;
    Int2Nat <"int2nat"> : Int <"arg"> -> Nat                  external defined_by_rewrite_rules;
    Pos2Int <"pos2int"> : Pos <"arg"> -> Int                  external defined_by_rewrite_rules;
    Int2Pos <"int2pos"> : Int <"arg"> -> Pos                  external defined_by_rewrite_rules;
    max <"maximum">:Pos <"left"> #Int <"right">->Pos          external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Pos <"right">->Pos          external defined_by_rewrite_rules;
    max <"maximum">:Nat <"left"> #Int <"right">->Nat          external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Nat <"right">->Nat          external defined_by_rewrite_rules;
    max <"maximum">:Int <"left"> #Int <"right">->Int          external defined_by_rewrite_rules;
    min <"minimum">:Int <"left"> #Int <"right">->Int          external defined_by_rewrite_rules;
    abs <"abs">:Int <"arg">->Nat                              external defined_by_rewrite_rules;
    - <"negate">:Pos <"arg">->Int                             external defined_by_rewrite_rules;
    - <"negate">:Nat <"arg">->Int                             external defined_by_rewrite_rules;
    - <"negate">:Int <"arg">->Int                             external defined_by_rewrite_rules;
    succ <"succ">:Int <"arg">->Int                            external defined_by_rewrite_rules;
    pred <"pred">:Nat <"arg">->Int                            external defined_by_rewrite_rules;
    pred <"pred">:Int <"arg">->Int                            external defined_by_rewrite_rules;
    + <"plus">:Int <"left"> #Int <"right">->Int               external defined_by_rewrite_rules;
    - <"minus">:Pos <"left"> # Pos <"right">->Int             external defined_by_rewrite_rules;
    - <"minus">:Nat <"left"> # Nat <"right">->Int             external defined_by_rewrite_rules;
    - <"minus">:Int <"left"> # Int <"right">->Int             external defined_by_rewrite_rules;
    * <"times">:Int <"left"> # Int <"right">->Int             external defined_by_rewrite_rules;
    div <"div">: Int <"left"> # Pos <"right"> -> Int          external defined_by_rewrite_rules;
    mod <"mod">:Int <"left"> # Pos <"right"> -> Nat           external defined_by_rewrite_rules;
    exp <"exp">:Int <"left"> # Nat <"right"> -> Int           external defined_by_rewrite_rules;

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
    less_total(x,y) = <(x,y);
    <=(@cInt(m),@cInt(n)) = <=(m,n);
    <=(@cInt(n),@cNeg(p)) = false;
    <=(@cNeg(p),@cInt(n)) = true;
    <=(@cNeg(p),@cNeg(q)) = <=(q,p);
    Nat2Int(n) = @cInt(n);
    Int2Nat(@cInt(n)) = n;
    Pos2Int(p) = @cInt(Pos2Nat(p));
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
    abs(@cNeg(p)) = Pos2Nat(p);
    -(p) = @cNeg(p);
    -(@c0) = @cInt(@c0);
%    -(n) = if(@equals_zero(n),@cInt(@most_significant_digitNat(@zero_word)),@cNeg(Nat2Pos(n)));
    @equals_zero(n) -> -(n) = @cInt(@most_significant_digitNat(@zero_word));
    @not_equals_zero(n) -> -(n) = @cNeg(Nat2Pos(n));
    -(@cInt(n)) = -(n);
    -(@cNeg(p)) = @cInt(Pos2Nat(p));
    succ(@cInt(n)) = @cInt(@succ_nat(n));  
    succ(@cNeg(p)) = -(pred(p));
    pred(n) = if(@equals_zero(n), @cNeg(@most_significant_digit(@one_word)), @cInt(@natpred(n)));  
    pred(@cInt(n)) = pred(n);
    pred(@cNeg(p)) = @cNeg(succ(p));
    +(@cInt(m),@cInt(n)) = @cInt(+(m,n));
    +(@cInt(n),@cNeg(p)) = -(n,Pos2Nat(p));
    +(@cNeg(p),@cInt(n)) = -(n,Pos2Nat(p));
    +(@cNeg(p),@cNeg(q)) = @cNeg(+(p,q));
    <=(q,p) -> -(p,q) = @cInt(@monus(Pos2Nat(p),Pos2Nat(q)));
    <(p,q) -> -(p,q) = -(@monus(Pos2Nat(q),Pos2Nat(p)));
    <=(n,m) -> -(m,n) = @cInt(@monus(m,n));
    <(m,n) -> -(m,n) = -(@monus(n,m));
    -(x,y) = +(x, -(y));
    *(@cInt(m),@cInt(n)) = @cInt(*(m,n));
    *(@cInt(n),@cNeg(p)) = -(*(Pos2Nat(p),n));
    *(@cNeg(p),@cInt(n)) = -(*(Pos2Nat(p),n));
    *(@cNeg(p),@cNeg(q)) = @cInt(Pos2Nat(*(p,q)));
    div(@cInt(n),p) = @cInt(div(n,p));
    div(@cNeg(p),q) = @cNeg(succ(div(pred(p),q)));
    mod(@cInt(n),p) = mod(n,p);
    mod(@cNeg(p),q) = Int2Nat(-(q,succ(mod(pred(p),q))));
    exp(@cInt(m),n) = @cInt(exp(m,n));
    exp(@cNeg(p),n) = if(@is_odd(n), @cNeg(exp(p,n)), @cInt(Pos2Nat(exp(p,n))));

