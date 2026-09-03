% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Real data sort.

#include bool.spec
#include pos1.spec
#include nat1.spec
#include int1.spec
#supertypeof Pos
#supertypeof Nat
#supertypeof Int

sort Real <"real_">;

map @cReal <"creal">: Int <"left"> # Pos <"right"> -> Real     internal defined_by_rewrite_rules;
    Pos2Real <"pos2real"> : Pos <"arg"> -> Real                external defined_by_rewrite_rules;
    Nat2Real <"nat2real"> : Nat <"arg"> -> Real                external defined_by_rewrite_rules;
    Int2Real <"int2real"> : Int <"arg"> -> Real                external defined_by_rewrite_rules;
    Real2Pos <"real2pos"> : Real <"arg"> -> Pos                external defined_by_rewrite_rules;
    Real2Nat <"real2nat"> : Real <"arg"> -> Nat                external defined_by_rewrite_rules;
    Real2Int <"real2int"> : Real <"arg"> -> Int                external defined_by_rewrite_rules;
    max <"maximum"> : Real <"left"> # Real <"right"> -> Real   external defined_by_rewrite_rules;
    min <"minimum"> : Real <"left"> # Real <"right"> -> Real   external defined_by_rewrite_rules;
    abs <"abs"> : Real <"arg"> -> Real                         external defined_by_rewrite_rules;
    - <"negate"> : Real <"arg"> -> Real                        external defined_by_rewrite_rules;
    succ <"succ"> : Real <"arg"> -> Real                       external defined_by_rewrite_rules;
    pred <"pred"> : Real <"arg"> -> Real                       external defined_by_rewrite_rules;
    + <"plus"> : Real <"left"> # Real <"right"> -> Real        external defined_by_rewrite_rules;
    - <"minus"> : Real <"left"> # Real <"right"> -> Real       external defined_by_rewrite_rules;
    * <"times"> : Real <"left"> # Real <"right"> -> Real       external defined_by_rewrite_rules;
    exp <"exp"> : Real <"left"> # Int <"right"> -> Real        external defined_by_rewrite_rules;
    / <"divides"> : Pos <"left"> # Pos <"right"> -> Real       external defined_by_rewrite_rules;
    / <"divides"> : Nat <"left"> # Nat <"right"> -> Real       external defined_by_rewrite_rules;
    / <"divides"> : Int <"left"> # Int <"right"> -> Real       external defined_by_rewrite_rules;
    / <"divides"> : Real <"left"> # Real <"right"> -> Real     external defined_by_rewrite_rules;
    floor <"floor"> : Real <"arg"> -> Int                      external defined_by_rewrite_rules;
    ceil <"ceil"> : Real <"arg"> -> Int                        external defined_by_rewrite_rules;
    round <"round"> : Real <"arg"> -> Int                      external defined_by_rewrite_rules;
    @redfrac <"reduce_fraction"> : Int <"left"> # Int <"right"> -> Real        
                                                               internal defined_by_rewrite_rules;
    @redfracwhr <"reduce_fraction_where"> : Pos <"arg1"> # Int <"arg2"> # Nat <"arg3"> -> Real
                                                               internal defined_by_rewrite_rules;
    @redfrachlp <"reduce_fraction_helper"> : Real <"left"> # Int <"right"> -> Real
                                                               internal defined_by_rewrite_rules;

var m:Nat;
    n:Nat;
    p:Pos;
    q:Pos;
    x:Int;
    y:Int;
    r:Real;
    s:Real;

eqn  ==(@cReal(x, p), @cReal(y, q))  =  ==(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p))));
     <(@cReal(x, p), @cReal(y, q))  =  <(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p))));
     <=(@cReal(x, p), @cReal(y, q))  =  <=(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p))));
     less_total(r,s) = <(r,s);
     Int2Real(x)  =  @cReal(x, @c1);
     Nat2Real(n)  =  @cReal(@cInt(n), @c1);
     Pos2Real(p)  =  @cReal(@cInt(@cNat(p)), @c1);
     Real2Int(@cReal(x, @c1))  =  x;
     Real2Nat(@cReal(x, @c1))  =  Int2Nat(x);
     Real2Pos(@cReal(x, @c1))  =  Int2Pos(x);
     min(r, s)  =  if(<(r, s), r, s);
     max(r, s)  =  if(<(r, s), s, r);
     abs(r)  =  if(<(r, @cReal(@cInt(@c0), @c1)), -(r), r);
     -(@cReal(x, p))  =  @cReal(-(x), p);
     succ(@cReal(x, p))  =  @cReal(+(x, @cInt(@cNat(p))), p);
     pred(@cReal(x, p))  =  @cReal(-(x, @cInt(@cNat(p))), p);
     +(@cReal(x, p), @cReal(y, q))  =  @redfrac(+(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p)))), @cInt(@cNat(*(p, q))));
     -(@cReal(x, p), @cReal(y, q))  =  @redfrac(-(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p)))), @cInt(@cNat(*(p, q))));
     *(@cReal(x, p), @cReal(y, q))  =  @redfrac(*(x, y), @cInt(@cNat(*(p, q))));
     *(r, @cReal(@cInt(@c0),p)) = @cReal(@cInt(@c0), @c1);
     *(@cReal(@cInt(@c0),p), r) = @cReal(@cInt(@c0), @c1);
     !=(y, @cInt(@c0))  ->  /(@cReal(x, p), @cReal(y, q))  =  @redfrac(*(x, @cInt(@cNat(q))), *(y, @cInt(@cNat(p))));
     /(p, q)  =  @redfrac(@cInt(@cNat(p)), @cInt(@cNat(q)));
     !=(n, @c0)  ->  /(m, n)  =  @redfrac(@cInt(m), @cInt(n));
     !=(y, @cInt(@c0))  ->  /(x, y)  =  @redfrac(x, y);
     exp(@cReal(x, p), @cInt(n))  =  @redfrac(exp(x, n), @cInt(@cNat(exp(p, n))));
     !=(x, @cInt(@c0))  ->  exp(@cReal(x, p), @cNeg(q))  =  @redfrac(@cInt(@cNat(exp(p, @cNat(q)))), exp(x, @cNat(q)));
     floor(@cReal(x, p))  =  div(x, p);
     ceil(r)  =  -(floor(-(r)));
     round(r)  =  floor(+(r, @cReal(@cInt(@cNat(@c1)), @cDub(false, @c1))));
     @redfrac(x, @cNeg(p))  =  @redfrac(-(x), @cInt(@cNat(p)));
     @redfrac(x, @cInt(@cNat(p)))  =  @redfracwhr(p, div(x, p), mod(x, p));
     @redfracwhr(p, x, @c0)  =  @cReal(x, @c1);
     @redfracwhr(p, x, @cNat(q))  =  @redfrachlp(@redfrac(@cInt(@cNat(p)), @cInt(@cNat(q))), x);
     @redfrachlp(@cReal(x, p), y)  =  @cReal(+(@cInt(@cNat(p)), *(y, x)), Int2Pos(x));
