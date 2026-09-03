% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the FBag data sort, denoting finite bags.
% Note that the specification relies on the underlying data type S to have a total ordering.
%
% The definition of an FBag originally had the shape
%
% sort FBag(S) <"fbag"> = struct {:} <"empty"> | @fbag_cons <"cons_"> : S <"arg1"> # Pos <"arg2"> # FBag(S) <"arg3">;
%
% However, this does not work as the automatically generated comparison operators <=, <, > and >= do not act
% as subset operators. Therefore, the constructors have been made explicit, as have the comparison operators.
% (April, 2017, Jan Friso Groote).
%
% Also changed @fbag_insert to become the constructor and @fbag_cons to become a map. All bags should have
% their elements in a list with @fbag_cons as head symbol be ordered. This is not the case in lists with @fbag_insert.
% If the @fbag_cons would be a constructor illegal lists would be constructed when evaluationg quantifications and
% sum operators. Now it is the case that too many bags will be generated when evaluating for instance a sum operator,
% but they are at least not incorrect.

#using S
#include bool.spec
#include pos1.spec
#include nat1.spec
#include fset1.spec

sort FBag(S) <"fbag">;

cons {:} <"empty"> : FBag(S)                                                               external defined_by_rewrite_rules;
     @fbag_insert <"insert"> : S <"arg1"> # Pos <"arg2"> # FBag(S) <"arg3"> -> FBag(S)     internal defined_by_rewrite_rules;

map @fbag_cons <"cons_"> : S <"arg1"> # Pos <"arg2"> # FBag(S) <"arg3"> -> FBag(S)         internal defined_by_rewrite_rules;
    @fbag_cinsert <"cinsert"> : S <"arg1"> # Nat <"arg2"> # FBag(S) <"arg3"> -> FBag(S)    internal defined_by_rewrite_rules;
    count <"count"> : S <"left"> # FBag(S) <"right"> -> Nat                                external defined_by_rewrite_rules;
    in <"in"> : S <"left"> # FBag(S) <"right"> -> Bool                                     external defined_by_rewrite_rules;
    @fset2fbag <"fset2fbag">: FSet(S) <"arg"> -> FBag(S)                                   internal defined_by_rewrite_rules;
    + <"union_"> : FBag(S) <"left"> # FBag(S) <"right"> -> FBag(S)                         external defined_by_rewrite_rules;
    * <"intersection"> : FBag(S) <"left"> # FBag(S) <"right"> -> FBag(S)                   external defined_by_rewrite_rules;
    - <"difference"> : FBag(S) <"left"> # FBag(S) <"right"> -> FBag(S)                     external defined_by_rewrite_rules;
    # <"count_all"> : FBag(S) <"arg"> -> Nat                                               external defined_by_rewrite_rules;
    pick <"pick">: FBag(S) <"arg"> -> S                                                    external defined_by_rewrite_rules;

var d: S;
    e: S;
    p: Pos;
    q: Pos;
    b: FBag(S);
    c: FBag(S);
    s: FSet(S);
eqn ==(@fbag_cons(d, p, b), {:})  =  false;
    ==({:}, @fbag_cons(d, p, b))  =  false;
    ==(@fbag_cons(d, p, b), @fbag_cons(e, q, c))  =  &&(==(p, q), &&(==(d,e),==(b, c)));
    <=(@fbag_cons(d, p, b), {:})  =  false;
    <=({:}, @fbag_cons(d, p, b))  =  true;
    <=(@fbag_cons(d, p, b), @fbag_cons(e, q, c))  =  if(<(d, e), false, if(==(d, e), &&(<=(p, q), <=(b, c)), <=(@fbag_cons(d, p, b), c)));
    <(@fbag_cons(d, p, b), {:})  =  false;
    <({:}, @fbag_cons(d, p, b))  =  true;
    <(@fbag_cons(d, p, b), @fbag_cons(e, q, c))  =  if(<(d, e), false, if(==(d, e), ||(&&(==(p, q), <(b,c)),  &&(<(p, q), <=(b, c))), <=(@fbag_cons(d, p, b), c)));
    less_total({:}, {:}) = false;
    less_total({:}, @fbag_cons(d, p, b))  =  true;
    less_total(@fbag_cons(d, p, b), {:})  =  false;
    less_total(@fbag_cons(d, p, b), @fbag_cons(e, q, c)) = ||(less_total(d, e), ||(&&(==(d, e),less_total(p, q)),&&(==(d, e),&&(==(p, q),less_total(b, c)))));
    @fbag_insert(d, p, {:})  =  @fbag_cons(d, p, {:});
    @fbag_insert(d, p, @fbag_cons(d, q, b))  =  @fbag_cons(d, @addc(false,p,q), b);
    less_total(d, e)  ->  @fbag_insert(d, p, @fbag_cons(e, q, b))  =  @fbag_cons(d, p, @fbag_cons(e, q, b));
    less_total(e, d)  ->  @fbag_insert(d, p, @fbag_cons(e, q, b))  =  @fbag_cons(e, q, @fbag_insert(d, p, b));
    @fbag_cinsert(d, @c0, b)  =  b;
    @fbag_cinsert(d, @cNat(p), b)  =  @fbag_insert(d, p, b);
    count(d, {:})  =  @c0;
    count(d, @fbag_cons(d, p, b))  =  @cNat(p);
    less_total(d, e)  ->  count(d, @fbag_cons(e, p, b))  =  @c0;
    less_total(e, d)  ->  count(d, @fbag_cons(e, p, b))  =  count(d, b);
    in(d, b)  =  >(count(d, b), @c0);
    @fset2fbag({})  =  {:};
    @fset2fbag(@fset_cons(d, s))  =  @fbag_cinsert(d, @cNat(@c1), @fset2fbag(s));
    -(b,{:}) = b;
    -({:},c) = {:};
    -(@fbag_cons(d,p,b),@fbag_cons(d,p,c)) = -(b,c);
    less_total(p,q) -> -(@fbag_cons(d,p,b),@fbag_cons(d,q,c)) = -(b,c);
    less_total(q,p) -> -(@fbag_cons(d,p,b),@fbag_cons(d,q,c)) = @fbag_cons(d,Nat2Pos(@gtesubtb(false,p,q)),-(b,c));
    less_total(d,e) -> -(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = @fbag_cons(d,p,-(b,@fbag_cons(e,q,c)));
    less_total(e,d) -> -(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = -(@fbag_cons(d,p,b),c);
    +(b,{:}) = b;
    +({:},c) = c;
    +(@fbag_cons(d,p,b),@fbag_cons(d,q,c)) = @fbag_cons(d,@addc(false,p,q),+(b,c));
    less_total(d,e) -> +(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = @fbag_cons(d,p,+(b,@fbag_cons(e,q,c)));
    less_total(e,d) -> +(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = @fbag_cons(e,q,+(@fbag_cons(d,p,b),c));
    *(b,{:}) = {:};
    *({:},c) = {:};
    *(@fbag_cons(d,p,b),@fbag_cons(d,q,c)) = @fbag_cons(d,min(p,q),*(b,c));
    less_total(d,e) -> *(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = *(b,@fbag_cons(e,q,c));
    less_total(e,d) -> *(@fbag_cons(d,p,b),@fbag_cons(e,q,c)) = *(@fbag_cons(d,p,b),c);
    #({:}) = @c0;
    #(@fbag_cons(d,p,{:})) = @cNat(p);
    #(@fbag_cons(d,p,@fbag_cons(e,q,b))) = @cNat(@addc(false,p,Nat2Pos(#(@fbag_cons(e,q,b)))));
    pick(@fbag_cons(d,p,b)) = d;
