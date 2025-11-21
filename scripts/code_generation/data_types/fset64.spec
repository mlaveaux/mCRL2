% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the FSet data sort, denoting finite sets.
% Note that the data type relies on the underlying data type S to have a total ordering.
%
% FSet(S) was initially specified using:
%
% sort FSet(S) <"fset"> = struct {} <"empty"> | @fset_cons <"cons_"> : S <"left"> # FSet(S) <"right">;
%
% But this does not work as the automatically generated relation <= is not the subset relation as
% would be expected. The same holds for all other comparison operators. Therefore an explicit definition
% is put into place (Jan Friso Groote, April 2017; problem reported by Tim Willemse).
%
% Also changed @fset_insert and @fset_cons. @fset_insert generates unordered lists. @fset_cons is used
% for ordered lists only. This has especially an effect when generating finite sets using quantifiers
% and sums. When @fset_cons is a constructor unordered lists are generated, for which the rewrite rules
% below do not work properly.

#using S
#include bool.spec
#include nat64.spec

sort FSet(S) <"fset">;

cons {} <"empty"> : FSet(S)                                                                                                                 external defined_by_rewrite_rules;
     @fset_insert <"insert">: S <"left"> # FSet(S) <"right"> -> FSet(S)                                                                     internal defined_by_rewrite_rules;

map  @fset_cons <"cons_"> : S <"left"> # FSet(S) <"right"> -> FSet(S)                                                                       internal defined_by_rewrite_rules;
     @fset_cinsert <"cinsert">: S <"arg1"> # Bool <"arg2"> # FSet(S) <"arg3"> -> FSet(S)                                                    internal defined_by_rewrite_rules;
     in <"in">: S <"left"> # FSet(S) <"right"> -> Bool                                                                                      external defined_by_rewrite_rules;
%     @fset_union <"fset_union"> : (S -> Bool) <"arg1"> # (S -> Bool) <"arg2"> # FSet(S) <"arg3"> # FSet(S) <"arg4"> -> FSet(S)             internal defined_by_rewrite_rules;
%     @fset_inter <"fset_intersection">: (S -> Bool) <"arg1"> # (S -> Bool) <"arg2"> # FSet(S) <"arg3"> # FSet(S) <"arg4"> -> FSet(S)       internal defined_by_rewrite_rules;
     - <"difference">: FSet(S) <"left"> # FSet(S) <"right"> -> FSet(S)                                                                      external defined_by_rewrite_rules;
     + <"union_"> : FSet(S) <"left"> # FSet(S) <"right"> -> FSet(S)                                                                         external defined_by_rewrite_rules;
     * <"intersection"> : FSet(S) <"left"> # FSet(S) <"right"> -> FSet(S)                                                                   external defined_by_rewrite_rules;
     # <"count"> : FSet(S) <"arg"> -> Nat                                                                                                   external defined_by_rewrite_rules;
     pick <"pick">: FSet(S) <"arg"> -> S                                                                                                    external defined_by_rewrite_rules;

var d:S;
    e:S;
    f:S->Bool;
    g:S->Bool;
    s:FSet(S);
    t:FSet(S);
eqn ==({}, @fset_cons(d, s))  =  false;
    ==(@fset_cons(d, s), {})  =  false;
    ==(@fset_cons(d, s), @fset_cons(e, t))  =  &&(==(d, e), ==(s, t));
    <=({}, @fset_cons(d, s))  =  true;
    <=(@fset_cons(d, s), {})  =  false;
    <=(@fset_cons(d, s), @fset_cons(e, t))  =  if(<(d,e), false, if(==(d, e), <=(s, t), <=(@fset_cons(d, s), t)));
    <({}, @fset_cons(d, s))  =  true;
    <(@fset_cons(d, s), {})  =  false;
    <(@fset_cons(d, s), @fset_cons(e, t))  =  if(<(d,e), false, if(==(d, e), <(s, t), <=(@fset_cons(d, s), t)));
    less_total({}, {}) = false;
    less_total({}, @fset_cons(d, s)) = true;
    less_total(@fset_cons(d, s), @fset_cons(e, t)) = ||(less_total(d, e), &&(==(d, e), less_total(s, t)));
    @fset_insert(d, {})  =  @fset_cons(d, {});
    @fset_insert(d, @fset_cons(d, s))  =  @fset_cons(d, s);
    less_total(d, e)  ->  @fset_insert(d, @fset_cons(e, s))  =  @fset_cons(d, @fset_cons(e, s));
    less_total(e, d)  ->  @fset_insert(d, @fset_cons(e, s))  =  @fset_cons(e, @fset_insert(d, s));
    @fset_cinsert(d, false, s)  =  s;
    @fset_cinsert(d, true, s)  =  @fset_insert(d, s);
    in(d, {})  =  false;
    in(d,@fset_cons(e,s)) = ||(==(d,e),in(d,s));
% The rule below is added such that set membership can still be calculated although the set elements cannot be effectively ordered.
    in(d,@fset_insert(e,s)) = ||(==(d,e),in(d,s));
    -(s,{}) = s;
    -({},t) = {};
    -(@fset_cons(d,s),@fset_cons(d,t)) = -(s,t);
    less_total(d,e) -> -(@fset_cons(d,s),@fset_cons(e,t)) = @fset_cons(d,-(s,@fset_cons(e,t)));
    less_total(e,d) -> -(@fset_cons(d,s),@fset_cons(e,t)) = -(@fset_cons(d,s),t);
    +(s,{}) = s;
    +({},t) = t;
    +(@fset_cons(d,s),@fset_cons(d,t)) = @fset_cons(d,+(s,t));
    less_total(d,e) -> +(@fset_cons(d,s),@fset_cons(e,t)) = @fset_cons(d,+(s,@fset_cons(e,t)));
    less_total(e,d) -> +(@fset_cons(d,s),@fset_cons(e,t)) = @fset_cons(e,+(@fset_cons(d,s),t));
    *(s,{}) = {};
    *({},t) = {};
    *(@fset_cons(d,s),@fset_cons(d,t)) = @fset_cons(d,*(s,t));
    less_total(d,e) -> *(@fset_cons(d,s),@fset_cons(e,t)) = *(s,@fset_cons(e,t));
    less_total(e,d) -> *(@fset_cons(d,s),@fset_cons(e,t)) = *(@fset_cons(d,s),t);
    #({}) = @c0;
    #(@fset_cons(d,s)) = @succ_nat(#(s));
% It is odd that the rule below has to be added separately.
    !=(s,t) = !(==(s,t));
    pick(@fset_cons(d,s)) = d;
