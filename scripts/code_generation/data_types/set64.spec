% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Set data sort.

#using S
#include bool.spec
#include fset64.spec
#supertypeof FSet

sort Set(S) <"set_">;

cons @set <"constructor"> : (S -> Bool) <"left"> # FSet(S) <"right"> -> Set(S)               internal defined_by_rewrite_rules;
% map {} <"empty"> : Set(S); Move this to FSet(S);
% I think that @setfset and @setcomp should not be part of the rewrite system, but
% become part of the internal generation of set representations. JFG
map @setfset <"set_fset"> : FSet(S) <"arg"> -> Set(S)                                        internal defined_by_rewrite_rules;
    @setcomp <"set_comprehension"> : (S -> Bool) <"arg"> -> Set(S)                           internal defined_by_rewrite_rules;
    in <"in"> : S <"left"> # Set(S) <"right"> -> Bool                                        external defined_by_rewrite_rules;
    ! <"complement"> : Set(S) <"arg"> -> Set(S)                                              external defined_by_rewrite_rules;
    + <"union_"> : Set(S) <"left"> # Set(S) <"right"> -> Set(S)                              external defined_by_rewrite_rules;
    * <"intersection"> : Set(S) <"left"> # Set(S) <"right"> -> Set(S)                        external defined_by_rewrite_rules;
    * <"intersection"> : FSet(S) <"left"> # Set(S) <"right"> -> FSet(S)                      external defined_by_rewrite_rules;
    * <"intersection"> : Set(S) <"left"> # FSet(S) <"right"> -> FSet(S)                      external defined_by_rewrite_rules;
    - <"difference"> : Set(S) <"left"> # Set(S) <"right"> -> Set(S)                          external defined_by_rewrite_rules;
    - <"difference"> : FSet(S) <"left"> # Set(S) <"right"> -> FSet(S)                        external defined_by_rewrite_rules;
    @false_ <"false_function"> : S <"arg"> -> Bool                                           internal defined_by_rewrite_rules;
    @true_ <"true_function"> : S <"arg"> -> Bool                                             internal defined_by_rewrite_rules;
    @not_ <"not_function"> : (S -> Bool) <"arg"> -> S -> Bool                                internal defined_by_rewrite_rules;
    @and_ <"and_function"> : (S -> Bool) <"left"> # (S -> Bool) <"right"> -> S -> Bool       internal defined_by_rewrite_rules;
    @or_ <"or_function"> : (S -> Bool) <"left"> # (S -> Bool) <"right"> -> S -> Bool         internal defined_by_rewrite_rules;
    @fset_union <"fset_union"> : (S -> Bool) <"arg1"> # (S -> Bool) <"arg2"> # FSet(S) <"arg3"> # FSet(S) <"arg4"> -> FSet(S)              
                                                                                             internal defined_by_rewrite_rules;
    @fset_inter <"fset_intersection">: (S -> Bool) <"arg1"> # (S -> Bool) <"arg2"> # FSet(S) <"arg3"> # FSet(S) <"arg4"> -> FSet(S)        
                                                                                             internal defined_by_rewrite_rules;

var d : S;
    e : S;
    s : FSet(S);
    t : FSet(S);
    f : S->Bool;
    g : S->Bool;
    x : Set(S);
    y : Set(S);
% eqn {}  =  @set(@false_, {});
eqn @setfset(s)  =  @set(@false_, s);
    @setcomp(f)  =  @set(f, {});
    in(e, @set(f, s))  =  !=(f(e), in(e, s));
    ==(@set(f, s), @set(g, t))  =  forall(c:S, ==(==(f(c),g(c)),==(in(c,s),in(c,t))));
    <(x, y)  =  &&(<=(x, y), !=(x, y));
    <=(x,y) = ==(*(x,y),x);
    less_total(@set(f, s), @set(g, t)) = ||(less_total(s, t), &&(==(s, t), less_total(f, s)));
    !(@set(f, s))  =  @set(@not_(f), s);
    +(x, x) = x;
    +(x, +(x, y)) = +(x, y);
    +(x, +(y, x)) = +(y, x);
    +(+(x, y), x) = +(x, y);
    +(+(y, x), x) = +(y, x);
    +(@set(f, s), @set(g, t))  =  @set(@or_(f, g), @fset_union(f, g, s, t));
    *(x, x) = x;
    *(x, *(x, y)) = *(x, y);
    *(x, *(y, x)) = *(y, x);
    *(*(x, y), x) = *(x, y);
    *(*(y, x), x) = *(y, x);
    *(@set(f, s), @set(g, t))  =  @set(@and_(f, g), @fset_inter(f, g, s, t));
    *({},x) = {};
    *(@fset_cons(d, s), x) = if(in(d, x), @fset_cons(d, *(s, x)), *(s, x));
    *(x,s) = *(s,x);
    -(x, y)  =  *(x, !(y));
    -(s, x) = *(s, !(x));
    @false_(e)  =  false;
    @true_(e)  =  true;
    ==(@false_, @true_)  =  false;
    ==(@true_, @false_)  =  false;
    @not_(f)(e)  =  !(f(e));
    @not_(@false_)  =  @true_;
    @not_(@true_)  =  @false_;
    @and_(f, g)(e)  =  &&(f(e), g(e));
    @and_(f, f)  =  f;
    @and_(f, @false_)  =  @false_;
    @and_(@false_, f)  =  @false_;
    @and_(f, @true_)  =  f;
    @and_(@true_, f)  =  f;
    @or_(f, f)  =  f;
    @or_(f, @false_)  =  f;
    @or_(@false_, f)  =  f;
    @or_(f, @true_)  =  @true_;
    @or_(@true_, f)  =  @true_;
    @or_(f, g)(e)  =  ||(f(e), g(e));
    @fset_union(@false_, @false_, s, t)  = +(s,t);
    @fset_union(f, g, {}, {})  =  {};
    @fset_union(f, g, @fset_cons(d, s), {})  =  @fset_cinsert(d, !(g(d)), @fset_union(f, g, s, {}));
    @fset_union(f, g, {}, @fset_cons(e, t))  =  @fset_cinsert(e, !(f(e)), @fset_union(f, g, {}, t));
    @fset_union(f, g, @fset_cons(d, s), @fset_cons(d, t))  =  @fset_cinsert(d, ==(f(d), g(d)), @fset_union(f, g, s, t));
    less_total(d, e)  ->  @fset_union(f, g, @fset_cons(d, s), @fset_cons(e, t))  =  @fset_cinsert(d, !(g(d)), @fset_union(f, g, s, @fset_cons(e, t)));
    less_total(e, d)  ->  @fset_union(f, g, @fset_cons(d, s), @fset_cons(e, t))  =  @fset_cinsert(e, !(f(e)), @fset_union(f, g, @fset_cons(d, s), t));
    @fset_inter(f, g, {}, {})  =  {};
    @fset_inter(f, g, @fset_cons(d, s), {})  =  @fset_cinsert(d, g(d), @fset_inter(f, g, s, {}));
    @fset_inter(f, g, {}, @fset_cons(e, t))  =  @fset_cinsert(e, f(e), @fset_inter(f, g, {}, t));
    @fset_inter(f, g, @fset_cons(d, s), @fset_cons(d, t))  =  @fset_cinsert(d, ==(f(d), g(d)), @fset_inter(f, g, s, t));
    less_total(d, e)  ->  @fset_inter(f, g, @fset_cons(d, s), @fset_cons(e, t))  =  @fset_cinsert(d, g(d), @fset_inter(f, g, s, @fset_cons(e, t)));
    less_total(e, d)  ->  @fset_inter(f, g, @fset_cons(d, s), @fset_cons(e, t))  =  @fset_cinsert(e, f(e), @fset_inter(f, g, @fset_cons(d, s), t));
