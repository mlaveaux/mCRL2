% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the List data sort.

#using S
#include bool.spec
#include pos64.spec
#include nat64.spec

sort List(S) <"list">;

cons [] <"empty"> : List(S)                                            external defined_by_rewrite_rules;
     |> <"cons_"> : S <"left"> # List(S) <"right"> -> List(S)          external defined_by_rewrite_rules;

map in <"in"> : S <"left"> # List(S) <"right"> -> Bool                 external defined_by_rewrite_rules;
    # <"count"> : List(S) <"arg"> -> Nat                               external defined_by_rewrite_rules;
    <| <"snoc"> : List(S) <"left"> # S <"right"> -> List(S)            external defined_by_rewrite_rules;
    ++ <"concat"> : List(S) <"left"> # List(S) <"right"> -> List(S)    external defined_by_rewrite_rules;
    . <"element_at"> : List(S) <"left"> # Nat <"right"> -> S           external defined_by_rewrite_rules;
    head <"head"> : List(S) <"arg"> -> S                               external defined_by_rewrite_rules;
    tail <"tail"> : List(S) <"arg"> -> List(S)                         external defined_by_rewrite_rules;
    rhead <"rhead"> : List(S) <"arg"> -> S                             external defined_by_rewrite_rules;
    rtail <"rtail"> : List(S) <"arg"> -> List(S)                       external defined_by_rewrite_rules;

var d:S;
    e:S;
    s:List(S);
    t:List(S);
    p:Pos;
    n:Nat;
eqn ==([], |>(d,s)) = false;
    ==(|>(d,s), []) = false;
    ==(|>(d,s), |>(e,t)) = &&(==(d,e), ==(s,t));
    <([],|>(d,s)) = true;
    <(|>(d,s),[]) = false;
    <(|>(d,s), |>(e,t)) = ||(&&(==(d,e),<(s,t)),<(d,e));
    less_total([],[]) = false;
    less_total([],|>(d,s)) = true;
    less_total(|>(d,s), []) = false;
    less_total(|>(d,s), |>(e,t)) = ||(&&(==(d,e),less_total(s,t)),less_total(d,e));
    <=([],|>(d,s)) = true;
    <=(|>(d,s),[]) = false;
    <=(|>(d,s), |>(e,t)) = ||(&&(==(d,e),<=(s,t)),<(d,e));
    in(d,[]) = false;
    in(d,|>(e,s)) = ||(==(d,e), in(d,s));
    #([]) = @c0;
    #(|>(d,s)) = @succ_nat(#(s));
    <|([],d) = |>(d,[]);
    <|(|>(d,s), e) = |>(d, <|(s,e));
    ++([],s) = s;
    ++(|>(d,s), t) = |>(d, ++(s,t));
    ++(s,[]) = s;
    .(|>(d,s),n) = if(==(n,@most_significant_digitNat(@zero_word)),d,.(s, @natpred(n)));
    head(|>(d,s)) = d;
    tail(|>(d,s)) = s;
    rhead(|>(d,[])) = d;
    rhead(|>(d,|>(e,s))) = rhead(|>(e,s));
    rtail(|>(d,[])) = [];
    rtail(|>(d,|>(e,s))) = |>(d,rtail(|>(e,s)));

