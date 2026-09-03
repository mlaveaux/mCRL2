% Author(s): Aad Mathijssen, Jeroen Keiren
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Nat data sort.

#include bool.spec
#include pos1.spec
#supertypeof Pos

sort Nat <"nat">;
% Auxiliary sort natpair, pair of natural numbers
     @NatPair <"natpair">;

cons @c0 <"c0"> : Nat                                                                                         internal defined_by_rewrite_rules;
     @cNat <"cnat"> : Pos <"arg"> -> Nat                                                                      internal defined_by_rewrite_rules;
% Constructor for natpair
     @cPair <"cpair"> : Nat <"arg1"> # Nat <"arg2"> -> @NatPair                                               internal defined_by_rewrite_rules;

map Pos2Nat <"pos2nat"> : Pos <"arg"> -> Nat                                                                  external defined_by_rewrite_rules;
    Nat2Pos <"nat2pos"> : Nat <"arg"> -> Pos                                                                  external defined_by_rewrite_rules;
    max <"maximum">:Pos <"left"> #Nat <"right">->Pos                                                          external defined_by_rewrite_rules;
    max <"maximum">:Nat <"left"> #Pos <"right">->Pos                                                          external defined_by_rewrite_rules;
    max <"maximum">:Nat <"left"> #Nat <"right">->Nat                                                          external defined_by_rewrite_rules;
    min <"minimum">:Nat <"left"> #Nat <"right">->Nat                                                          external defined_by_rewrite_rules;
    succ <"succ">:Nat <"arg">->Pos                                                                            external defined_by_rewrite_rules;
    pred <"pred">:Pos <"arg">->Nat                                                                            external defined_by_rewrite_rules;
% Double and conditionally add 1
    @dub <"dub">:Bool <"left"> # Nat <"right"> -> Nat                                                         internal defined_by_rewrite_rules;
% Double and add 1
    @dubsucc <"dubsucc">:Nat <"arg"> -> Pos                                                                   internal defined_by_rewrite_rules;
    + <"plus">:Pos <"left"> #Nat <"right">->Pos                                                               external defined_by_rewrite_rules;
    + <"plus">:Nat <"left"> #Pos <"right">->Pos                                                               external defined_by_rewrite_rules;
    + <"plus">:Nat <"left"> #Nat <"right">->Nat                                                               external defined_by_rewrite_rules;
    @gtesubtb <"gte_subtract_with_borrow">: Bool <"arg1"> # Pos <"arg2"> # Pos <"arg3"> -> Nat                internal defined_by_rewrite_rules;
    * <"times">:Nat <"left"> #Nat <"right">->Nat                                                              external defined_by_rewrite_rules;
    div <"div">: Nat <"left"> # Pos <"right"> -> Nat                                                          external defined_by_rewrite_rules;
    mod <"mod">:Nat <"left"> # Pos <"right"> -> Nat                                                           external defined_by_rewrite_rules;
    exp <"exp">:Pos <"left"> # Nat <"right"> -> Pos                                                           external defined_by_rewrite_rules;
    exp <"exp">:Nat <"left"> # Nat <"right"> -> Nat                                                           external defined_by_rewrite_rules;
    @even <"even">:Nat <"arg"> -> Bool                                                                        internal defined_by_rewrite_rules;
    @monus <"monus">:Nat <"left"> # Nat <"right"> -> Nat                                                      internal defined_by_rewrite_rules;
    @swap_zero <"swap_zero">:Nat <"left"> # Nat <"right"> -> Nat                                              internal defined_by_rewrite_rules;
    @swap_zero_add <"swap_zero_add">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat         internal defined_by_rewrite_rules;
    @swap_zero_min <"swap_zero_min">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat         internal defined_by_rewrite_rules;
    @swap_zero_monus <"swap_zero_monus">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat     internal defined_by_rewrite_rules;
    sqrt <"sqrt">:Nat <"arg"> -> Nat                                                                          external defined_by_rewrite_rules;
    @sqrt_nat <"sqrt_nat_aux_func">:Nat <"arg1"> # Nat <"arg2"> # Pos <"arg3"> -> Nat                         internal defined_by_rewrite_rules;

% functions for natpair
    @first <"first"> : @NatPair <"arg"> -> Nat                                                                internal defined_by_rewrite_rules;
    @last <"last"> : @NatPair <"arg"> -> Nat                                                                  internal defined_by_rewrite_rules;
    @divmod <"divmod"> : Pos <"left"> # Pos <"right"> -> @NatPair                                             internal defined_by_rewrite_rules;
    @gdivmod <"generalised_divmod"> : @NatPair <"arg1"> # Bool <"arg2"> # Pos <"arg3"> -> @NatPair            internal defined_by_rewrite_rules;
    @ggdivmod <"doubly_generalised_divmod"> : Nat <"arg1"> # Nat <"arg2"> # Pos <"arg3"> -> @NatPair          internal defined_by_rewrite_rules;


var b:Bool;
    c:Bool;
    p:Pos;
    q:Pos;
    n:Nat;
    m:Nat;
    u:Nat;
    v:Nat;
eqn ==(@c0, @cNat(p)) = false;
    ==(@cNat(p), @c0) = false;
    ==(@cNat(p), @cNat(q)) = ==(p,q);
    <(n, @c0) = false;
    <(@c0, @cNat(p)) = true;
    <(@cNat(p), @cNat(q)) = <(p,q);
    less_total(m,n) = <(m,n);
    <=(@c0, n) = true;
    <=(@cNat(p), @c0) = false;
    <=(@cNat(p), @cNat(q)) = <=(p,q);
    Pos2Nat(p) = @cNat(p);
    Nat2Pos(@cNat(p)) = p;
    max(p, @c0) = p;
    max(p, @cNat(q)) = if(<=(p,q),q,p);
    max(@c0,p) = p;
    max(@cNat(p),q) = if(<=(p,q),q,p);
    max(m,n) = if(<=(m,n),n,m);
    min(m,n) = if(<=(m,n),m,n);
    succ(@c0) = @c1;
    succ(@cNat(p)) = succ(p);
% The rule below is essential for the enumeration of lists.
%    succ(succ(n)) = @cDub(!(@even(n)), succ(div(n,@cDub(false,@c1))));
    pred(@c1) = @c0;
% Pulling the @cNat outside of the if, instead of using pattern matching,
% helps the enumerator
    pred(@cDub(b,p)) = @cNat(if(b,@cDub(false,p),@dubsucc(pred(p))));
    @dubsucc(@c0) = @c1;
    @dubsucc(@cNat(p)) = @cDub(true,p);
    @dub(false,@c0) = @c0;
    @dub(true,@c0) = @cNat(@c1);
    @dub(b,@cNat(p)) = @cNat(@cDub(b,p));
    +(p,@c0) = p;
    +(p,@cNat(q)) = @addc(false,p,q);
    +(@c0,p) = p;
    +(@cNat(p),q) = @addc(false,p,q);
    +(@c0,n) = n;
    +(n,@c0) = n;
    +(@cNat(p),@cNat(q)) = @cNat(@addc(false,p,q));
    @gtesubtb(false,p,@c1) = pred(p);
    @gtesubtb(true,p,@c1) = pred(Nat2Pos(pred(p)));
    @gtesubtb(b,@cDub(c,p),@cDub(c,q)) = @dub(b, @gtesubtb(b,p,q));
    @gtesubtb(b,@cDub(false,p),@cDub(true,q)) = @dub(!(b),@gtesubtb(true,p,q));
    @gtesubtb(b,@cDub(true,p),@cDub(false,q)) = @dub(!(b),@gtesubtb(false,p,q));
    *(@c0,n) = @c0;
    *(n,@c0) = @c0;
    *(@cNat(p),@cNat(q)) = @cNat(*(p,q));
    exp(p,@c0) = @c1;
    exp(p,@cNat(@c1)) = p;
    exp(p,@cNat(@cDub(false,q))) = exp(*(p,p),@cNat(q));
    exp(p,@cNat(@cDub(true,q))) = *(p,exp(*(p,p),@cNat(q)));
    exp(n,@c0) = @cNat(@c1);
    exp(@c0, @cNat(p)) = @c0;
    exp(@cNat(p),n) = @cNat(exp(p,n));
    @even(@c0) = true;
    @even(@cNat(@c1)) = false;
    @even(@cNat(@cDub(b,p))) = !(b);
    div(@c0,p) = @c0;
    div(@cNat(p),q) = @first(@divmod(p,q));
    mod(@c0,p) = @c0;
    mod(@cNat(p),q) = @last(@divmod(p,q));
    @monus(@c0,n) = @c0;
    @monus(n,@c0) = n;
    @monus(@cNat(p),@cNat(q)) = @gtesubtb(false,p,q);
    @swap_zero(m,@c0) = m;
    @swap_zero(@c0,n) = n;
    @swap_zero(@cNat(p),@cNat(p)) = @c0;
    !=(p,q) -> @swap_zero(@cNat(p),@cNat(q)) = @cNat(q);
    @swap_zero_add(@c0, @c0, m, n) = +(m,n);
    @swap_zero_add(@cNat(p), @c0, m, @c0) = m;
    @swap_zero_add(@cNat(p), @c0, m, @cNat(q)) = @swap_zero(@cNat(p), +(@swap_zero(@cNat(p),m), @cNat(q)));
    @swap_zero_add(@c0, @cNat(p), @c0, n) = n;
    @swap_zero_add(@c0, @cNat(p), @cNat(q), n) = @swap_zero(@cNat(p), +(@cNat(q), @swap_zero(@cNat(p), n)));
    @swap_zero_add(@cNat(p), @cNat(q), m, n) = @swap_zero(+(@cNat(p), @cNat(q)), +(@swap_zero(@cNat(p),m),@swap_zero(@cNat(q),n)));
    @swap_zero_min(@c0, @c0, m, n) = min(m,n);
    @swap_zero_min(@cNat(p), @c0, m, @c0) = @c0;
    @swap_zero_min(@cNat(p), @c0, m, @cNat(q)) = min(@swap_zero(@cNat(p),m), @cNat(q));
    @swap_zero_min(@c0, @cNat(p), @c0, n) = @c0;
    @swap_zero_min(@c0, @cNat(p), @cNat(q), n) = min(@cNat(q), @swap_zero(@cNat(p),n));
    @swap_zero_min(@cNat(p), @cNat(q), m, n) = @swap_zero(min(@cNat(p), @cNat(q)), min(@swap_zero(@cNat(p), m), @swap_zero(@cNat(q), n)));
    @swap_zero_monus(@c0, @c0, m, n) = @monus(m,n);
    @swap_zero_monus(@cNat(p), @c0, m, @c0) = m;
    @swap_zero_monus(@cNat(p), @c0, m, @cNat(q)) = @swap_zero(@cNat(p), @monus(@swap_zero(@cNat(p), m), @cNat(q)));
    @swap_zero_monus(@c0, @cNat(p), @c0, n) = @c0;
    @swap_zero_monus(@c0, @cNat(p), @cNat(q), n) = @monus(@cNat(q), @swap_zero(@cNat(p), n));
    @swap_zero_monus(@cNat(p), @cNat(q), m, n) = @swap_zero(@monus(@cNat(p),@cNat(q)),@monus(@swap_zero(@cNat(p),m), @swap_zero(@cNat(q),n)));
    sqrt(@c0) = @c0;
    sqrt(@cNat(p)) = @sqrt_nat(@cNat(p),@c0,@powerlog2(p));
    @sqrt_nat(n,m,@c1) = if(<=(n,m),@c0,@cNat(@c1));
    @sqrt_nat(n,m,@cDub(b,p)) =
              if(>(*(+(@cNat(@cDub(b,p)),m),@cNat(@cDub(b,p))),n),
                   @sqrt_nat(n,m,p),
                   +(@cNat(@cDub(b,p)),@sqrt_nat(@monus(n,*(+(@cNat(@cDub(b,p)),m),@cNat(@cDub(b,p)))),+(m,@cNat(@cDub(false,@cDub(b,p)))),p)));

% equations for natpair
    ==(@cPair(m,n), @cPair(u,v)) = &&(==(m,u),==(n,v));
    <(@cPair(m,n), @cPair(u,v)) = ||(<(m,u), &&(==(m,u), <(n,v)));
    <=(@cPair(m,n), @cPair(u,v)) = ||(<(m,u), &&(==(m,u), <=(n,v)));
    @first(@cPair(m,n)) = m;
    @last(@cPair(m,n)) = n;
    @divmod(@c1,@c1) = @cPair(@cNat(@c1),@c0);
    @divmod(@c1,@cDub(b,p)) = @cPair(@c0, @cNat(@c1));
    @divmod(@cDub(b,p),q) = @gdivmod(@divmod(p,q),b,q);
    @gdivmod(@cPair(m,n),b,p) = @ggdivmod(@dub(b,n),m,p);
    @ggdivmod(@c0,n,p) = @cPair(@dub(false,n),@c0);
    <(p,q) -> @ggdivmod(@cNat(p),n,q) = @cPair(@dub(false,n),@cNat(p));
    <=(q,p) -> @ggdivmod(@cNat(p),n,q) = @cPair(@dub(true,n),@gtesubtb(false,p,q));
