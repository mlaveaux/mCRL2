% Author(s): Jan Friso Groote 
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Nat data sort.

#include machine_word.spec
#include bool.spec
#include pos64.spec
#supertypeof Pos

sort Nat <"nat">;
% Auxiliary sort natpair, pair of natural numbers
     @NatNatPair <"natnatpair">;

cons @c0 <"c0"> : Nat                                                                                            internal defined_by_rewrite_rules;
     @succ_nat <"succ_nat"> : Nat <"arg">->Nat                                                                   internal defined_by_rewrite_rules;
% Is the constructor below needed?
% Constructor for natpair
     @nnPair <"nnpair"> : Nat <"arg1"> # Nat <"arg2"> -> @NatNatPair                                             internal defined_by_rewrite_rules;

map  @most_significant_digitNat <"most_significant_digit_nat">: @word <"arg"> -> Nat                             internal defined_by_rewrite_rules;   
% concat_digit(p,w) represents (2^N)*p + w.
     @concat_digit <"concat_digit"> : Nat <"arg1"> # @word <"arg2"> -> Nat                                       internal defined_by_rewrite_rules;
     @equals_zero <"equals_zero">: Nat <"arg"> -> Bool                                                           internal defined_by_rewrite_rules;
     @not_equals_zero <"not_equals_zero">: Nat <"arg"> -> Bool                                                   internal defined_by_rewrite_rules;
     @equals_one <"equals_one">: Nat <"arg"> -> Bool                                                             internal defined_by_rewrite_rules;
     Pos2Nat <"pos2nat"> : Pos <"arg"> -> Nat                                                                    external defined_by_rewrite_rules;
     Nat2Pos <"nat2pos"> : Nat <"arg"> -> Pos                                                                    external defined_by_rewrite_rules;
     succ <"succ"> : Nat <"arg">->Pos                                                                            external defined_by_rewrite_rules;
     max <"maximum">:Pos <"left"> # Nat <"right">->Pos                                                           external defined_by_rewrite_rules;
     max <"maximum">:Nat <"left"> # Pos <"right">->Pos                                                           external defined_by_rewrite_rules;
     max <"maximum">:Nat <"left"> # Nat <"right">->Nat                                                           external defined_by_rewrite_rules;
     min <"minimum">:Nat <"left"> # Nat <"right">->Nat                                                           external defined_by_rewrite_rules;
     pred <"pred">:Pos <"arg">->Nat                                                                              external defined_by_rewrite_rules;
     @pred_whr <"pred_whr">:Nat <"arg">->Nat                                                                     internal defined_by_rewrite_rules;
     + <"plus">:Pos <"left"> # Nat <"right">->Pos                                                                external defined_by_rewrite_rules;
     + <"plus">:Nat <"left"> # Pos <"right">->Pos                                                                external defined_by_rewrite_rules;
     + <"plus">:Nat <"left"> # Nat <"right">->Nat                                                                external defined_by_rewrite_rules;
     @add_with_carry <"add_with_carry">:Nat <"left"> #Nat <"right">->Nat                                         external defined_by_rewrite_rules;
% The following function is used when the symbol + is overloaded, such as in fbags.
     @plus_nat <"auxiliary_plus_nat">: Nat <"left"> # Nat <"right"> -> Nat                                       external defined_by_rewrite_rules;
     * <"times">:Nat <"left"> # Nat <"right">->Nat                                                               external defined_by_rewrite_rules;
     @times_ordered <"times_ordered">:Nat <"left"> # Nat <"right">->Nat                                          internal defined_by_rewrite_rules;
     @times_overflow <"times_overflow">: Nat <"arg1"> # @word <"arg2"> # @word <"arg3"> -> Nat                   external defined_by_rewrite_rules;
     div <"div">: Nat <"left"> # Pos <"right"> -> Nat                                                            external defined_by_rewrite_rules;
     mod <"mod">:Nat <"left"> # Pos <"right"> -> Nat                                                             external defined_by_rewrite_rules;
     exp <"exp">:Pos <"left"> # Nat <"right"> -> Pos                                                             external defined_by_rewrite_rules;
     exp <"exp">:Nat <"left"> # Nat <"right"> -> Nat                                                             external defined_by_rewrite_rules;
     sqrt <"sqrt">:Nat <"arg"> -> Nat                                                                            external defined_by_rewrite_rules;
     @natpred <"natpred">: Nat <"arg"> -> Nat                                                                    internal defined_by_rewrite_rules;
     @is_odd <"is_odd">: Nat <"arg"> -> Bool                                                                     internal defined_by_rewrite_rules;
     @div2 <"div2">: Nat <"arg"> -> Nat                                                                          internal defined_by_rewrite_rules;
     @monus <"monus">:Nat <"left"> # Nat <"right"> -> Nat                                                        internal defined_by_rewrite_rules;
     @monus_whr <"monus_whr">:Nat <"arg1"> # @word <"arg2"> # Nat <"arg3"> # @word <"arg4"> # Nat <"arg5">-> Nat internal defined_by_rewrite_rules;
     @exp_aux3p <"exp_aux3p">: Bool <"arg1"> # Pos <"arg2"> # @word <"arg3"> -> Pos                              internal defined_by_rewrite_rules;
     @exp_aux4p <"exp_aux4p">: Bool <"arg1"> # Pos <"arg2"> # Nat <"arg3"> # @word <"arg4"> -> Pos               internal defined_by_rewrite_rules;
     @exp_aux3n <"exp_aux3n">: Bool <"arg1"> # Nat <"arg2"> # @word <"arg3"> -> Nat                              internal defined_by_rewrite_rules;
     @exp_aux4n <"exp_aux4n">: Bool <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # @word <"arg4"> -> Nat               internal defined_by_rewrite_rules;
     @exp_auxtruep <"exp_auxtruep">: Pos <"arg1"> # Nat <"arg2"> # @word <"arg3"> -> Nat                         internal defined_by_rewrite_rules;
     @exp_auxtruen <"exp_auxtruen">: Nat <"arg1"> # Nat <"arg2"> # @word <"arg3"> -> Nat                         internal defined_by_rewrite_rules;
     @exp_auxfalsep <"exp_auxfalsep">: Pos <"arg1"> # Nat <"arg2"> # @word <"arg3"> -> Nat                       internal defined_by_rewrite_rules;
     @exp_auxfalsen <"exp_auxfalsen">: Nat <"arg1"> # Nat <"arg2"> # @word <"arg3"> -> Nat                       internal defined_by_rewrite_rules;
     @div_bold <"div_bold">: Nat <"arg1"> # Pos <"arg2"> -> @word                                                internal defined_by_rewrite_rules;
     @div_bold_whr <"div_bold_whr">: Nat <"arg1"> # @word <"arg2"> # Pos <"arg3"> # @word <"arg4"> # @word <"arg5"> # @word <"arg6"> -> @word   internal defined_by_rewrite_rules;
     @div_whr1 <"div_whr1">: Nat <"arg1"> # @word <"arg2"> # @word <"arg3"> # @NatNatPair <"arg4"> -> Nat        internal defined_by_rewrite_rules;
     @div_whr2 <"div_whr2">: Nat <"arg1"> # @word <"arg2"> # Pos <"arg3"> # @word <"arg4"> # @NatNatPair <"arg5"> -> Nat                        internal defined_by_rewrite_rules;
     @mod_whr1 <"mod_whr1">: @word <"arg1"> # Pos <"arg2"> # @word <"arg3"> # Nat <"arg4"> -> Nat               internal defined_by_rewrite_rules; 
     @divmod_aux <"divmod_aux">: Nat <"arg1"> # Pos <"arg2"> -> @NatNatPair                                      internal defined_by_rewrite_rules;
     @divmod_aux_whr1 <"divmod_aux_whr1">: Nat <"arg1"> # @word <"arg2"> # @word <"arg3"> # @NatNatPair <"arg4"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @divmod_aux_whr2 <"divmod_aux_whr2">: Nat <"arg1"> # @word <"arg2"> # Pos <"arg3"> # @word <"arg4"> # Nat <"arg5"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @divmod_aux_whr3 <"divmod_aux_whr3">: Nat <"arg1"> # @word <"arg2"> # Pos <"arg3"> # @word <"arg4"> # Nat <"arg5"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @divmod_aux_whr4 <"divmod_aux_whr4">: @word <"arg1"> # Pos <"arg2"> # @word <"arg3"> # @NatNatPair <"arg4"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @divmod_aux_whr5 <"divmod_aux_whr5">: Pos <"arg1"> # @word <"arg2"> # @NatNatPair <"arg3"> # Nat <"arg4"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @divmod_aux_whr6 <"divmod_aux_whr6">: Pos <"arg1"> # @word <"arg2"> # @NatNatPair <"arg3"> # Nat <"arg4"> # Nat <"arg5"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @msd <"msd">:  Nat <"arg"> -> @word                                                                         internal defined_by_rewrite_rules;
     @swap_zero <"swap_zero">:Nat <"left"> # Nat <"right"> -> Nat                                                internal defined_by_rewrite_rules;
     @swap_zero_add <"swap_zero_add">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat           internal defined_by_rewrite_rules;
     @swap_zero_min <"swap_zero_min">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat           internal defined_by_rewrite_rules;
     @swap_zero_monus <"swap_zero_monus">:Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> -> Nat       internal defined_by_rewrite_rules;
     @sqrt_whr1 <"sqrt_whr1">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> -> Nat          internal defined_by_rewrite_rules;
     @sqrt_whr2 <"sqrt_whr2">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> # @word <"arg5"> -> Nat          internal defined_by_rewrite_rules;
     @sqrt_pair <"sqrt_pair">: Nat <"arg"> -> @NatNatPair                                                        internal defined_by_rewrite_rules;
     @sqrt_pair_whr1 <"sqrt_pair_whr1">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # Nat <"arg4"> -> @NatNatPair  internal defined_by_rewrite_rules;
     @sqrt_pair_whr2 <"sqrt_pair_whr2">: @word <"arg1"> # @word <"arg2"> # @word <"arg3"> # @word <"arg4"> # Nat <"arg5"> -> @NatNatPair          internal defined_by_rewrite_rules;
     @sqrt_pair_whr3 <"sqrt_pair_whr3">: @word <"arg1"> # @word <"arg2"> # @NatNatPair <"arg3"> -> @NatNatPair   internal defined_by_rewrite_rules;
     @sqrt_pair_whr4 <"sqrt_pair_whr4">: Nat <"arg1"> # @word <"arg2"> # @NatNatPair <"arg3"> # Nat <"arg4"> # Nat <"arg5"> # Nat <"arg6"> -> @NatNatPair          internal defined_by_rewrite_rules;
     @sqrt_pair_whr5 <"sqrt_pair_whr5">: @NatNatPair <"arg1"> # Nat <"arg2"> # Nat <"arg3"> # Nat <"arg4"> # Nat <"arg5"> -> @NatNatPair          internal defined_by_rewrite_rules;
     @sqrt_pair_whr6 <"sqrt_pair_whr6">: Nat <"arg1"> # Nat <"arg2"> # Nat <"arg3"> -> @NatNatPair               internal defined_by_rewrite_rules;
% functions for pairs.
     @first <"first"> : @NatNatPair <"arg"> -> Nat                                                               internal defined_by_rewrite_rules;
     @last <"last"> : @NatNatPair <"arg"> -> Nat                                                                 internal defined_by_rewrite_rules;


var  b:Bool;
     p:Pos;
     p1:Pos;
     p2:Pos;
     n:Nat;
     n1:Nat;
     n2:Nat;
     m:Nat;
     m1:Nat;
     m2:Nat;
     m3:Nat;
     m4:Nat;
     m5:Nat;
     predp:Nat;
     diff:Nat;
     shift_n1:Nat;
     solution:Nat;
     pq:Nat;
     y:Nat;
     y_guess:Nat;
     pair_:@NatNatPair;
     lp:Nat;
     w:@word;
     w1:@word;
     w2:@word;
     w3:@word;
     w4:@word;
     shift_w:@word;
     overflow:@word;

eqn  @c0 = @most_significant_digitNat(@zero_word);
     @equals_zero(@most_significant_digitNat(w)) = @equals_zero_word(w);
     @equals_zero(@concat_digit(n,w)) = false;
     @equals_zero(@succ_nat(n)) = false;
     @not_equals_zero(@most_significant_digitNat(w)) = @not_equals_zero_word(w);
     @not_equals_zero(@concat_digit(n,w)) = true;
     @not_equals_zero(@succ_nat(n)) = true;
     @equals_one(@most_significant_digitNat(w)) = @equals_one_word(w);
     @equals_one(@concat_digit(n,w)) = false;
     @succ_nat(@most_significant_digitNat(w)) = if(@equals_max_word(w),
                                                    @concat_digit(@most_significant_digitNat(@one_word),@zero_word),
                                                    @most_significant_digitNat(@succ_word(w)));
     @succ_nat(@concat_digit(n,w)) = if(@equals_max_word(w),
                                           @concat_digit(@succ_nat(n),@zero_word),
                                           @concat_digit(n,@succ_word(w)));
 
     succ(@most_significant_digitNat(w)) = if(@equals_max_word(w),
                                                    @concat_digit(@most_significant_digit(@one_word),@zero_word),
                                                    @most_significant_digit(@succ_word(w)));
     succ(@concat_digit(n,w)) = if(@equals_max_word(w),
                                           @concat_digit(succ(n),@zero_word),
                                           @concat_digit(Nat2Pos(n),@succ_word(w)));
 
% The rules for comparison operators in conjunction with @succ_nat are required in enumerations. 
     ==(@most_significant_digitNat(w1), @most_significant_digitNat(w2)) = @equal(w1,w2);
     ==(@concat_digit(n,w1), @most_significant_digitNat(w2)) = false;
     ==(@most_significant_digitNat(w1),@concat_digit(n,w2)) = false;
     ==(@concat_digit(n1,w1), @concat_digit(n2,w2)) = &&(@equal(w1,w2), ==(n1,n2));
     (@equals_zero(n2)) -> ==(@succ_nat(n1),n2) = false;
     (@not_equals_zero(n2)) -> ==(@succ_nat(n1),n2) = ==(n1,@natpred(n2));
     (@equals_zero(n1)) -> ==(n1, @succ_nat(n2)) = false;
     (@not_equals_zero(n1)) -> ==(n1, @succ_nat(n2)) = ==(@natpred(n1),n2);
 
     <(@most_significant_digitNat(w1), @most_significant_digitNat(w2)) = @less(w1,w2);
     <(@concat_digit(n,w1), @most_significant_digitNat(w2)) = false;
     <(@most_significant_digitNat(w1),@concat_digit(n,w2)) = true;
     <(@concat_digit(n1,w1), @concat_digit(n2,w2)) = if(@less(w1,w2),<=(n1,n2),<(n1,n2));
     <(@succ_nat(n1),n2) = &&(<(@most_significant_digitNat(@one_word),n2),<(n1,@natpred(n2)));
     <(n1, @succ_nat(n2)) = <=(n1,n2);
     less_total(n,m) = <(n,m);
     @equals_zero_word(w1) ->  <(n,@most_significant_digitNat(w1)) = false;
 
     <=(@most_significant_digitNat(w1), @most_significant_digitNat(w2)) = @less_equal(w1,w2);
     <=(@concat_digit(n,w1), @most_significant_digitNat(w2)) = false;
     <=(@most_significant_digitNat(w1),@concat_digit(n,w2)) = true;
     <=(@concat_digit(n1,w1), @concat_digit(n2,w2)) = if(@less_equal(w1,w2),<=(n1,n2),<(n1,n2));
     <=(@succ_nat(n1),n2) = <(n1,n2);
     <=(n1, @succ_nat(n2)) = <=(@natpred(n1),n2);
     @equals_zero_word(w1) ->  <=(@most_significant_digitNat(w1), n) = true;
 
     Pos2Nat(@most_significant_digit(w)) = @most_significant_digitNat(w);
     Pos2Nat(@concat_digit(p,w)) = @concat_digit(Pos2Nat(p),w);
     Pos2Nat(@succ_pos(p)) = @succ_nat(Pos2Nat(p));
     Pos2Nat(succ(n)) = @succ_nat(n);
     @not_equals_zero_word(w) -> Nat2Pos(@most_significant_digitNat(w)) = @most_significant_digit(w);
     Nat2Pos(@concat_digit(n,w)) = @concat_digit(Nat2Pos(n),w);
% If important the 2 max functions below could be made more efficient by introducing a <=:Pos#Nat and <=:Nat#Pos.
     max(p, n) = if(<=(n,Pos2Nat(p)), p, Nat2Pos(n));
     max(n, p) = if(<=(n,Pos2Nat(p)), p, Nat2Pos(n));
     max(m,n) = if(<=(m,n),n,m);
     min(m,n) = if(<=(m,n),m,n);
     pred(@most_significant_digit(w)) = @most_significant_digitNat(@pred_word(w));
     pred(@concat_digit(p,w)) = if(@equals_zero_word(w),
                                   @pred_whr(pred(p)),
                                   @concat_digit(Pos2Nat(p),@pred_word(w)));
     @pred_whr(predp) = if(@equals_zero(predp),
                                    @most_significant_digitNat(@max_word),
                                    @concat_digit(predp,@max_word));

     +(p,n) = +(n,p);
     +(@most_significant_digitNat(w1),@most_significant_digit(w2)) =
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@most_significant_digit(@one_word),@add_word(w1,w2)),
                                              @most_significant_digit(@add_word(w1,w2)));
     +(@concat_digit(n1,w1),@most_significant_digit(w2)) = 
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(succ(n1),@add_word(w1,w2)),
                                              @concat_digit(Nat2Pos(n1),@add_word(w1,w2)));
     +(@most_significant_digitNat(w1),@concat_digit(p,w2)) = 
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@succ_pos(p),@add_word(w1,w2)),
                                              @concat_digit(p,@add_word(w1,w2)));
     +(@concat_digit(n1,w1),@concat_digit(p,w2)) =
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(+(succ(n1),p),@add_word(w1,w2)),
                                              @concat_digit(+(n1,p), @add_word(w1,w2)));


     +(@most_significant_digitNat(w1),@most_significant_digitNat(w2)) =
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@most_significant_digitNat(@one_word),@add_word(w1,w2)),
                                              @most_significant_digitNat(@add_word(w1,w2)));
     @add_with_carry(@most_significant_digitNat(w1),@most_significant_digitNat(w2)) =
                                     if(@add_with_carry_overflow_word(w1,w2),
                                              @concat_digit(@most_significant_digitNat(@one_word),@add_with_carry_word(w1,w2)),
                                              @most_significant_digitNat(@add_with_carry_word(w1,w2)));
     +(@concat_digit(n1,w1),@most_significant_digitNat(w2)) = 
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@succ_nat(n1),@add_word(w1,w2)),
                                              @concat_digit(n1,@add_word(w1,w2)));
     @add_with_carry(@concat_digit(n1,w1),@most_significant_digitNat(w2)) = 
                                     if(@add_with_carry_overflow_word(w1,w2),
                                              @concat_digit(@succ_nat(n1),@add_with_carry_word(w1,w2)),
                                              @concat_digit(n1,@add_with_carry_word(w1,w2)));
     +(@most_significant_digitNat(w1),@concat_digit(n2,w2)) = 
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@succ_nat(n2),@add_word(w1,w2)),
                                              @concat_digit(n2,@add_word(w1,w2)));
     @add_with_carry(@most_significant_digitNat(w1),@concat_digit(n2,w2)) = 
                                     if(@add_with_carry_overflow_word(w1,w2),
                                              @concat_digit(@succ_nat(n2),@add_with_carry_word(w1,w2)),
                                              @concat_digit(n2,@add_with_carry_word(w1,w2)));
     +(@concat_digit(n1,w1),@concat_digit(n2,w2)) =
                                     if(@add_overflow_word(w1,w2),
                                              @concat_digit(@add_with_carry(n1,n2),@add_word(w1,w2)),
                                              @concat_digit(+(n1,n2),@add_word(w1,w2)));
     @add_with_carry(@concat_digit(n1,w1),@concat_digit(n2,w2)) =
                                     if(@add_with_carry_overflow_word(w1,w2),
                                              @concat_digit(@add_with_carry(n1,n2),@add_with_carry_word(w1,w2)),
                                              @concat_digit(+(n1,n2),@add_with_carry_word(w1,w2)));

% The rules below are useful in solving expressions with plus and quantifiers.
     +(@succ_nat(n1),n2) = @succ_nat(+(n1,n2));
     +(n1,@succ_nat(n2)) = @succ_nat(+(n1,n2));
     +(@succ_nat(n1),p2) = @succ_pos(+(n1,p2));
     +(n1,@succ_pos(p2)) = @succ_pos(+(n1,p2));
     +(@succ_pos(p1),n2) = @succ_pos(+(p1,n2));
     +(p1,@succ_nat(n2)) = @succ_pos(+(p1,n2));
     +(@most_significant_digitNat(@zero_word),n) = n;
     +(n, @most_significant_digitNat(@zero_word)) = n;
     +(@most_significant_digitNat(@zero_word),p) = p;
     +(p, @most_significant_digitNat(@zero_word)) = p;

     @plus_nat(n1,n2) = +(n1,n2);
 
     @natpred(@most_significant_digitNat(w)) = if(@equals_zero_word(w), @most_significant_digitNat(@zero_word), @most_significant_digitNat(@pred_word(w)));
     @natpred(@concat_digit(n,w)) = if(@equals_zero_word(w),
                                       if(@equals_one(n),
                                               @most_significant_digitNat(@max_word),
                                               @concat_digit(@natpred(n),@max_word)),
                                       @concat_digit(n,@pred_word(w)));
     @natpred(@succ_nat(n)) = n;



     @monus(@most_significant_digitNat(w1),@most_significant_digitNat(w2)) = @most_significant_digitNat(@monus_word(w1,w2));

     @monus(@concat_digit(n1,w1),@most_significant_digitNat(w2)) = if(@less(w1,w2),
                                                                      if(@equals_one(n1),
                                                                               @most_significant_digitNat(@minus_word(w1,w2)),
                                                                               @concat_digit(@natpred(n1),@minus_word(w1,w2))),
                                                                      @concat_digit(n1,@minus_word(w1,w2)));
     @monus(@most_significant_digitNat(w1),@concat_digit(n2,w2)) = @most_significant_digitNat(@zero_word);
     @monus(@concat_digit(n1,w1),@concat_digit(n2,w2)) = @monus_whr(n1,w1,n2,w2,@monus(n1,n2));
     @monus_whr(n1,w1,n2,w2,diff) = if(@less(w1,w2),
                                        if(@equals_zero(diff),
                                               @most_significant_digitNat(@zero_word),
                                               if(@equals_one(diff),
                                                        @most_significant_digitNat(@minus_word(w1,w2)),
                                                        @concat_digit(@natpred(diff),@minus_word(w1,w2)))),
                                        if(@equals_zero(diff),
                                               @most_significant_digitNat(@minus_word(w1,w2)),
                                               @concat_digit(diff,@minus_word(w1,w2))));

     *(@most_significant_digitNat(w1),@most_significant_digitNat(w2)) =
                           if(@equals_zero_word(@times_overflow_word(w1,w2)),
                                     @most_significant_digitNat(@times_word(w1,w2)),
                                     @concat_digit(@most_significant_digitNat(@times_overflow_word(w1,w2)),@times_word(w1,w2)));
     *(@most_significant_digitNat(w1),@concat_digit(n2,w2)) =
                           if(@equals_zero_word(w1),
                                @most_significant_digitNat(@zero_word),
                                @concat_digit(@times_overflow(n2,w1,@times_overflow_word(w1,w2)),@times_word(w1,w2)));
     *(@concat_digit(n1,w1),@most_significant_digitNat(w2)) = 
                                if(@equals_zero_word(w2),
                                         @most_significant_digitNat(@zero_word),
                                         @concat_digit(@times_overflow(n1,w2,@times_overflow_word(w1,w2)),@times_word(w1,w2)));
     *(@concat_digit(n1,w1),@concat_digit(n2,w2)) = 
                                if(<(n1,n2),
                                         @times_ordered(@concat_digit(n1,w1),@concat_digit(n2,w2)),
                                         @times_ordered(@concat_digit(n2,w2),@concat_digit(n1,w1)));

% In @times_ordered, the lhs is not equal to zero and the rhs has more digits than the rhs, always at least two.
%     @times_ordered(@most_significant_digitNat(w1),@most_significant_digitNat(w2)) =
%                           if(@equals_zero_word(@times_overflow_word(w1,w2)),
%                                     @most_significant_digitNat(@times_word(w1,w2)),
%                                     @concat_digit(@most_significant_digitNat(@times_overflow_word(w1,w2)),@times_word(w1,w2)));
     @times_ordered(@most_significant_digitNat(w1),@concat_digit(n2,w2)) =
                                @concat_digit(@times_overflow(n2,w1,@times_overflow_word(w1,w2)),@times_word(w1,w2));
     @times_ordered(@concat_digit(n1,w1),n2) = +(@concat_digit(@times_ordered(n1,n2),@zero_word), @times_overflow(n2,w1,@zero_word));

     @times_overflow(@most_significant_digitNat(w1),w2,overflow) =
                           if(@equals_zero_word(@times_with_carry_overflow_word(w1,w2,overflow)),
                                     @most_significant_digitNat(@times_with_carry_word(w1,w2,overflow)),
                                     @concat_digit(@most_significant_digitNat(@times_with_carry_overflow_word(w1,w2,overflow)),
                                                   @times_with_carry_word(w1,w2,overflow)));
     @times_overflow(@concat_digit(n1,w1),w2,overflow) =
                           if(@equals_zero_word(w2),
                                @most_significant_digitNat(overflow),
                                @concat_digit(@times_overflow(n1,w2,@times_with_carry_overflow_word(w1,w2,overflow)),
                                              @times_with_carry_word(w1,w2,overflow)));

     @is_odd(@most_significant_digitNat(w)) = @rightmost_bit(w);
     @is_odd(@concat_digit(n,w)) = @rightmost_bit(w);

     @div2(@most_significant_digitNat(w)) = @most_significant_digitNat(@shift_right(false,w));
     @div2(@concat_digit(n,w)) = if(@equals_zero(n),
                                     @most_significant_digitNat(@shift_right(@is_odd(n),w)),
                                     @concat_digit(@div2(n),@shift_right(@is_odd(n),w)));

     @msd(@most_significant_digitNat(w)) = w;
     @msd(@concat_digit(n,w)) = @msd(n);

     exp(n,@most_significant_digitNat(w)) = @exp_aux3n(@rightmost_bit(w),n,w);
     exp(n,@concat_digit(n1,w1)) = @exp_aux4n(@rightmost_bit(w1),n,n1,w1);
 
     @exp_aux3n(true,n,w) = if(@equals_one_word(w),
                                n,
                                *(n,@exp_aux3n(@rightmost_bit(@shift_right(false,w)),*(n,n),@shift_right(false,w))));
 
     @exp_aux3n(false,n,w) = if(@equals_zero_word(w),
                                @most_significant_digitNat(@one_word),
                                @exp_aux3n(@rightmost_bit(@shift_right(false,w)),*(n,n),@shift_right(false,w)));
 
     @exp_aux4n(true,n,n1,w) = @exp_auxtruen(n,@div2(n1),@shift_right(@is_odd(n1), w));
     @exp_auxtruen(n,shift_n1,shift_w) =
                            if(@equals_zero(shift_n1),
                                *(n,@exp_aux3n(@rightmost_bit(shift_w),*(n,n),shift_w)),
                                *(n,@exp_aux4n(@rightmost_bit(shift_w),*(n,n),shift_n1,shift_w)));

     @exp_aux4n(false,n,n1,w) = @exp_auxfalsen(n,@div2(n1),@shift_right(@is_odd(n1),w));
     @exp_auxfalsen(n,shift_n1,shift_w) =
                            if(@equals_zero(shift_n1),
                                @exp_aux3n(@rightmost_bit(shift_w),*(n,n),shift_w),
                                @exp_aux4n(@rightmost_bit(shift_w),*(n,n),shift_n1,shift_w));
 
     exp(p,@most_significant_digitNat(w)) = @exp_aux3p(@rightmost_bit(w),p,w);
     exp(p,@concat_digit(n1,w1)) = @exp_aux4p(@rightmost_bit(w1),p,n1,w1);
 
     @exp_aux3p(true,p,w) = if(@equals_one_word(w),
                                p,
                                *(p,@exp_aux3p(@rightmost_bit(@shift_right(false,w)),*(p,p),@shift_right(false,w))));
                                
     @exp_aux3p(false,p,w) = if(@equals_zero_word(w),
                                 @most_significant_digit(@one_word),
                                 @exp_aux3p(@rightmost_bit(@shift_right(false,w)),*(p,p),@shift_right(false,w)));

     @exp_aux4p(true,p,n1,w) = @exp_auxtruep(p,@div2(n1),@shift_right(@is_odd(n1), w));
     @exp_auxtruep(p,shift_n1,shift_w) =
                           if(@equals_zero(shift_n1),
                                *(p,@exp_aux3p(@rightmost_bit(shift_w),*(p,p),shift_w)),
                                *(p,@exp_aux4p(@rightmost_bit(shift_w),*(p,p),shift_n1,shift_w)));

     @exp_aux4p(false,p,n1,w) = @exp_auxfalsep(p,@div2(n1),@shift_right(@is_odd(n1),w));
     @exp_auxfalsep(p,shift_n1,shift_w) =
                           if(@equals_zero(shift_n1),
                                @exp_aux3p(@rightmost_bit(shift_w),*(p,p),shift_w),
                                @exp_aux4p(@rightmost_bit(shift_w),*(p,p),shift_n1,shift_w));
 

    div(@most_significant_digitNat(w1),@most_significant_digit(w2)) = @most_significant_digitNat(@div_word(w1,w2));
    mod(@most_significant_digitNat(w1),@most_significant_digit(w2)) = @most_significant_digitNat(@mod_word(w1,w2));

    div(@most_significant_digitNat(w1),@concat_digit(p,w2)) = @most_significant_digitNat(@zero_word);
    mod(@most_significant_digitNat(w1),@concat_digit(p,w2)) = @most_significant_digitNat(w1);

%    >(n,@most_significant_digitNat(@zero_word)) ->    Only checks the proper shape of a natural number.
    div(@concat_digit(n,w1),@most_significant_digit(w2)) = @div_whr1(n,w1,w2,@divmod_aux(n,@most_significant_digit(w2)));

    @div_whr1(n,w1,w2,pair_) = 
                          if(<(n,@most_significant_digitNat(w2)),
                                    @most_significant_digitNat(@div_bold(@concat_digit(n,w1),@most_significant_digit(w2))),
                                    if (@equals_zero(@first(pair_)),
                                        @most_significant_digitNat(
                                           @div_bold(
                                               if(@equals_zero(@last(pair_)), 
                                                   @most_significant_digitNat(w1),
                                                   @concat_digit(@last(pair_),w1)),
                                               @most_significant_digit(w2))),
                                        @concat_digit(@first(pair_),
                                           @div_bold(
                                               if(@equals_zero(@last(pair_)), 
                                                   @most_significant_digitNat(w1),
                                                   @concat_digit(@last(pair_),w1)),
                                               @most_significant_digit(w2)))));

%    >(n,@most_significant_digitNat(@zero_word)) ->     Only checks the proper shape of a natural number.
    mod(@concat_digit(n,w1),@most_significant_digit(w2)) = @most_significant_digitNat(@mod_doubleword(@msd(mod(n,@most_significant_digit(w2))),w1,w2));

%    >(n,@most_significant_digitNat(@zero_word)) ->     Only checks the proper shape of a natural number
    div(@concat_digit(n,w1),@concat_digit(p,w2)) =
                          if(<(n,Pos2Nat(@concat_digit(p,w2))),
                                  @most_significant_digitNat(@div_bold(@concat_digit(n,w1),@concat_digit(p,w2))),
                                  @div_whr2(n,w1,p,w2,@divmod_aux(n,@concat_digit(p,w2))));
    @div_whr2(n,w1,p,w2,pair_) = 
                          +(if(@equals_zero(@first(pair_)), 
                                    @most_significant_digitNat(@zero_word),
                                    @concat_digit(@first(pair_),@zero_word)),
                               @most_significant_digitNat(@div_bold(if(@equals_zero(@last(pair_)), 
                                                                         @most_significant_digitNat(w1),
                                                                         @concat_digit(@last(pair_),w1)),
                                                          @concat_digit(p,w2))));

%    >(n,@most_significant_digitNat(@zero_word)) ->    Only checks the proper shape of a natural number. 
    mod(@concat_digit(n,w1),@concat_digit(p,w2)) = @mod_whr1(w1,p,w2,mod(n,@concat_digit(p,w2)));

    @mod_whr1(w1,p,w2,m1) = 
                          @monus(if(<(@most_significant_digitNat(@zero_word),m1),@concat_digit(m1,w1),@most_significant_digitNat(w1)),
                                 *(@concat_digit(Pos2Nat(p),w2),
                                   @most_significant_digitNat(
                                        @div_bold(
                                             if(<(@most_significant_digitNat(@zero_word),m1),@concat_digit(m1,w1),@most_significant_digitNat(w1)),
                                             @concat_digit(p,w2)))));

    @divmod_aux(@most_significant_digitNat(w1),@most_significant_digit(w2)) = 
                                      @nnPair(@most_significant_digitNat(@div_word(w1,w2)),
                                              @most_significant_digitNat(@mod_word(w1,w2)));

    @divmod_aux(@most_significant_digitNat(w1),@concat_digit(p,w2)) = @nnPair(@most_significant_digitNat(@zero_word),@most_significant_digitNat(w1));

%    >(n,@most_significant_digitNat(@zero_word)) ->    Only checks the proper shape of a natural number.
    @divmod_aux(@concat_digit(n,w1),@most_significant_digit(w2)) = @divmod_aux_whr1(n,w1,w2,@divmod_aux(n,@most_significant_digit(w2)));
    
    @divmod_aux_whr1(n,w1,w2,pair_) =
                          @nnPair(if(<(n,@most_significant_digitNat(w2)),
                                    @most_significant_digitNat(@div_bold(@concat_digit(n,w1),@most_significant_digit(w2))),
                                    if (@equals_zero(@first(pair_)),
                                        @most_significant_digitNat(@div_bold(
                                            if(@equals_zero(@last(pair_)), 
                                                 @most_significant_digitNat(w1),
                                                 @concat_digit(@last(pair_),w1)),@most_significant_digit(w2))),
                                        @concat_digit(@first(pair_),@div_bold(
                                            if(@equals_zero(@last(pair_)), 
                                                 @most_significant_digitNat(w1),
                                                 @concat_digit(@last(pair_),w1)),@most_significant_digit(w2))))),
                                  @most_significant_digitNat(@mod_doubleword(@msd(@last(pair_)),w1,w2)));

%    >(n,@most_significant_digitNat(@zero_word)) ->    Only checks whether natural numbers are of the proper shape, which must always be the case.
    @divmod_aux(@concat_digit(n,w1),@concat_digit(p,w2)) =
                          if(<(n,Pos2Nat(@concat_digit(p,w2))),
                             @divmod_aux_whr2(n,w1,p,w2,mod(n,@concat_digit(p,w2))),
                             @divmod_aux_whr4(w1,p,w2,@divmod_aux(n,@concat_digit(p,w2))));

    @divmod_aux_whr2(n,w1,p,w2,lp) = @divmod_aux_whr3(n,w1,p,w2,if(@equals_zero(lp),@most_significant_digitNat(w1),@concat_digit(lp,w1)));
    @divmod_aux_whr3(n,w1,p,w2,m) = @nnPair(@most_significant_digitNat(@div_bold(@concat_digit(n,w1),@concat_digit(p,w2))),
                                            @monus(m,*(@concat_digit(Pos2Nat(p),w2),@most_significant_digitNat(@div_bold(m,@concat_digit(p,w2))))));
    @divmod_aux_whr4(w1,p,w2,pair_) = 
                         @divmod_aux_whr5(p,w2,pair_,
                                          if(@equals_zero(@last(pair_)),@most_significant_digitNat(w1),@concat_digit(@last(pair_),w1)));
    @divmod_aux_whr5(p,w2,pair_,m) = @divmod_aux_whr6(p,w2,pair_,m,@most_significant_digitNat(@div_bold(m,@concat_digit(p,w2))));

    @divmod_aux_whr6(p,w2,pair_,m,m1) = 
                         @nnPair(+(@concat_digit(@first(pair_),@zero_word),m1),
                                 @monus(m,*(@concat_digit(Pos2Nat(p),w2),m1)));

    @div_bold(@most_significant_digitNat(w1),@most_significant_digit(w2)) = @div_word(w1,w2);

    @div_bold(@most_significant_digitNat(w1),@concat_digit(p,w2)) = @zero_word;

    @div_bold(@concat_digit(@most_significant_digitNat(w1),w2),@most_significant_digit(w3)) = @div_doubleword(w1,w2,w3);

    @div_bold(@concat_digit(@most_significant_digitNat(w),w1),@concat_digit(@most_significant_digit(w2),w3)) =
                                               @div_double_doubleword(w,w1,w2,w3);

    @div_bold(@concat_digit(@concat_digit(@most_significant_digitNat(w),w1),w2),@concat_digit(@most_significant_digit(w3),w4)) =
                                               @div_triple_doubleword(w,w1,w2,w3,w4);

%    >(n,@most_significant_digitNat(@zero_word)) ->   Only checks if term has proper shape, which must always be the case. 
    @div_bold(@concat_digit(n,w1),@concat_digit(@concat_digit(p,w2),w3)) =
                          if(<=(@concat_digit(Pos2Nat(@concat_digit(p,w2)),@zero_word),n),
                             @max_word,
                             @div_bold_whr(n,w1,p,w2,w3,@div_bold(n,@concat_digit(p,w2))));

    @div_bold_whr(n,w1,p,w2,w3,w) = 
                             if(<(@concat_digit(n,w1),*(@most_significant_digitNat(w),Pos2Nat(@concat_digit(@concat_digit(p,w2),w3)))),
                                @pred_word(w),
                                w);


% equations for pairs
     ==(@nnPair(n1,n2), @nnPair(m1,m2)) = &&(==(n1,m1),==(n2,m2));
     <(@nnPair(n1,n2), @nnPair(m1,m2)) = ||(<(n1,m1), &&(==(n1,m1), <(n2,m2)));
     <=(@nnPair(n1,n2), @nnPair(m1,m2)) = ||(<(n1,m1), &&(==(n1,m1), <=(n2,m2)));
     @first(@nnPair(m,n)) = m;
     @last(@nnPair(m,n)) = n;

%Residues. 

%     @swap_zero(m,@c0) = m;
%     @swap_zero(@c0,n) = n;
%     @swap_zero(@cNat(p),@cNat(p)) = @c0;
%     !=(p,q) -> @swap_zero(@cNat(p),@cNat(q)) = @cNat(q);
     @swap_zero(m,n) = if(@equals_zero(n), m,
                       if(@equals_zero(m), n,
                       if(==(n,m), @most_significant_digitNat(@zero_word), n)));

%     @swap_zero_add(@c0, @c0, m, n) = +(m,n);
%     @swap_zero_add(@c0, @cNat(p), @c0, n) = n;
%     @swap_zero_add(@c0, @cNat(p), @cNat(q), n) = @swap_zero(@cNat(p), +(@cNat(q), @swap_zero(@cNat(p), n)));
%     @swap_zero_add(@cNat(p), @c0, m, @c0) = m;
%     @swap_zero_add(@cNat(p), @c0, m, @cNat(q)) = @swap_zero(@cNat(p), +(@swap_zero(@cNat(p),m), @cNat(q)));
%     @swap_zero_add(@cNat(p), @cNat(q), m, n) = @swap_zero(+(@cNat(p), @cNat(q)), +(@swap_zero(@cNat(p),m),@swap_zero(@cNat(q),n)));
     @swap_zero_add(n1,n2,m1,m2) = 
                       if(@equals_zero(n1),
                            if(@equals_zero(n2), 
                                 +(m1,m2),
                                 if(@equals_zero(m1),
                                    m2,
                                    @swap_zero(n2, +(m1, @swap_zero(n2, m2))))),
                            if(@equals_zero(n2), 
                                 if(@equals_zero(m2),
                                      n1,
                                      @swap_zero(n1, +(@swap_zero(n1,m1), m2))),
                                 @swap_zero(+(n1, n2), +(@swap_zero(n1,m1),@swap_zero(n2,m2)))));

%     @swap_zero_min(@c0, @c0, m, n) = min(m,n);
%     @swap_zero_min(@c0, @cNat(p), @c0, n) = @c0;
%     @swap_zero_min(@c0, @cNat(p), @cNat(q), n) = min(@cNat(q), @swap_zero(@cNat(p),n));
%     @swap_zero_min(@cNat(p), @c0, m, @c0) = @c0;
%     @swap_zero_min(@cNat(p), @c0, m, @cNat(q)) = min(@swap_zero(@cNat(p),m), @cNat(q));
%     @swap_zero_min(@cNat(p), @cNat(q), m, n) = @swap_zero(min(@cNat(p), @cNat(q)), min(@swap_zero(@cNat(p), m), @swap_zero(@cNat(q), n)));
     @swap_zero_min(n1, n2, m1, m2) =
                       if(@equals_zero(n1),
                            if(@equals_zero(n2), 
                                 min(m1,m2),
                                 if(@equals_zero(m1),
                                    @most_significant_digitNat(@zero_word),
                                    min(m1, @swap_zero(n2,m2)))),
                            if(@equals_zero(n2), 
                                 if(@equals_zero(m2),
                                      @most_significant_digitNat(@zero_word),
                                      min(@swap_zero(n1,m1), m2)),
                                 @swap_zero(min(n1, n2), min(@swap_zero(n1, m1), @swap_zero(n2, m2)))));

%     @swap_zero_monus(@c0, @c0, m, n) = @monus(m,n);
%     @swap_zero_monus(@c0, @cNat(p), @c0, n) = @c0;
%     @swap_zero_monus(@c0, @cNat(p), @cNat(q), n) = @monus(@cNat(q), @swap_zero(@cNat(p), n));
%     @swap_zero_monus(@cNat(p), @c0, m, @c0) = m;
%     @swap_zero_monus(@cNat(p), @c0, m, @cNat(q)) = @swap_zero(@cNat(p), @monus(@swap_zero(@cNat(p), m), @cNat(q)));
%     @swap_zero_monus(@cNat(p), @cNat(q), m, n) = @swap_zero(@monus(@cNat(p),@cNat(q)),@monus(@swap_zero(@cNat(p),m), @swap_zero(@cNat(q),n)));
     @swap_zero_monus(n1, n2, m1, m2) =
                       if(@equals_zero(n1),
                            if(@equals_zero(n2), 
                                 @monus(m1,m2),
                                 if(@equals_zero(m1),
                                    @most_significant_digitNat(@zero_word),
                                    @monus(m1, @swap_zero(n2, m2)))),
                            if(@equals_zero(n2), 
                                 if(@equals_zero(m2),
                                      m1,
                                      @swap_zero(n1, @monus(@swap_zero(n1, m1), m2))),
                                 @swap_zero(@monus(n1,n2),@monus(@swap_zero(n1,m1), @swap_zero(n2,m2)))));

     sqrt(@most_significant_digitNat(w)) = @most_significant_digitNat(@sqrt_word(w));
     sqrt(@concat_digit(@most_significant_digitNat(w1),w2)) = @most_significant_digitNat(@sqrt_doubleword(w1,w2));
     sqrt(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3)) = @sqrt_whr1(w1,w2,w3,@sqrt_tripleword_overflow(w1,w2,w3));
     @sqrt_whr1(w1,w2,w3,overflow) =
              if (@equals_zero_word(overflow),
                  @most_significant_digitNat(@sqrt_tripleword(w1,w2,w3)),
                  @concat_digit(@most_significant_digitNat(overflow),@sqrt_tripleword(w1,w2,w3)));

     sqrt(@concat_digit(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3),w4)) =
                         @sqrt_whr2(w1,w2,w3,w4,@sqrt_quadrupleword_overflow(w1,w2,w3,w4));
     @sqrt_whr2(w1,w2,w3,w4,overflow) =
              if (@equals_zero_word(overflow),
                  @most_significant_digitNat(@sqrt_quadrupleword(w1,w2,w3,w4)),
                  @concat_digit(@most_significant_digitNat(overflow),@sqrt_quadrupleword(w1,w2,w3,w4)));

     sqrt(@concat_digit(@concat_digit(@concat_digit(@concat_digit(n,w1),w2),w3),w4)) = @first(@sqrt_pair(@concat_digit(@concat_digit(@concat_digit(@concat_digit(n,w1),w2),w3),w4)));

     @sqrt_pair(@most_significant_digitNat(w)) = @nnPair(@most_significant_digitNat(@sqrt_word(w)), @most_significant_digitNat(@minus_word(w,@times_word(@sqrt_word(w),@sqrt_word(w)))));
     @sqrt_pair(@concat_digit(@most_significant_digitNat(w1),w2)) =
                 @nnPair(@most_significant_digitNat(@sqrt_doubleword(w1,w2)),
                         @monus(@concat_digit(@most_significant_digitNat(w1),w2),exp(@most_significant_digitNat(@sqrt_doubleword(w1,w2)),@most_significant_digitNat(@two_word))));

     @sqrt_pair(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3)) =
                 @sqrt_pair_whr1(w1,w2,w3,sqrt(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3)));
     @sqrt_pair_whr1(w1,w2,w3,solution) =      
                    @nnPair(solution,
                            @monus(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3),*(solution,solution)));

     @sqrt_pair(@concat_digit(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3),w4)) =
                  @sqrt_pair_whr2(w1,w2,w3,w4,sqrt(@concat_digit(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3),w4))); 
     @sqrt_pair_whr2(w1,w2,w3,w4,solution) =
                    @nnPair(solution,
                            @monus(@concat_digit(@concat_digit(@concat_digit(@most_significant_digitNat(w1),w2),w3),w4), *(solution,solution)));

     @sqrt_pair(@concat_digit(@concat_digit(@concat_digit(@concat_digit(n,w1),w2),w3),w4)) =
                    @sqrt_pair_whr3(w3,w4,@sqrt_pair(@concat_digit(@concat_digit(n,w1),w2)));

% pq = @nnPair(p,q) = sqrt(n).
     @sqrt_pair_whr3(w3,w4,pq) =
                 @sqrt_pair_whr4(w3,w4,pq,
                                 if(@equals_zero(@first(pq)), @most_significant_digitNat(@zero_word), @concat_digit(@first(pq),@zero_word)),
                                 if(<(@most_significant_digitNat(@zero_word),@last(pq)),
                                      @concat_digit(@concat_digit(@last(pq),w3),w4),
                                      if(@not_equals_zero_word(w3),
                                            @concat_digit(@most_significant_digitNat(w3),w4),
                                            @most_significant_digitNat(w4))),
                                 +(*(@first(pq),@first(pq)),@last(pq)));


% m3 = p*base.
% m2 = q*base^2 + w3*base + w4.
% m5 = p^2 + q.
     @sqrt_pair_whr4(w3,w4,pq,m3,m2,m5) = 
                @sqrt_pair_whr5(pq,m3,m2,
                            if(<(@most_significant_digitNat(@zero_word),m5),
                                 @concat_digit(@concat_digit(m5,w3),w4),
                                 if(@not_equals_zero_word(w3),
                                      @concat_digit(@most_significant_digitNat(w3),w4),
                                      @most_significant_digitNat(w4))),
                             div(m2,Nat2Pos(*(@most_significant_digitNat(@two_word),m3))));
                               

% m4 = (p^2 + q)*base^2 + w3*base + w4.
% y_guess = (q*base^2 + w3*base + w4) div (2*p*base).
     @sqrt_pair_whr5(pq,m3,m2,m4,y_guess) = 
                @sqrt_pair_whr6(m2,
                                if(<(*(@most_significant_digitNat(@four_word),@first(pq)),
                                     @concat_digit(@most_significant_digitNat(@three_word), @zero_word)),
                                   @monus(sqrt(m4),m3),
                                   if(<(m2,*(+(*(@most_significant_digitNat(@two_word),m3),y_guess),y_guess)),
                                       @natpred(y_guess),
                                       y_guess)),
                                if(@equals_zero(@first(pq)),@most_significant_digitNat(@zero_word),@concat_digit(@first(pq),@zero_word)));


% m1 = p*base.
     @sqrt_pair_whr6(m2,y,m1) = 
                    @nnPair(+(m1,y), @monus(m2, *(+(*(m1,@most_significant_digitNat(@two_word)),y),y)));
