% Author(s): Jan Friso Groote
% Copyright: see the accompanying file COPYING or copy at
% https://github.com/mCRL2org/mCRL2/blob/master/COPYING
%
% Distributed under the Boost Software License, Version 1.0.
% (See accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)
%
% Specification of the Pos data sort that uses machine numbers for efficiency.

#include bool.spec
#include machine_word.spec

sort Pos <"pos">;

cons @c1 <"c1">: Pos                                                                                  internal defined_by_rewrite_rules;
% The successor constructor should be merged with the successor below, by removing the latter.
% Currently, this does not work, as the translator to C code does not see that there are multiple
% successor functions of different types, as this one is a constructor. 
     @succ_pos <"succpos">:Pos <"arg"> -> Pos                                                         internal defined_by_rewrite_rules;

map  @most_significant_digit <"most_significant_digit">: @word <"arg"> -> Pos                         internal defined_by_rewrite_rules;
     @concat_digit <"concat_digit">: Pos <"arg1"> # @word <"arg2"> -> Pos                             internal defined_by_rewrite_rules;
     @equals_one <"equals_one">: Pos <"arg"> -> Bool                                                  internal defined_by_rewrite_rules;
     max <"maximum">: Pos <"left"> # Pos <"right"> -> Pos                                             external defined_by_rewrite_rules;
     min <"minimum">: Pos <"left"> # Pos <"right"> -> Pos                                             external defined_by_rewrite_rules;
% There is a special mapping succ, as overloading a constructor is not possible. Therefore the constructor @succ_pos has a unique name. 
     succ <"succ">: Pos <"arg"> -> Pos                                                                external defined_by_rewrite_rules;
     @pospred <"pos_predecessor">: Pos <"arg"> -> Pos                                                 internal defined_by_rewrite_rules;
     + <"plus">: Pos <"left"> # Pos <"right"> -> Pos                                                  external defined_by_rewrite_rules;
     @add_with_carry <"add_with_carry">: Pos <"left"> # Pos <"right"> -> Pos                          internal defined_by_rewrite_rules;
% The following function is used when the symbol + is overloaded, such as in fbags. 
     @plus_pos <"auxiliary_plus_pos">: Pos <"left"> # Pos <"right"> -> Pos                            external defined_by_rewrite_rules;
     * <"times">: Pos <"left"> # Pos <"right"> -> Pos                                                 external defined_by_rewrite_rules;
     @times_overflow <"times_overflow">: Pos <"arg1"> # @word <"arg2"> # @word <"arg3"> -> Pos        external defined_by_rewrite_rules;
     @times_ordered <"times_ordered">: Pos <"left"> # Pos <"right"> -> Pos                            external defined_by_rewrite_rules;
% Auxiliary function to implement multiplication that uses where clauses.
     @times_whr_mult_overflow <"times_whr_mult_overflow">: @word <"arg1"> # @word <"arg2"> -> Pos     internal defined_by_rewrite_rules;

var  b: Bool;
     p: Pos;
     p1: Pos;
     p2: Pos;
     w: @word;
     w1: @word;
     w2: @word;
     overflow: @word;
eqn  @c1 = @most_significant_digit(@one_word);
     @equals_one(@most_significant_digit(w)) = @equals_one_word(w);
     @equals_one(@concat_digit(p,w)) = false;
     @equals_one(@succ_pos(p)) = false;
     succ(p) = @succ_pos(p);
     @succ_pos(@most_significant_digit(w1)) = if(@equals_max_word(w1),
                                                 @concat_digit(@most_significant_digit(@one_word),@zero_word),
                                                 @most_significant_digit(@succ_word(w1)));
     @succ_pos(@concat_digit(p,w1)) = if(@equals_max_word(w1),
                                             @concat_digit(@succ_pos(p),@zero_word),
                                             @concat_digit(p,@succ_word(w1)));
 
% The rules for comparison operators must be explicitly be defined on succ (= @succ_pos) to allow enumeration of positive numbers.
     ==(@most_significant_digit(w1),@most_significant_digit(w2)) = @equal(w1,w2);
     ==(@concat_digit(p,w1),@most_significant_digit(w2)) = false;
     ==(@most_significant_digit(w1),@concat_digit(p,w2)) = false;
     ==(@concat_digit(p1,w1),@concat_digit(p2,w2)) = &&(@equal(w1,w2),==(p1,p2));
     (@equals_one(p2)) -> ==(@succ_pos(p1),p2) = false;
     !(@equals_one(p2)) -> ==(@succ_pos(p1),p2) = ==(p1,@pospred(p2));
     (@equals_one(p1)) -> ==(p1, @succ_pos(p2)) = false;
     !(@equals_one(p1)) -> ==(p1, @succ_pos(p2)) = ==(@pospred(p1),p2);
 
     <(@most_significant_digit(w1),@most_significant_digit(w2)) = @less(w1,w2);
     <(@concat_digit(p,w1),@most_significant_digit(w2)) = false;
     <(@most_significant_digit(w1),@concat_digit(p,w2)) = true;
     <(@concat_digit(p1,w1),@concat_digit(p2,w2)) = if(@less(w1,w2),<=(p1,p2),<(p1,p2));
     <(@succ_pos(p1),p2) = &&(<(@most_significant_digit(@two_word),p2),<(p1,@pospred(p2)));
     <(p1, @succ_pos(p2)) = <=(p1,p2);
     less_total(p,p1) = <(p,p1);
     @equals_one_word(w1) ->  <(p,@most_significant_digit(w1)) = false;
 
     <=(@most_significant_digit(w1),@most_significant_digit(w2)) = @less_equal(w1,w2);
     <=(@concat_digit(p,w1),@most_significant_digit(w2)) = false;
     <=(@most_significant_digit(w1),@concat_digit(p,w2)) = true;
     <=(@concat_digit(p1,w1),@concat_digit(p2,w2)) = if(@less_equal(w1,w2),<=(p1,p2),<(p1,p2));
     <=(@succ_pos(p1),p2) = <(p1,p2);
     <=(p1, @succ_pos(p2)) = ||(@equals_one(p1),<=(@pospred(p1),p2));
     @equals_one_word(w1) ->  <=(@most_significant_digit(w1), p) = true;
 
     max(p1,p2) = if(<=(p1,p2),p2,p1);
     min(p1,p2) = if(<=(p1,p2),p1,p2);

     @pospred(@most_significant_digit(w1)) = if(@equals_one_word(w1),
                                                         @most_significant_digit(@one_word),
                                                         @most_significant_digit(@pred_word(w1)));
     @pospred(@concat_digit(p,w1)) = if(@equals_zero_word(w1),
                                            if(@equals_one(p),
                                                         @most_significant_digit(@max_word),
                                                         @concat_digit(@pospred(p),@max_word)),
                                            @concat_digit(p,@pred_word(w1)));
     @pospred(@succ_pos(p)) = p;

     +(@most_significant_digit(w1),@most_significant_digit(w2)) =
                                     if(@add_overflow_word(w1,w2),
                                                    @concat_digit(@most_significant_digit(@one_word),@add_word(w1,w2)),
                                                    @most_significant_digit(@add_word(w1,w2)));
     @add_with_carry(@most_significant_digit(w1),@most_significant_digit(w2)) =
                                     if(@add_with_carry_overflow_word(w1,w2),
                                                    @concat_digit(@most_significant_digit(@one_word),(@add_with_carry_word(w1,w2))),
                                                    @most_significant_digit(@add_with_carry_word(w1,w2)));

% The rules below are not efficient as a zero carry ripples through the whole term. 
%     +(@concat_digit(p1,w1),@most_significant_digit(w2)) = @concat_digit(+(@most_significant_digit(@add_overflow_word(w1,w2)),p1),
%                                                                         @add_word(w1,w2));
%     +(@most_significant_digit(w1),@concat_digit(p2,w2)) = @concat_digit(+(@most_significant_digit(@add_overflow_word(w1,w2)), p2),
%                                                                         @add_word(w1,w2));
%     +(@concat_digit(p1,w1),@concat_digit(p2,w2)) = @concat_digit(+(@most_significant_digit(@add_overflow_word(w1,w2)), +(p1,p2)),
%                                                                       @add_word(w1,w2));
     +(@concat_digit(p1,w1),@most_significant_digit(w2)) = 
                   if(@add_overflow_word(w1,w2),
                                       @concat_digit(@succ_pos(p1),@add_word(w1,w2)),
                                       @concat_digit(p1, @add_word(w1,w2)));
     @add_with_carry(@concat_digit(p1,w1),@most_significant_digit(w2)) = 
                   if(@add_with_carry_overflow_word(w1,w2),
                                       @concat_digit(@succ_pos(p1),@add_with_carry_word(w1,w2)),
                                       @concat_digit(p1, @add_with_carry_word(w1,w2)));

     +(@most_significant_digit(w1),@concat_digit(p2,w2)) = 
                   if(@add_overflow_word(w1,w2),
                                       @concat_digit(@succ_pos(p2),@add_word(w1,w2)),
                                       @concat_digit(p2, @add_word(w1,w2)));
     @add_with_carry(@most_significant_digit(w1),@concat_digit(p2,w2)) = 
                   if(@add_with_carry_overflow_word(w1,w2),
                                       @concat_digit(@succ_pos(p2),@add_with_carry_word(w1,w2)),
                                       @concat_digit(p2, @add_with_carry_word(w1,w2)));

     +(@concat_digit(p1,w1),@concat_digit(p2,w2)) = 
                   if(@add_overflow_word(w1,w2),
                                       @concat_digit(@add_with_carry(p1,p2), @add_word(w1,w2)),
                                       @concat_digit(+(p1,p2), @add_word(w1,w2)));
     @add_with_carry(@concat_digit(p1,w1),@concat_digit(p2,w2)) = 
                   if(@add_with_carry_overflow_word(w1,w2),
                                       @concat_digit(@add_with_carry(p1,p2), @add_with_carry_word(w1,w2)),
                                       @concat_digit(+(p1,p2), @add_with_carry_word(w1,w2)));

% The rules below are useful in solving expressions with plus and quantifiers.
     +(@succ_pos(p1),p2) = @succ_pos(+(p1,p2));
     +(p1,@succ_pos(p2)) = @succ_pos(+(p1,p2));

     @plus_pos(p1,p2) = +(p1,p2);

% The definition below uses where clauses. The where clauses are translated away by the introduction of the 
% function @times_whr_mult_overflow. 
%      *(@most_significant_digit(w1),@most_significant_digit(w2)) =
%                                                 if(==(overflow,@zero_word),
%                                                       @most_significant_digit(mult),
%                                                       @concat_digit(@most_significant_digit(overflow),mult))
%                                                 whr mult=@times_word(w1,w2),
%                                                     overflow=@times_overflow_word(w1,w2) end;
 
     *(@most_significant_digit(w1),@most_significant_digit(w2)) =
                @times_whr_mult_overflow(@times_word(w1,w2),@times_overflow_word(w1,w2));
     *(@most_significant_digit(w1),@concat_digit(p2,w2)) = @concat_digit(
                                                                  @times_overflow(p2,w1,@times_overflow_word(w1,w2)),
                                                                  @times_word(w1,w2));

     *(@concat_digit(p1,w1),@most_significant_digit(w2)) = 
                                     @concat_digit(
                                           @times_overflow(p1,w2,@times_overflow_word(w1,w2)),
                                           @times_word(w1,w2));

     *(@concat_digit(p1,w1),@concat_digit(p2,w2)) = 
                              if(<(p1, p2),
                                     @times_ordered(@concat_digit(p1,w1),@concat_digit(p2,w2)),
                                     @times_ordered(@concat_digit(p2,w2),@concat_digit(p1,w1)));

% The following case is not needed as the second argument of @times_ordered has always more than one digit. 
%     @times_ordered(@most_significant_digit(w1),@most_significant_digit(w2)) =
%                @times_whr_mult_overflow(@times_word(w1,w2),@times_overflow_word(w1,w2));

     @times_ordered(@most_significant_digit(w1),@concat_digit(p2,w2)) = 
                                                           @concat_digit(
                                                                  @times_overflow(p2,w1,@times_overflow_word(w1,w2)),
                                                                  @times_word(w1,w2));

     @times_ordered(@concat_digit(p1,w1),p2) = +(@concat_digit(@times_ordered(p1,p2),@zero_word), @times_overflow(p2,w1,@zero_word));

     @times_whr_mult_overflow(w1,w2) = if(@equals_zero_word(w2),
                                              @most_significant_digit(w1),
                                              @concat_digit(@most_significant_digit(w2),w1));
 
     @times_overflow(@most_significant_digit(w1),w2,overflow) =
                @times_whr_mult_overflow(@times_with_carry_word(w1,w2,overflow),@times_with_carry_overflow_word(w1,w2,overflow));

     @times_overflow(@concat_digit(p1,w1),w2,overflow) = 
                                              @concat_digit(
                                                         @times_overflow(p1,w2,@times_with_carry_overflow_word(w1,w2,overflow)),
                                                         @times_with_carry_word(w1,w2,overflow));


