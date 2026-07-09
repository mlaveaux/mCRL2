# lts library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


A review of the `mcrl2_lts` library (`libraries/lts/`). The library is large
(~52k lines); the core data structures, all I/O, the builders/utilities, and the
trace and dispatch layers were reviewed in depth, and the large bisimulation
algorithm files were skimmed for structural concerns.

## Thread-safety (multithreading is ON by default)

### Non-RAII locking → deadlock on exception

- [libraries/lts/include/mcrl2/lts/lts_builder.h](libraries/lts/include/mcrl2/lts/lts_builder.h#L99)
- [libraries/lts/include/mcrl2/lts/stochastic_lts_builder.h](libraries/lts/include/mcrl2/lts/stochastic_lts_builder.h#L141)
- [libraries/lts/include/mcrl2/lts/state_space_generator.h](libraries/lts/include/mcrl2/lts/state_space_generator.h#L443)

All the builders do `m_exclusive_transition_access.lock(); … work …;
m_exclusive_transition_access.unlock();` manually. If the work between (e.g.
`add_action`, `lps::pp`, `out <<`, allocation) throws, the mutex is never released
→ deadlock. Replace with `std::lock_guard`/`std::scoped_lock`.

### Static algorithm state (non-reentrant)

- [libraries/lts/include/mcrl2/lts/detail/liblts_bisim_gjkw.h](libraries/lts/include/mcrl2/lts/detail/liblts_bisim_gjkw.h#L463):
  `static state_type nr_of_blocks;` and `static …_ptr s_i_begin/s_i_end/perm_begin;`

These make the GJKW algorithm non-reentrant and unsafe to run concurrently or
nested. If they are debug-only, isolate them under `#ifndef NDEBUG`; otherwise they
are shared mutable state. (The `simple_list` pool is now a function-local static,
which is fine.)

## Design / modernisation / minor

- [libraries/lts/source/liblts_lts.cpp](libraries/lts/source/liblts_lts.cpp#L269):
  `reinterpret_cast<const state_label_lts&>(term)` on an aterm — project guidelines
  say use `atermpp::down_cast`, never `reinterpret_cast`, on aterm types. The same
  function also has a dead inner `if (filename.empty())` nested inside an
  `if (!filename.empty())` block.
- [transition.h](libraries/lts/include/mcrl2/lts/transition.h#L135) and
  [detail/transition.h](libraries/lts/include/mcrl2/lts/detail/transition.h): the
  hand-written `operator<`/`==`/`!=` are good candidates for C++20 `operator<=>` +
  defaulted `==` (the codebase is actively adopting this). The comparator functors'
  `operator()` should be `const`. The `detail/transition.h` `#endif`
  ([detail/transition.h](libraries/lts/include/mcrl2/lts/detail/transition.h#L192))
  is missing its trailing guard comment and the file lacks the `\file`/`\brief`
  preamble required by the project guidelines;
  [state_label_empty.h](libraries/lts/include/mcrl2/lts/state_label_empty.h#L69)
  `#endif` also lacks the comment.
- Dead/leftover code: unreachable `return;` after `throw` in `save()` in
  [liblts_fsm.cpp](libraries/lts/source/liblts_fsm.cpp#L191).
- [liblts_dot.cpp](libraries/lts/source/liblts_dot.cpp#L56): DOT output does not
  escape `"` in state/action labels, producing malformed DOT if a label contains a
  quote.
- [lts_builder.h](libraries/lts/include/mcrl2/lts/lts_builder.h#L171):
  `lts_aut_disk_builder::finalize` rewrites a fixed-width dummy header via
  `seekp(0)`; for extremely large transition/state counts the formatted header can
  exceed the reserved padding and overwrite the first transition line.

## Overall assessment

The library is well-structured and the newer code (`liblts_lts.cpp`, the builders)
makes good use of modern C++ (`constexpr if`, `std::optional`, structured
bindings, RAII streams). The main remaining risk is **thread-safety**: several
"avoid re-allocating" statics in the bisimulation algorithms and the manual
lock/unlock patterns in the builders are unsafe under the default multithreaded
exploration. Hardening the hand-written `.aut` stream parser (input-size
validation, `unsigned char` for ctype functions) remains worthwhile as well.

## Resolved since previous review (removed)

Verified fixed in the current source and dropped from this document: the
`probabilistic_state` hash/equality contract violation (the length-1 distribution
case now hashes via the same `hash_combine` path as a single state); the
out-of-bounds `std::array` access in `liblts.cpp` (`string_for_type`/
`extension_for_type`/`mime_type_for_type` now use `.at()`, `type_desc_strings` has
the matching 5 entries, and `supported_lts_formats_text` uses a `remaining == 2`
distance check instead of `types.end() - 2`); the `simple_list` static `my_pool`
(now a function-local static); the pointless `try { load } catch(...) { throw; }`
in the `trace.h` constructor; the `lts::clear_actions()` doc/behaviour mismatch
(comment now matches the reset-to-one-tau behaviour); the uninitialised
`m_init_state` (now `= 0`); the `// YYYYYY TODO FINISH.` marker in `lts.h`; and the
unreachable `return;` in `liblts_aut.cpp::save()`.
