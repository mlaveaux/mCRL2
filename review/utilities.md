# utilities library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `utilities` library is the lowest-level support layer of the toolset: the
tool/CLI framework, threading primitives, hand-rolled containers
(`indexed_set`, `unordered_set/map`, `hashtable`, `fixed_size_cache`), custom
allocators (`block_allocator`, `detail/bucket_list`, `detail/free_list`),
arbitrary-precision arithmetic (`big_numbers.h`), logging, and string/file/number
helpers. It is depended on by every other library, so correctness and
thread-safety here have toolset-wide impact.

## 1. Correctness

### 1.5 `hashtable::find` is O(n) — *performance bug*

`include/mcrl2/utilities/detail/hashtable.h` (`find()` body)

Unlike `insert`/`erase` (which probe from `get_index(key)`), `find` does a **full
linear scan** of the whole table, defeating the point of a hash table:

```cpp
for (auto it = begin(); it != end(); ++it)
  if (*it != nullptr && m_equals(*it, key)) return it;
return end();
```

The equality functor is now used correctly (`m_equals`, previously a raw `==`), so
the table no longer disagrees with itself — but the lookup is still `O(n)`.
**Fix:** mirror `erase` and probe from `get_index(key)`, stopping at the first
`nullptr`.

## 2. Thread-safety / hygiene

*(No open items — see "Resolved since previous review".)*

## 3. Performance

### 3.1 `indexed_set` carries a dead per-instance heap allocation

`include/mcrl2/utilities/indexed_set.h:37`

`mutable std::shared_ptr<std::mutex> m_mutex;` is allocated in every constructor
(`detail/indexed_set.h`: `m_mutex(new std::mutex())`) but **never used** — all
locking goes through `m_shared_mutexes`. Every `indexed_set` pays a heap
allocation and an atomic refcount for nothing. **Fix:** delete the member.

### 3.3 `fixed_size_cache` FIFO policy stores a second copy of every key — *memory note*

`include/mcrl2/utilities/cache_policy.h`

`fifo_policy` keeps a `std::forward_list<key_type> m_queue` of every key in
addition to the map, doubling key memory. For the common case (LRU/FIFO over
aterm pointers) the key is a single pointer, so this is acceptable; flagged only
as a memory note.

## 4. Maintainability / architecture

### 4.1 `indexed_set`'s `ThreadSafe` template parameter is dead

`include/mcrl2/utilities/indexed_set.h:25`

The class is templated on `bool ThreadSafe = false`, but the value is never used
to gate behaviour — it only appears in the `INDEXED_SET_TEMPLATE`/`INDEXED_SET`
macro parameter lists (`detail/indexed_set.h:49-50`). Locking is unconditional
(gated only at runtime by `GlobalThreadSafe` inside `shared_mutex`). Either honour
the flag (skip the `shared_mutex` machinery when `ThreadSafe == false`) or remove
the parameter.

### 4.2 Three near-identical open-addressing hash tables

`indexed_set`, `hashtable`, and `unordered_set` each re-implement open
addressing / chaining with their own probing, resize, and load-factor logic. They
share the `detail::minimal_hashtable_size`/`PRIME_NUMBER` constants but little
else, and they diverge in quality (1.5 affects only `hashtable`). Sharing at least
the probe loop would reduce the surface for bugs like 1.5.

### 4.3 Wrong `\file` comment in `detail/indexed_set.h`

`include/mcrl2/utilities/detail/indexed_set.h:8` tags the file
`\file utilities/detail/indexed_set.cpp` (wrong extension — it is a header).

### 4.4 `std::aligned_storage_t` is deprecated in C++23

`include/mcrl2/utilities/thread_local.h:81` uses
`std::aligned_storage_t<sizeof(T), alignof(T)>` for the entry storage. It is
deprecated in C++23; the modern form is `alignas(T) std::byte[sizeof(T)]`. The
codebase targets C++20 today but is actively modernising, so this will start
warning on a compiler bump.

### 4.5 Header-guard trailing-underscore inconsistencies

Several headers use a trailing-underscore guard (`..._H_`) against the project
convention of `..._H`: `mutex.h`, `noncopyable.h`, `thread_local.h`,
`shared_reference.h`, `tagged_pointer.h`, `block_allocator.h`,
`detail/atomic_wrapper.h`, `detail/bucket_list.h`, `detail/free_list.h`. None
affect behaviour.

### 4.7 `block_allocator` is not a conforming STL allocator

`include/mcrl2/utilities/block_allocator.h`

The allocator derives from `noncopyable` and defines no `operator==`/`operator!=`,
so it cannot be used with `std::vector`/`std::unordered_map`; it works only with
the in-house `unordered_set`/`bucket_list`. That is a deliberate restriction but
undocumented — a `\details` note ("intended only for the mCRL2 in-house
containers; not a conforming std allocator") would set expectations. Also
`allocate(n, hint)` throws `std::bad_alloc` for any `n != 1`, i.e. it is strictly a
node allocator.

### 4.8 `ceil_log2` is actually `bit_width` (misnamed, off by one on powers of two)

`include/mcrl2/utilities/math.h:25`

`ceil_log2(n)` counts how many right-shifts reduce `n` to zero, i.e. it returns
`floor(log2(n)) + 1` = `std::bit_width(n)`, **not** `ceil(log2(n))`. For exact
powers of two it is one too large (`ceil_log2(8)` returns 4, true `ceil(log2 8)`
is 3). The call sites actually want *bit-width* semantics, so they are correct
today — but the name actively misleads. **Fix:** rename to `bit_width` or use C++20
`std::bit_width`. The file also carries a wrong `\file general_utilities.h` tag
(`math.h:9`).

## 5. Minor / stylistic

- **PascalCase types in `thread_local.h`.** `ThreadLocal`, `Entry`,
  `ThreadBucket`, `Iter` (and `block_allocator.h`'s `ThreadLocalAllocState`) break
  the library's `snake_case` convention (Rust-port heritage).
- **`atomic_wrapper` default ctor assumes a numeric `T`.**
  `detail/atomic_wrapper.h` initialises `std::atomic<T>(0)`, so the default
  constructor only compiles for types constructible from `0`.
- **`unordered_set_iterator : std::iterator_traits<Key>`** (`unordered_set.h:91`)
  privately inherits from `iterator_traits` for no reason (it redeclares all five
  member typedefs itself); the base is inert and can be dropped.
- **`fifo_function_cache` alias takes a single `Args`** (`fixed_size_cache.h`)
  while the underlying `function_cache` is variadic (`Args...`); the alias cannot
  express multi-argument cached functions.
- **`parse_numbers.h` passes raw `char` to `std::isspace`/`std::isdigit`.**
  `parse_next_natural_number`/`parse_natural_number` call the `<cctype>` functions
  with a `char`; on signed-char platforms a byte ≥ 0x80 is a negative `int`, which
  is UB. Cast to `unsigned char` first. The file brief is also the
  `add your file description here.` placeholder (one of many across the library).
- **`free_list`/`bucket_list` `operator Element&()`/`operator Key&()`** are marked
  `explicit` yet documented "implicit conversion".
- **`detail/iota.h`** is a one-line reimplementation of `std::iota` (`<numeric>`);
  prefer the standard algorithm unless there is a reason to avoid the include.
- **`probabilistic_arbitrary_precision_fraction` `operator+`/`operator/` can break
  the `enumerator <= denominator` invariant.** The constructors `assert` the
  fraction is in [0,1] (it models a probability), yet `a + b` (e.g. `1/2 + 3/4`)
  and `a / b` can exceed 1 and assert-fail in debug. Document these operators as
  partial or drop the invariant assert.

## 6. Things that look good

- **The `shared_mutex` busy/forbidden protocol** is a careful, correct asymmetric
  readers–writer lock with explicit reasoning about nesting (`m_lock_depth`).
- **`indexed_set` lock-free fast path.** `put_in_hashtable` uses a
  `compare_exchange_strong` on the slot and a `RESERVED` sentinel so readers spin
  only briefly; the `reserve_indices` batching is a thoughtful contention
  optimisation.
- **`unordered_set` lock-free `emplace`.** The `EnableLockfreeInsertion` path uses
  `emplace_front_unique` on the bucket so concurrent inserts into the same bucket
  are resolved without a mutex; the `static_assert(!(ThreadSafe && Resize))`
  correctly forbids the one combination the design cannot support.
- **`ThreadLocal` per-thread ownership** makes the `relaxed` load of `present`
  safe: only the owning thread ever writes its entry.
- **`fixed_size_cache` find/emplace split** is correctly reasoned (eviction must
  not pick the just-inserted key).
- **`container_utility.h`** is a tidy, dependency-light set of `std::set`/`std::map`
  helpers with the right specialisations.

## Resolved since previous review (removed)

Verified fixed in the current source and dropped from this document:
`1.1` (ThreadLocal iterator no longer mutates on deref), `1.2` (unordered_map
`insert`/`insert_or_assign` overloads), `1.3` (unordered_set post-increment),
`1.4` (hashtable `clear()` resets count/table), `1.6` (indexed_set reverse-iterator
consistency), `1.7` (stack_array `empty()`), `1.8` (free_list compile errors —
now compiled and exercised by `free_list_test.cpp`), `1.9` (block_allocator
`consolidate()` now marks the un-bumped tail as sentinel), `1.10` (bucket_list
`emplace_front_unique` compares the constructed key), `2.1` (interference-size
symbols moved to `hardware_interference_size.h` in `mcrl2::utilities`),
`2.2` (resolved together with 1.1), `2.3` (tagged_pointer `tag`/`untag` reworked
as pure free functions), `2.4` (`big_natural_number::operator<<` now uses a
`thread_local`), `4.3a` (hashtable `#endif` comment), `4.6` (free_list dead/
non-compiling code), and the `stack_array`/`power_of_two` guard names,
`shared_reference.h` double include, `configuration.h` preamble order, and
`big_numbers.h` `isdigit` `unsigned char` cast (minor items).
