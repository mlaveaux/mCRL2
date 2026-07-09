# smt library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `smt` library is the modern bridge to an external SMT solver (Z3): it translates
mCRL2 data specifications, sorts and expressions into SMT-LIB 2, spawns the solver as
a child process, and exchanges problems and answers over pipes. It is distinct from —
and much cleaner than — the experimental in-tree prover reviewed under `data`. The two
source files (`child_process.cpp`, `solver.cpp`) and the translation headers were read
file-by-file.

## 1. Robustness

### 1.1 `read(timeout)` packs the whole microsecond count into `tv_usec`

[child_process.cpp](libraries/smt/source/child_process.cpp#L211)

```cpp
struct timeval tv{};
tv.tv_sec = 0;
tv.tv_usec = timeout.count();   // timeout is std::chrono::microseconds
...
int result = ::select(m_pimpl->pipe_stdout[0] + 1, &readfds, nullptr, nullptr, &tv);
```

`tv_usec` must hold a sub-second value (`< 1'000'000`); any timeout of a second or
more makes it out of range, and POSIX `select` is entitled to fail with `EINVAL`,
which this code turns into a thrown "Error while waiting for input" rather than a
timeout. The fix is the usual split — `tv_sec = count / 1'000'000`, `tv_usec = count %
1'000'000`.

### 1.2 Windows `read()` underflows its index at end-of-stream

[child_process.cpp](libraries/smt/source/child_process.cpp#L60)

```cpp
DWORD dwRead, totalRead = 0;
do {
  bSuccess = ReadFile(m_pimpl->g_hChildStd_OUT_Rd, output + totalRead, 512 - totalRead, &dwRead, NULL);
  totalRead += dwRead;
} while (bSuccess && output[totalRead-1] != '\n');
```

If the first `ReadFile` succeeds with `dwRead == 0` (the child closed the pipe),
`totalRead` is `0` and `output[totalRead-1]` indexes `output[(DWORD)-1]` — a wild
out-of-bounds read. The loop also has no room past 512 bytes: once `totalRead`
reaches 512, `512 - totalRead` is `0` and a solver reply longer than the buffer with
no embedded newline spins reading nothing. (POSIX `read()` uses a bounded
`std::array` and is fine; this is the Windows path only.)

### 1.3 Destructor reaps an arbitrary child

[child_process.cpp](libraries/smt/source/child_process.cpp#L368)

```cpp
child_process::~child_process()
{
  ::close(...);
  int return_status;
  ::wait(&return_status);   // waits for ANY child
}
```

`::wait` reaps whichever child happens to exit first, not necessarily this solver's
`child_pid`. When more than one child process exists (e.g. several solver instances,
or any other forked helper), this can reap the wrong one and leave the actual solver
as a zombie. `::waitpid(m_pimpl->child_pid, …)` targets the right process.

## 2. Security note (calibrated)

`initialize()` launches the solver with `execvp("z3", …)`
([child_process.cpp](libraries/smt/source/child_process.cpp#L273)), a `PATH`-based
lookup. Crucially the SMT *problem* is delivered over the child's **stdin**
(`child_process::write`), never through `argv` or a shell, so there is **no command
injection** surface. The only residual concern is the usual one for invoking an
external tool by bare name: a caller-controlled `PATH` could substitute a malicious
`z3`. This is a minor hardening note, not a vulnerability — the same calibration as
the experimental prover in `data`.

## 3. Maintainability / minor

- **Process-global signal change per init.**
  [child_process.cpp](libraries/smt/source/child_process.cpp#L238) calls
  `signal(SIGPIPE, SIG_IGN)` on every `initialize()`, mutating process-wide signal
  disposition (and `signal` is less portable than `sigaction`).
- **Name/binary case mismatch.** `smt_solver` constructs its `child_process` as
  `z3("Z3")` so diagnostics say "Z3", while the executable launched is `"z3"` — purely
  cosmetic, but worth aligning.
- **Child keeps the original pipe fds.** After `dup2` the forked child closes the
  opposite pipe ends but not the now-redundant originals (no `FD_CLOEXEC`); harmless
  here but slightly untidy fd hygiene.
- A handful of legitimate `// TODO`s mark incomplete SMT translations
  ([solver.cpp](libraries/smt/source/solver.cpp#L233),
  [utilities.h](libraries/smt/include/mcrl2/smt/utilities.h#L136),
  [unfold_pattern_matching.h](libraries/smt/include/mcrl2/smt/unfold_pattern_matching.h#L139)).

## 4. Things that look good

- **Clean answer parsing.** `smt_solver::execute_and_check`
  ([solver.cpp](libraries/smt/source/solver.cpp#L19)) reads into a `std::string` and
  uses `starts_with("sat")`/`"unsat"`/`"unknown"` (correctly ordered, since "unsat"
  does not start with "sat"), throwing a `runtime_error` and logging the problem on
  any unexpected reply — a marked contrast with the experimental `data` prover's
  uninitialised-buffer `strncmp`.
- **No shell, correct fork discipline.** The problem goes over stdin, and the forked
  child calls `_exit(errno)` (not `exit`) after a failed `execvp`, which is the
  correct practice for a fork child; pipe-creation and fork failures throw with the
  already-opened descriptors cleaned up.
- **Tidy production code.** No placeholder file briefs, no stray `std::cout`/`cerr`/
  `clog` (only the examples print), no `system`/`popen`, and the Windows/POSIX
  implementations are cleanly separated behind a `platform_impl` pimpl.
