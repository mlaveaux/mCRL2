# gui library

_Part of the [mCRL2 Library Reviews](LIBRARY_REVIEWS.md); see the index for the consolidated priority list._


The `gui` library is the small Qt6 helper layer shared by the graphical tools
(`mcrl2ide`, `mcrl2xi`, `mcrl2-gui`, `lpsxsim`, `ltsgraph`, `diagraphica`). It is
not a verification component: it provides the `qt_tool` base that wraps a tool in a
`QApplication`, a syntax-highlighting `CodeEditor` for mCRL2 specifications and
mu-calculus formulae, a `LogWidget`/`LogRelay` that pipes the toolset logger into a
Qt text view, an `ExtendedTabWidget`, a directory-remembering `PersistentFileDialog`,
a handful of `Setting*` value wrappers, and the arcball mouse-rotation helper for the
OpenGL views. The whole library (~1.7k lines across nine headers and seven sources)
is gated behind `MCRL2_ENABLE_GUI_TOOLS` and was read file-by-file. It contains **no
untrusted-input deserialisation** — unlike `atermpp`, `lts`, `pg` or `symbolic`,
nothing here parses a hostile binary file — so there is no security section; the
findings are a latent crash, a resource accumulation bug, and a cluster of Qt-idiom
and style issues.

## 1. Correctness / robustness

### 1.1 Syntax highlighters accumulate on every palette change

[codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L251)

`changeHighlightingRules()` allocates a fresh highlighter and overwrites the member
without disposing of the previous one:

```cpp
void CodeEditor::changeHighlightingRules()
{
  highlighter = new CodeHighlighter(isSpecificationEditor, lightPalette,
                                    this->document());
}
```

The `CodeHighlighter` is constructed with `this->document()` as its parent, so each
instance installs itself on the document and stays connected to its
`contentsChange` signal. The method is called once from `setPurpose()` and then
**again on every `QEvent::PaletteChange`**
([:442](libraries/gui/source/codeeditor.cpp#L442)). Every theme/palette switch
therefore leaves the old highlighter alive and attached: after *n* palette changes
the document has *n+1* highlighters all re-highlighting on every keystroke. They are
parented to the document so they are eventually freed when the editor is destroyed
(no permanent leak), but during the editor's lifetime this is unbounded growth of
redundant highlighting work and a latent source of flicker. The fix is to `delete
highlighter;` (or reuse a single instance and only rebuild its rules) before
reassigning.

### 1.2 `~CodeEditor` dereferences a possibly-null highlighter

[codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L218)

```cpp
CodeEditor::~CodeEditor()
{
  lineNumberArea->deleteLater();
  highlighter->deleteLater();   // highlighter may still be nullptr
}
```

`highlighter` is initialised to `nullptr`
([codeeditor.h](libraries/gui/include/mcrl2/gui/codeeditor.h#L215)) and is only
assigned inside `changeHighlightingRules()`, which the constructor does **not**
call. A `CodeEditor` that is constructed and destroyed without `setPurpose()` ever
being invoked (and without a palette-change event arriving) destructs through a null
`QObject::deleteLater()` — undefined behaviour / crash. Every current call site
happens to call `setPurpose()` immediately after construction
([mcrl2xi documentmanager.cpp](tools/release/mcrl2xi/documentmanager.cpp#L36),
[mcrl2ide mainwindow.cpp](tools/release/mcrl2ide/mainwindow.cpp#L27),
[addeditpropertydialog.cpp](tools/release/mcrl2ide/addeditpropertydialog.cpp#L44)),
so the bug does not trigger today, but the class is public library API and its
destructor is not self-consistent. Both `deleteLater()` calls are in any case
redundant: `lineNumberArea` is parented to `this` and `highlighter` to the document,
so Qt's parent-child ownership already reclaims them. Guarding with `if (highlighter)`
— or simply removing both lines — fixes it.

## 2. Maintainability / style

- **`PersistentFileDialog` inherits `QObject` privately.**
  [persistentfiledialog.h](libraries/gui/include/mcrl2/gui/persistentfiledialog.h#L18)
  declares `class PersistentFileDialog : QObject` (no access specifier ⇒ private
  inheritance) while carrying the `Q_OBJECT` macro. A `QObject` subclass is expected
  to inherit publicly; private inheritance hides the *is-a* relationship so the type
  cannot be passed as a `QObject*`, used with `qobject_cast`, or have its
  signals/slots connected from outside. It compiles and works only because the class
  is used purely through its own concrete methods and just needs the parent pointer
  for lifetime management. Should be `public QObject`.
- **Broken `<br>` in the About box.**
  [qt_tool.h](libraries/gui/include/mcrl2/gui/qt_tool.h#L88) builds the about-dialog
  HTML with `message += "<\br>";`
  ([qt_tool.h](libraries/gui/include/mcrl2/gui/qt_tool.h#L90)). `\b` is the C++
  escape for backspace (0x08), so
  the string is `<`, backspace, `r>` rather than a `<br>` line break. Cosmetic, but
  it injects a stray control character into the rendered HTML.
- **Header guards name the wrong library.** Several headers that now live in `gui`
  still use `MCRL2_UTILITIES_*` guards:
  [setting.h](libraries/gui/include/mcrl2/gui/setting.h#L9) (`MCRL2_UTILITIES_SETTING_H`),
  [extendedtabwidget.h](libraries/gui/include/mcrl2/gui/extendedtabwidget.h#L10),
  [persistentfiledialog.h](libraries/gui/include/mcrl2/gui/persistentfiledialog.h#L10),
  [logwidget.h](libraries/gui/include/mcrl2/gui/logwidget.h#L9), and
  [arcball.h](libraries/gui/include/mcrl2/gui/arcball.h#L9)
  (`MCRL2_UTILITIES_ARCBALL_H`). They should be `MCRL2_GUI_*` to match the directory,
  as `codeeditor.h`, `utilities.h` and `glu.h` already are. The same trailing-guard
  inconsistency was flagged in the `utilities` review.
- **Context-menu actions are dispatched by display text.**
  [extendedtabwidget.cpp](libraries/gui/source/extendedtabwidget.cpp#L54) decides what
  to do by comparing the triggered action's `text()` against English string literals
  (`act->text() == "Close tab"`, `"Close all tabs"`, …). This couples behaviour to the
  exact label and breaks the moment the strings are translated or reworded; storing an
  enum/`data()` on each `QAction`, or comparing pointers, is the robust idiom.
- **Legacy fixed-function OpenGL.**
  [arcball.cpp](libraries/gui/source/arcball.cpp#L77) drives the rotation through
  `glRotatef`/`glGetIntegerv` and [glu.h](libraries/gui/include/mcrl2/gui/glu.h#L1)
  pulls in the deprecated GLU headers (`<GL/glu.h>` / `<OpenGL/glu.h>`). This is the
  immediate-mode pipeline that has been removed from OpenGL core profiles and from
  recent macOS; it is pre-existing tech debt shared with the older OpenGL tools, not a
  defect, but it is a portability liability worth tracking.
- **Unrecognised CMake parameter.**
  [CMakeLists.txt](libraries/gui/CMakeLists.txt#L4) passes `INSTALL_HEADERS TRUE` to
  `mcrl2_add_library`, which does not declare that keyword
  ([cmake/MCRL2AddTarget.cmake](cmake/MCRL2AddTarget.cmake#L7) only parses
  `EXCLUDE_HEADERTEST DPARSER_SOURCES INCLUDE_DIRS SOURCES DEPENDS`), so it is silently
  ignored. Harmless today but misleading.
- **Minor.** `setExtraSelections({})` plus the two highlight passes run *inside*
  `CodeEditor::paintEvent`
  ([codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L420)); this mutates document
  state during painting, which is normally discouraged, but it is guarded by the
  `lastCursor` check so it converges rather than looping
  ([codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L383)). `arcballVector` divides the
  mouse coordinates by the viewport width/height without checking for a zero-sized
  viewport, relying on the later `normalized()` to swallow the resulting NaN.

## 3. Things that look good

- **Cross-thread logging is done the right way.**
  [logwidget.cpp](libraries/gui/source/logwidget.cpp#L15) has `LogRelay::output` (an
  `mcrl2::log::output_policy` that the toolset logger can call from any worker thread)
  merely `emit` a Qt signal; the `LogWidget` connects it with the default
  auto/queued connection so the actual GUI update is marshalled onto the GUI thread.
  The relay registers itself in the constructor and **unregisters in the destructor**
  ([:36](libraries/gui/source/logwidget.cpp#L36)), so the logger never calls a
  dangling policy.
- **`QApplication` creation is postponed for headless `--help`.**
  [qt_tool.h](libraries/gui/include/mcrl2/gui/qt_tool.h#L175) deliberately defers
  constructing the `QApplication` from `execute` to `pre_run`, with a comment
  explaining that doing it earlier would make the `--help` text impossible to show on
  a machine with no display server; the construction is wrapped in `try/catch` and
  reported through `mCRL2log`. The Windows `AttachConsole` dance and the custom
  `qInstallMessageHandler` that routes Qt diagnostics to `stderr` are equally careful.
- **The arcball math is numerically guarded.**
  [arcball.cpp](libraries/gui/source/arcball.cpp#L71) clamps the `acos` argument to
  `[-1, 1]` before taking the angle and falls back to a zero vector when the projected
  length is not `std::isnormal`, avoiding the classic NaN-from-rounding rotation glitch.
- **Resources use RAII / smart pointers.** The `QApplication`
  ([qt_tool.h](libraries/gui/include/mcrl2/gui/qt_tool.h#L130)) and the generated
  the `Ui::LogWidget` pimpl ([logwidget.h](libraries/gui/include/mcrl2/gui/logwidget.h#L54))
  are held in `std::unique_ptr`, and the standard context menu in `CodeEditor`
  ([codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L259)) is owned by a
  `std::unique_ptr<QMenu>`.
- **Monospaced tab width is computed honestly.** `CodeEditor::setFontSize`
  ([codeeditor.cpp](libraries/gui/source/codeeditor.cpp#L230)) averages the advance of
  a long repeated string in `double` precision rather than trusting the integer
  `maxWidth`/`averageCharWidth`, with a comment explaining why — a nice attention to a
  fiddly Qt detail.
