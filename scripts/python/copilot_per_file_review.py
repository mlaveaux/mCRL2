#!/usr/bin/env python3

"""
Run per-file C++ reviews with Copilot CLI using the repository's critical review skill.

Modes:
- file: review one file (clean content or diff)
- diff: review changed source/header files from a diff range
- all: review all source/header files in the repository

The script requires the code-review-graph Python package for transitive dependency context.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

CPP_SUFFIXES = {".h", ".hh", ".hpp", ".hxx", ".c", ".cc", ".cpp", ".cxx"}
DEFAULT_SKILL_PATH = Path(".github/skills/mcrl2-critical-cpp20-review/SKILL.md")
DEFAULT_OUTPUT_DIR = Path("review_reports")
DEFAULT_CRG_CLI_CANDIDATES = (".venv/bin/code-review-graph", "code-review-graph")
DEFAULT_CLANGD_MCP_COMMAND = "iflow-mcp-felipeerias-clangd-mcp-server"
DEFAULT_COMPILE_COMMANDS_DIR = "build"


class CommandError(RuntimeError):
    pass


@dataclass
class ReviewTarget:
    path: Path
    mode: str  # "clean" or "diff"
    base_ref: str | None


@dataclass
class CopilotRunResult:
    returncode: int
    stdout: str
    stderr: str


SESSION_LIMIT_PATTERNS = (
    "session limit",
    "session limits",
    "rate limit",
    "quota",
    "usage limit",
    "too many requests",
)


DEFAULT_SAFE_ALLOW_TOOLS = (
    "read",
    "write",
    "shell(cmake:*)",
    "shell(ninja:*)",
    "shell(ctest:*)",
    "shell(make:*)",
)


def run_command(cmd: Sequence[str], cwd: Path, check: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        list(cmd),
        cwd=str(cwd),
        text=True,
        capture_output=True,
        check=False,
    )
    if check and result.returncode != 0:
        raise CommandError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def run_streaming_command(cmd: Sequence[str], cwd: Path) -> CopilotRunResult:
    proc = subprocess.Popen(
        list(cmd),
        cwd=str(cwd),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
    )

    assert proc.stdout is not None
    assert proc.stderr is not None

    stdout_lines: list[str] = []
    stderr_lines: list[str] = []

    def consume(pipe: subprocess.PIPE, sink: list[str], mirror) -> None:
        try:
            for line in iter(pipe.readline, ""):
                sink.append(line)
                mirror.write(line)
                mirror.flush()
        finally:
            pipe.close()

    t_out = threading.Thread(target=consume, args=(proc.stdout, stdout_lines, sys.stdout), daemon=True)
    t_err = threading.Thread(target=consume, args=(proc.stderr, stderr_lines, sys.stderr), daemon=True)
    t_out.start()
    t_err.start()

    returncode = proc.wait()
    t_out.join()
    t_err.join()

    return CopilotRunResult(
        returncode=returncode,
        stdout="".join(stdout_lines),
        stderr="".join(stderr_lines),
    )


def repo_root_from_git(cwd: Path) -> Path:
    out = run_command(["git", "rev-parse", "--show-toplevel"], cwd=cwd)
    return Path(out.stdout.strip()).resolve()


def is_cpp_source_or_header(path: Path) -> bool:
    return path.suffix.lower() in CPP_SUFFIXES


def is_ignored_path(path: Path) -> bool:
    ignored_prefixes = (
        ".git/",
        "build/",
        ".venv/",
        "3rd-party/",
    )
    text = path.as_posix()
    return any(text.startswith(prefix) for prefix in ignored_prefixes)


def read_text_full(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def git_file_diff(repo_root: Path, rel_path: Path, base_ref: str) -> str:
    out = run_command(
        ["git", "diff", "--no-color", f"{base_ref}", "--", rel_path.as_posix()],
        cwd=repo_root,
    )
    return out.stdout


def list_diff_files(repo_root: Path, base_ref: str) -> list[Path]:
    out = run_command(
        ["git", "diff", "--name-only", "--diff-filter=ACMRTUXB", f"{base_ref}", "--"],
        cwd=repo_root,
    )
    files: list[Path] = []
    for line in out.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        rel = Path(line)
        if is_cpp_source_or_header(rel) and not is_ignored_path(rel):
            files.append(rel)
    return sorted(set(files))


def list_all_cpp_files(repo_root: Path) -> list[Path]:
    files: list[Path] = []
    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(repo_root)
        if is_ignored_path(rel):
            continue
        if is_cpp_source_or_header(rel):
            files.append(rel)
    return sorted(files)


def transitive_dependencies(
    repo_root: Path,
    rel_path: Path,
    max_depth: int,
) -> list[Path]:
    """Get transitive impacted files for rel_path using code_review_graph."""
    try:
        from code_review_graph.tools.query import get_impact_radius  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "The code-review-graph package is required. Install dependencies from "
            "scripts/python/review_requirements.txt."
        ) from exc
    result = get_impact_radius(
        changed_files=[rel_path.as_posix()],
        max_depth=max_depth,
        repo_root=str(repo_root),
        detail_level="standard",
    )

    impacted = result.get("impacted_files", []) if isinstance(result, dict) else []
    normalized: set[Path] = set()
    for item in impacted:
        p = Path(str(item))
        if p.is_absolute():
            try:
                p = p.resolve().relative_to(repo_root)
            except Exception:
                continue
        if p != rel_path:
            normalized.add(p)
    return sorted(normalized)


def resolve_crg_cli(repo_root: Path, explicit_path: str | None) -> str:
    if explicit_path:
        candidate = Path(explicit_path)
        if candidate.is_absolute():
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
        else:
            rel = (repo_root / candidate).resolve()
            if rel.is_file() and os.access(rel, os.X_OK):
                return str(rel)
            found = shutil.which(explicit_path)
            if found:
                return found
        raise FileNotFoundError(f"code-review-graph CLI not found/executable: {explicit_path}")

    for candidate in DEFAULT_CRG_CLI_CANDIDATES:
        p = Path(candidate)
        if p.is_absolute() or "/" in candidate:
            resolved = (repo_root / p).resolve() if not p.is_absolute() else p
            if resolved.is_file() and os.access(resolved, os.X_OK):
                return str(resolved)
        else:
            found = shutil.which(candidate)
            if found:
                return found

    raise FileNotFoundError(
        "Could not find code-review-graph CLI. Install dependencies from "
        "scripts/python/review_requirements.txt."
    )


def build_or_update_dependency_graph_via_cli(
    repo_root: Path,
    graph_cli: str,
    full_rebuild: bool,
    base_ref: str,
) -> dict[str, object]:
    if full_rebuild:
        cmd = [graph_cli, "build", "--repo", str(repo_root)]
    else:
        cmd = [graph_cli, "update", "--repo", str(repo_root), "--base", base_ref]

    out = run_command(cmd, cwd=repo_root)
    return {
        "status": "ok",
        "command": cmd,
        "stdout": out.stdout,
    }


def initialize_mcp_config(
    output_dir: Path,
    repo_root: Path,
    graph_cli: str,
    include_clangd: bool,
    clangd_command: str,
    clangd_args: Sequence[str],
    clangd_env: dict[str, str],
    config_name: str,
) -> Path:
    servers: dict[str, dict[str, object]] = {
        "code-review-graph": {
            "type": "stdio",
            "command": graph_cli,
            "args": ["serve", "--repo", str(repo_root)],
            "cwd": str(repo_root),
        }
    }

    if include_clangd:
        clangd_entry: dict[str, object] = {
            "type": "stdio",
            "command": clangd_command,
            "args": list(clangd_args),
            "cwd": str(repo_root),
        }
        if clangd_env:
            clangd_entry["env"] = clangd_env
        servers["clangd"] = clangd_entry

    mcp_config = {
        "mcpServers": servers,
    }

    config_path = output_dir / config_name
    config_path.write_text(json.dumps(mcp_config, indent=2) + "\n", encoding="utf-8")
    return config_path


def resolve_clangd_env(
    repo_root: Path,
    compile_commands_dir: str,
    clangd_path: str | None,
) -> dict[str, str]:
    """Build environment for the clangd MCP server.

    The felipeerias clangd MCP server reads PROJECT_ROOT, COMPILE_COMMANDS_DIR
    and CLANGD_PATH from its environment.
    """
    ccdir = Path(compile_commands_dir)
    if not ccdir.is_absolute():
        ccdir = (repo_root / ccdir).resolve()

    env: dict[str, str] = {
        "PROJECT_ROOT": str(repo_root),
        "COMPILE_COMMANDS_DIR": str(ccdir),
    }

    resolved_clangd = clangd_path or shutil.which("clangd")
    if resolved_clangd:
        env["CLANGD_PATH"] = resolved_clangd

    return env


def load_skill_text(repo_root: Path, skill_path: Path) -> str:
    path = skill_path if skill_path.is_absolute() else (repo_root / skill_path)
    if not path.is_file():
        raise FileNotFoundError(f"Skill file not found: {path}")
    return path.read_text(encoding="utf-8", errors="replace")


def build_prompt(
    *,
    skill_text: str,
    rel_path: Path,
    review_input_mode: str,
    content: str,
    dependencies: Sequence[Path],
    thinking: str,
    effort: str,
) -> str:
    dep_lines = "\n".join(f"- {p.as_posix()}" for p in dependencies) if dependencies else "- (none found)"

    return (
        "You are reviewing mCRL2 code using the project review skill.\n"
        "Follow the exact markdown output requirements from the skill and be critical and skeptical.\n\n"
        f"Model guidance:\n- thinking mode: {thinking}\n- reasoning effort: {effort}\n\n"
        "Skill text:\n"
        "```markdown\n"
        f"{skill_text}\n"
        "```\n\n"
        f"Target file: {rel_path.as_posix()}\n"
        f"Input mode: {review_input_mode}\n\n"
        "Transitive dependency context (from code-review-graph):\n"
        f"{dep_lines}\n\n"
        "Review input:\n"
        "```\n"
        f"{content}\n"
        "```\n"
    )


def contains_session_limit(text: str) -> bool:
    lower = text.lower()
    return any(pattern in lower for pattern in SESSION_LIMIT_PATTERNS)


def extract_markdown_review(stdout: str) -> str:
    text = stdout.strip()
    if not text:
        return "No markdown review content produced by Copilot."

    match = re.search(r"(^#\s+Review\s+Findings.*)", text, flags=re.IGNORECASE | re.DOTALL | re.MULTILINE)
    if match:
        return match.group(1).strip()
    return text


def build_copilot_command_prefix(args: argparse.Namespace) -> list[str]:
    cmd: list[str] = [args.copilot_bin]
    cmd.extend(["--model", args.model])
    cmd.extend(["--stream", "on"])

    if args.copilot_log_level:
        cmd.extend(["--log-level", args.copilot_log_level])

    if args.approval_mode == "all":
        cmd.append("--allow-all")
    elif args.approval_mode == "safe-defaults":
        for tool in DEFAULT_SAFE_ALLOW_TOOLS:
            cmd.extend(["--allow-tool", tool])
    
    for tool in args.allow_tool:
        cmd.extend(["--allow-tool", tool])

    if getattr(args, "copilot_mcp_config", None):
        cmd.extend(["--additional-mcp-config", f"@{args.copilot_mcp_config}"])

    return cmd


def collect_targets(args: argparse.Namespace, repo_root: Path) -> list[ReviewTarget]:
    if args.command == "file":
        rel = Path(args.file)
        if rel.is_absolute():
            rel = rel.resolve().relative_to(repo_root)
        if not (repo_root / rel).is_file():
            raise FileNotFoundError(f"File not found: {rel}")
        if not is_cpp_source_or_header(rel):
            raise ValueError(f"Not a source/header file: {rel}")
        return [ReviewTarget(rel, args.input_mode, args.base_ref)]

    if args.command == "diff":
        files = list_diff_files(repo_root, args.base_ref)
        return [ReviewTarget(path=f, mode="diff", base_ref=args.base_ref) for f in files]

    if args.command == "all":
        files = list_all_cpp_files(repo_root)
        return [ReviewTarget(path=f, mode="clean", base_ref=None) for f in files]

    raise ValueError(f"Unknown command: {args.command}")


def output_path(output_dir: Path, rel_path: Path, input_mode: str) -> Path:
    safe = rel_path.as_posix().replace("/", "__")
    return output_dir / f"{safe}.{input_mode}.review.md"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run per-file C++/header reviews via Copilot CLI using mCRL2 review skill.",
    )

    parser.add_argument(
        "--copilot-bin",
        default="copilot",
        help="Copilot CLI executable (default: copilot)",
    )
    parser.add_argument(
        "--model",
        default="gpt-5.3-codex",
        help="Copilot model passed to --model (default: gpt-5.3-codex)",
    )
    parser.add_argument(
        "--thinking",
        choices=["off", "on", "deep"],
        default="on",
        help="Thinking hint injected into the review prompt.",
    )
    parser.add_argument(
        "--effort",
        choices=["low", "medium", "high"],
        default="high",
        help="Reasoning effort hint injected into the review prompt.",
    )
    parser.add_argument(
        "--approval-mode",
        choices=["safe-defaults", "all", "none"],
        default="safe-defaults",
        help=(
            "Copilot tool approval profile. safe-defaults allows read/write and "
            "common build tools (cmake/ninja/ctest/make); all uses --allow-all."
        ),
    )
    parser.add_argument(
        "--allow-tool",
        action="append",
        default=[],
        help="Additional --allow-tool entries passed to Copilot CLI (repeatable).",
    )
    parser.add_argument(
        "--graph-cli-bin",
        default=None,
        help=(
            "Path or executable name for the code-review-graph CLI. "
            "Defaults to .venv/bin/code-review-graph or PATH lookup."
        ),
    )
    parser.add_argument(
        "--mcp-config-name",
        default="copilot.mcp.json",
        help="Filename for generated Copilot MCP config inside output directory.",
    )
    parser.add_argument(
        "--disable-clangd-mcp",
        action="store_true",
        help="Do not include clangd MCP server in generated Copilot MCP config.",
    )
    parser.add_argument(
        "--clangd-mcp-command",
        default=DEFAULT_CLANGD_MCP_COMMAND,
        help=f"Command for clangd MCP server (default: {DEFAULT_CLANGD_MCP_COMMAND}).",
    )
    parser.add_argument(
        "--clangd-mcp-arg",
        action="append",
        default=[],
        help="Additional arg for clangd MCP server (repeatable).",
    )
    parser.add_argument(
        "--compile-commands-dir",
        default=DEFAULT_COMPILE_COMMANDS_DIR,
        help=(
            "Directory containing compile_commands.json for the clangd MCP "
            f"server (default: {DEFAULT_COMPILE_COMMANDS_DIR})."
        ),
    )
    parser.add_argument(
        "--clangd-path",
        default=None,
        help="Path to the clangd binary for the clangd MCP server (default: PATH lookup).",
    )
    parser.add_argument(
        "--copilot-log-level",
        choices=["none", "error", "warning", "info", "debug", "all", "default"],
        default="info",
        help="Copilot CLI log level for verbose terminal output.",
    )
    parser.add_argument(
        "--skill-path",
        type=Path,
        default=DEFAULT_SKILL_PATH,
        help="Path to the review skill markdown file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for markdown review outputs.",
    )
    parser.add_argument(
        "--dependency-depth",
        type=int,
        default=2,
        help="Transitive impact depth used by code-review-graph (default: 2).",
    )
    parser.add_argument(
        "--graph-full-rebuild",
        action="store_true",
        help="Force a full code-review-graph rebuild before reviews.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit number of reviewed files (0 means no limit).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="Skip files that already have a markdown review output (default: true).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run reviews even when output markdown already exists.",
    )
    parser.add_argument(
        "--stop-on-session-limit",
        action="store_true",
        default=True,
        help="Stop processing remaining files when Copilot reports session/rate limits.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected targets and exit without calling Copilot CLI.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    file_parser = subparsers.add_parser("file", help="Review a single file")
    file_parser.add_argument("file", help="Path to the file (repo-relative or absolute)")
    file_parser.add_argument(
        "--input-mode",
        choices=["clean", "diff"],
        default="clean",
        help="Review clean file contents or git diff.",
    )
    file_parser.add_argument(
        "--base-ref",
        default="HEAD",
        help="Git ref/range for diff collection in file mode (default: HEAD).",
    )

    diff_parser = subparsers.add_parser("diff", help="Review changed files from git diff")
    diff_parser.add_argument(
        "--base-ref",
        default="master...HEAD",
        help="Git diff reference passed to `git diff --name-only <base-ref>`.",
    )

    subparsers.add_parser("all", help="Review all source/header files")

    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    repo_root = repo_root_from_git(Path.cwd())

    # Quick validation of copilot CLI before heavy work.
    run_command([args.copilot_bin, "--version"], cwd=repo_root)

    skill_text = load_skill_text(repo_root, args.skill_path)
    targets = collect_targets(args, repo_root)

    if not targets:
        print("No files selected for review.")
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Repository: {repo_root}")
    print(f"Targets: {len(targets)}")
    print(f"Model: {args.model}")
    print(f"Thinking: {args.thinking} | Effort: {args.effort}")
    print(f"Approval mode: {args.approval_mode}")

    graph_base = args.base_ref if hasattr(args, "base_ref") else "HEAD~1"
    graph_cli = resolve_crg_cli(repo_root, args.graph_cli_bin)
    graph_result = build_or_update_dependency_graph_via_cli(
        repo_root=repo_root,
        graph_cli=graph_cli,
        full_rebuild=args.graph_full_rebuild,
        base_ref=graph_base,
    )
    print(f"code-review-graph update: {graph_result.get('status', 'ok')}")

    include_clangd = not args.disable_clangd_mcp
    clangd_env: dict[str, str] = {}
    if include_clangd:
        clangd_env = resolve_clangd_env(
            repo_root=repo_root,
            compile_commands_dir=args.compile_commands_dir,
            clangd_path=args.clangd_path,
        )
        ccjson = Path(clangd_env["COMPILE_COMMANDS_DIR"]) / "compile_commands.json"
        if not ccjson.is_file():
            print(
                f"  warning: {ccjson} not found; configure a CMake build with "
                "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON so clangd MCP works.",
                file=sys.stderr,
            )

    mcp_config = initialize_mcp_config(
        output_dir=args.output_dir,
        repo_root=repo_root,
        graph_cli=graph_cli,
        include_clangd=include_clangd,
        clangd_command=args.clangd_mcp_command,
        clangd_args=args.clangd_mcp_arg,
        clangd_env=clangd_env,
        config_name=args.mcp_config_name,
    )
    args.copilot_mcp_config = mcp_config
    print(f"Copilot MCP config: {mcp_config}")

    if args.dry_run:
        shown = 0
        for target in targets:
            out_path = output_path(args.output_dir, target.path, target.mode)
            if args.skip_existing and (not args.force) and out_path.exists():
                continue
            print(f"- {target.path.as_posix()} ({target.mode})")
            shown += 1
            if args.limit > 0 and shown >= args.limit:
                break
        return 0

    failures = 0
    reviewed = 0
    session_limited = False
    for index, target in enumerate(targets, start=1):
        abs_path = repo_root / target.path
        output_mode = target.mode
        out_path = output_path(args.output_dir, target.path, output_mode)

        if args.skip_existing and (not args.force) and out_path.exists():
            print(f"[{index}/{len(targets)}] Skipping {target.path.as_posix()} (already reviewed)")
            continue

        print(f"[{index}/{len(targets)}] Reviewing {target.path.as_posix()} ({target.mode})")

        if target.mode == "clean":
            review_input = read_text_full(abs_path)
        else:
            diff_text = git_file_diff(repo_root, target.path, target.base_ref or "HEAD")
            if not diff_text.strip():
                review_input = read_text_full(abs_path)
                output_mode = "clean-fallback"
            else:
                review_input = diff_text
                output_mode = "diff"

        try:
            deps = transitive_dependencies(
                repo_root,
                target.path,
                max_depth=args.dependency_depth,
            )
        except Exception as exc:
            print(f"  dependency analysis failed: {exc}", file=sys.stderr)
            failures += 1
            continue

        prompt = build_prompt(
            skill_text=skill_text,
            rel_path=target.path,
            review_input_mode=target.mode,
            content=review_input,
            dependencies=deps,
            thinking=args.thinking,
            effort=args.effort,
        )

        cmd_prefix = build_copilot_command_prefix(args)
        cmd = cmd_prefix + ["-p", prompt]

        try:
            response = run_streaming_command(cmd, cwd=repo_root)
        except Exception as exc:
            print(f"  review failed: {exc}", file=sys.stderr)
            failures += 1
            continue

        combined_output = (response.stdout or "") + "\n" + (response.stderr or "")
        limit_hit = contains_session_limit(combined_output)

        review_markdown = extract_markdown_review(response.stdout)

        timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
        header = {
            "file": target.path.as_posix(),
            "mode": target.mode,
            "output_mode": output_mode,
            "base_ref": target.base_ref,
            "model": args.model,
            "thinking": args.thinking,
            "effort": args.effort,
            "generated_at_utc": timestamp,
            "dependency_count": len(deps),
            "copilot_return_code": response.returncode,
            "session_limit_detected": limit_hit,
        }

        out_path = output_path(args.output_dir, target.path, output_mode)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        report = (
            "---\n"
            f"{json.dumps(header, indent=2)}\n"
            "---\n\n"
            "# File Review\n\n"
            f"{review_markdown}\n\n"
            "## Copilot Raw Stdout\n\n"
            "```text\n"
            f"{response.stdout}"
            "\n```\n\n"
            "## Copilot Raw Stderr\n\n"
            "```text\n"
            f"{response.stderr}"
            "\n```\n"
        )
        out_path.write_text(
            report,
            encoding="utf-8",
        )
        print(f"  wrote {out_path}")
        reviewed += 1

        if response.returncode != 0:
            failures += 1
            print(
                f"  copilot invocation failed with exit code {response.returncode}",
                file=sys.stderr,
            )

        if limit_hit and args.stop_on_session_limit:
            print(
                "  Copilot session/rate limit detected. Stopping now so a later run can resume from remaining files.",
                file=sys.stderr,
            )
            session_limited = True
            break

        if args.limit > 0 and reviewed >= args.limit:
            print(f"  Reached review limit of {args.limit} file(s). Stopping.")
            break

    if session_limited:
        return 2

    if failures:
        print(f"Finished with {failures} failure(s).", file=sys.stderr)
        return 1

    print("Finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
