# (C) Copyright 2024- ECMWF and individual contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Validate and tag a thermofeel release.

Usage (normally via the Makefile)::

    make release X.Y.Z              # plan only: validate and show what would happen
    make release X.Y.Z CONFIRM=1    # apply: bump if needed, commit, create the tag

Environment toggles:
    CONFIRM=1       perform local changes (bump/commit/tag); without it, plan only
    ALLOW_MAJOR=1   permit a MAJOR version increase (guarded per AGENTS.md)
    ALLOW_BRANCH=1  permit running off the 'main' branch (testing only)
    SKIP_GATE=1     skip running `make all` before tagging

The version argument is the single source of truth for the release: the script
checks it against ``thermofeel/__init__.py`` ``__version__`` and bumps that
forward to match if needed -- never backward (SemVer). The script NEVER pushes;
it prints the push commands. Publishing is done by CI
(.github/workflows/cd.yml) when the bare tag is pushed.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path
from typing import NoReturn

INIT_FILE = Path("thermofeel/__init__.py")
CHANGELOG = Path("CHANGELOG.md")
RELEASE_BRANCH = "main"
REMOTE = "origin"

SEMVER = re.compile(r"^\d+\.\d+\.\d+$")
VERSION_RE = re.compile(r"""__version__\s*=\s*["']([^"']+)["']""")


def flag(name: str) -> bool:
    return os.environ.get(name, "").strip() not in ("", "0", "false", "no")


CONFIRM = flag("CONFIRM")
ALLOW_MAJOR = flag("ALLOW_MAJOR")
ALLOW_BRANCH = flag("ALLOW_BRANCH")
SKIP_GATE = flag("SKIP_GATE")


def die(msg: str) -> NoReturn:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(1)


def git(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(["git", *args], text=True, capture_output=True)
    if check and proc.returncode != 0:
        die(f"`git {' '.join(args)}` failed: {proc.stderr.strip()}")
    return proc


def git_out(*args: str) -> str:
    return git(*args).stdout.strip()


def read_version() -> str:
    if not INIT_FILE.is_file():
        die(f"cannot find {INIT_FILE}")
    match = VERSION_RE.search(INIT_FILE.read_text())
    if not match:
        die(f"no __version__ found in {INIT_FILE}")
    return match.group(1)


def parse(version: str) -> tuple[int, int, int]:
    major, minor, micro = (int(part) for part in version.split("."))
    return major, minor, micro


def main() -> int:
    requested = sys.argv[1].strip() if len(sys.argv) > 1 else ""

    print(f"thermofeel release -- target version: {requested or '<none>'}\n")

    if not requested:
        die("no version given. Usage: make release X.Y.Z  (or VERSION=X.Y.Z)")
    if not SEMVER.match(requested):
        die(
            f"'{requested}' is not a bare SemVer MAJOR.MINOR.MICRO "
            "(no leading 'v', no suffix)"
        )
    if not Path(".git").exists() and git("rev-parse", check=False).returncode != 0:
        die("not inside a git repository")

    current = read_version()
    if not SEMVER.match(current):
        die(f"current __version__ '{current}' is not bare SemVer")

    print(f"current __version__: {current}")
    print("checks:")

    blockers: list[str] = []

    def ok(msg: str) -> None:
        print(f"  ok     {msg}")

    def block(msg: str) -> None:
        blockers.append(msg)
        print(f"  BLOCK  {msg}")

    # 1. SemVer relation between requested and current.
    req_t, cur_t = parse(requested), parse(current)
    bump = False
    if req_t == cur_t:
        ok(f"version already at {requested} (no bump needed)")
    elif req_t > cur_t:
        bump = True
        if req_t[0] > cur_t[0] and not ALLOW_MAJOR:
            block(f"MAJOR increase {current} -> {requested} requires ALLOW_MAJOR=1")
        else:
            ok(f"will bump __version__ {current} -> {requested}")
    else:
        block(f"refusing to move version backward: {current} -> {requested}")

    # 2. Tag must not already exist (locally or on the remote).
    git("fetch", "--quiet", "--tags", REMOTE, check=False)
    if (
        git(
            "rev-parse", "-q", "--verify", f"refs/tags/{requested}", check=False
        ).returncode
        == 0
    ):
        block(f"tag '{requested}' already exists locally")
    elif (
        git(
            "ls-remote", "--tags", "--exit-code", REMOTE, requested, check=False
        ).returncode
        == 0
    ):
        block(f"tag '{requested}' already exists on {REMOTE}")
    else:
        ok(f"tag '{requested}' does not exist yet")

    # 3. CHANGELOG has a section for the release.
    changelog_ok = CHANGELOG.is_file() and re.search(
        rf"^##\s+{re.escape(requested)}(\s|$)", CHANGELOG.read_text(), re.MULTILINE
    )
    if changelog_ok:
        ok(f"{CHANGELOG} has a '## {requested}' section")
    else:
        block(f"{CHANGELOG} has no '## {requested}' section")

    # 4. On the release branch.
    branch = git_out("rev-parse", "--abbrev-ref", "HEAD")
    if branch == RELEASE_BRANCH:
        ok(f"on branch '{RELEASE_BRANCH}'")
    elif ALLOW_BRANCH:
        ok(f"on branch '{branch}' (ALLOW_BRANCH=1 override)")
    else:
        block(
            f"must be on '{RELEASE_BRANCH}' (current: '{branch}'); "
            "set ALLOW_BRANCH=1 to override"
        )

    # 5. Working tree clean (no tracked changes).
    dirty = (
        git("diff", "--quiet", check=False).returncode != 0
        or git("diff", "--cached", "--quiet", check=False).returncode != 0
    )
    if dirty:
        block("working tree has uncommitted changes")
    else:
        ok("working tree clean (no tracked changes)")

    # 6. In sync with the remote release branch.
    git("fetch", "--quiet", REMOTE, RELEASE_BRANCH, check=False)
    if (
        git(
            "rev-parse", "-q", "--verify", f"{REMOTE}/{RELEASE_BRANCH}", check=False
        ).returncode
        != 0
    ):
        block(f"cannot find {REMOTE}/{RELEASE_BRANCH}")
    else:
        local_sha = git_out("rev-parse", "HEAD")
        remote_sha = git_out("rev-parse", f"{REMOTE}/{RELEASE_BRANCH}")
        if local_sha == remote_sha:
            ok(f"HEAD is in sync with {REMOTE}/{RELEASE_BRANCH}")
        else:
            base = git_out("merge-base", "HEAD", f"{REMOTE}/{RELEASE_BRANCH}")
            if local_sha == base:
                block(f"HEAD is behind {REMOTE}/{RELEASE_BRANCH} (pull first)")
            elif remote_sha == base:
                block(f"HEAD is ahead of {REMOTE}/{RELEASE_BRANCH} (unpushed commits)")
            else:
                block(f"HEAD has diverged from {REMOTE}/{RELEASE_BRANCH}")

    # Plan summary.
    print("\nplan:")
    if bump:
        print(f"  - set __version__ = {requested} in {INIT_FILE}")
        print(f"  - commit: 'Release {requested}'")
    print(f"  - create annotated tag: {requested}")
    push = (
        f"git push {REMOTE} {RELEASE_BRANCH} && git push {REMOTE} {requested}"
        if bump
        else f"git push {REMOTE} {requested}"
    )
    print(f"  - push (manual): {push}")
    print("  - pushing the tag triggers .github/workflows/cd.yml -> PyPI\n")

    if blockers:
        die(f"{len(blockers)} blocking issue(s); resolve them before releasing")

    if not CONFIRM:
        print("plan only (no changes made). Re-run with CONFIRM=1 to apply:")
        print(f"    make release {requested} CONFIRM=1")
        return 0

    # --- apply ---
    if not SKIP_GATE:
        print("running gate: make all ...")
        subprocess.run([os.environ.get("MAKE", "make"), "all"], check=True)

    if bump:
        text = INIT_FILE.read_text()
        match = VERSION_RE.search(text)
        assert match is not None
        INIT_FILE.write_text(text[: match.start(1)] + requested + text[match.end(1) :])
        if read_version() != requested:
            die(f"failed to update __version__ to {requested}")
        git("add", str(INIT_FILE))
        git("commit", "-m", f"Release {requested}")
        print(f"committed version bump to {requested}")

    git("tag", "-a", requested, "-m", f"Release {requested}")
    print(f"created annotated tag {requested}\n")
    print("Done. NOTHING has been pushed. To publish, run:")
    print(f"    {push}")
    print("(pushing the tag triggers the cd workflow which publishes to PyPI)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
