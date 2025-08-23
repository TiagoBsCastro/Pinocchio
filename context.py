#!/usr/bin/env python3
import argparse, os, sys, subprocess, textwrap, unicodedata
from pathlib import Path

DEFAULT_EXCLUDE_DIRS = {
    ".git","node_modules","dist","build",".venv","venv","__pycache__",
    ".tox",".mypy_cache",".pytest_cache",".idea",".vscode",".DS_Store"
}
DEFAULT_EXCLUDE_EXT = {
    ".png",".jpg",".jpeg",".gif",".webp",".svg",".pdf",".zip",".tar",".gz",".bz2",".7z",
    ".mp3",".wav",".flac",".mp4",".mov",".avi",".mkv",
    ".exe",".dll",".so",".dylib",".bin",".class",".o",".a",".jar",
    ".otf",".ttf",".woff",".woff2",".ico",".lock"
}
DEFAULT_SECRET_BASENAMES = {
    ".env","id_rsa","id_ed25519","credentials.json","service-account.json","secret.json",
    "secrets.json","firebase-service-account.json","google-credentials.json"
}
DEFAULT_PRIORITY_FILES = [
    "README.md","README","README.rst","README.txt",
    "pyproject.toml","requirements.txt","environment.yml","Pipfile","Pipfile.lock",
    "package.json","yarn.lock","pnpm-lock.yaml","go.mod","go.sum",
    "Cargo.toml","Cargo.lock","Gemfile","pom.xml","Makefile","Dockerfile",".env.example"
]

def is_binary(path: Path, sample_bytes=8192) -> bool:
    try:
        with path.open("rb") as f:
            chunk = f.read(sample_bytes)
        if b"\x00" in chunk:
            return True
        # Heuristic: too many non-text characters
        text_ratio = sum(32 <= b <= 126 or b in (9,10,13) for b in chunk) / max(1, len(chunk))
        return text_ratio < 0.80
    except Exception:
        return True

def git_ls_files(repo: Path):
    try:
        out = subprocess.check_output(["git","ls-files"], cwd=repo, text=True)
        return [repo / p for p in out.splitlines()]
    except Exception:
        # Fallback to walk
        return [p for p in repo.rglob("*") if p.is_file()]

def within_size_limits(path: Path, max_bytes_per_file: int) -> bool:
    try:
        return path.stat().st_size <= max_bytes_per_file
    except Exception:
        return False

def looks_texty(path: Path) -> bool:
    # Allow files without extensions if not binary
    ext = path.suffix.lower()
    if ext in DEFAULT_EXCLUDE_EXT:
        return False
    return True

def normalize(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    # Strip trailing spaces to save tokens
    return "\n".join(line.rstrip() for line in s.splitlines())

def read_head(path: Path, max_lines: int) -> str:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            lines = []
            for i, line in enumerate(f, 1):
                lines.append(line)
                if i >= max_lines:
                    break
        return normalize("".join(lines))
    except Exception as e:
        return f"[Could not read file: {e}]"

def main():
    ap = argparse.ArgumentParser(description="Build a compact context.md from a repo.")
    ap.add_argument("--repo", default=".", help="Path to repo root")
    ap.add_argument("--out", default="context.md", help="Output file")
    ap.add_argument("--max-lines-per-file", type=int, default=400)
    ap.add_argument("--max-bytes-per-file", type=int, default=200_000)
    ap.add_argument("--max-total-bytes", type=int, default=5_000_000)
    ap.add_argument("--depth", type=int, default=3, help="Tree depth")
    args = ap.parse_args()

    root = Path(args.repo).resolve()
    files = git_ls_files(root)

    # Filter
    candidates = []
    for p in files:
        rel = p.relative_to(root)
        parts = set(rel.parts)
        if any(d in parts for d in DEFAULT_EXCLUDE_DIRS):
            continue
        if p.name in DEFAULT_SECRET_BASENAMES:
            continue
        if not looks_texty(p):
            continue
        if is_binary(p):
            continue
        if not within_size_limits(p, args.max_bytes_per_file):
            continue
        candidates.append(p)

    # Put priority files first
    priority = []
    others = []
    priority_names = set(DEFAULT_PRIORITY_FILES)
    for p in candidates:
        if p.name in priority_names or str(p.relative_to(root)) in DEFAULT_PRIORITY_FILES:
            priority.append(p)
        else:
            others.append(p)
    ordered = priority + sorted(others, key=lambda x: str(x).lower())

    # Build header
    def make_tree(root: Path, depth: int) -> str:
        lines = []
        def walk(base: Path, level: int):
            if level > depth: 
                return
            entries = sorted([p for p in base.iterdir() if p.name not in DEFAULT_EXCLUDE_DIRS], key=lambda x: (x.is_file(), x.name.lower()))
            for e in entries:
                indent = "  " * (level - 1)
                lines.append(f"{indent}- {e.name}{'/' if e.is_dir() else ''}")
                if e.is_dir():
                    walk(e, level + 1)
        lines.append(f"{root.name}/")
        walk(root, 1)
        return "\n".join(lines)

    header = [
        "# Repository context for LLM",
        "",
        "## How this was built",
        "- Files are text only and truncated to a head limit.",
        "- Common binary and build artifacts are excluded.",
        "- Likely secret files are excluded.",
        "",
        "## Directory tree (truncated)",
        "```",
        make_tree(root, args.depth),
        "```",
        ""
    ]

    total_bytes = sum(len(s)+1 for s in header)
    chunks = ["\n".join(header)]

    for p in ordered:
        if ("/tests/" in str(p)) or ("/CAMBFiles/" in str(p)) or (".example." in str(p)) or ("HMF_Validation" in str(p)) or ("example/log" in str(p)):
            continue
        rel = str(p.relative_to(root))
        body = read_head(p, args.max_lines_per_file)
        block = f"\n\n----- FILE: {rel} -----\n```text\n{body}\n```"
        b = len(block.encode("utf-8"))
        if total_bytes + b > args.max_total_bytes:
            chunks.append("\n\n[Truncated due to max_total_bytes limit]\n")
            break
        chunks.append(block)
        total_bytes += b

    out = "\n".join(chunks)
    Path(args.out).write_text(out, encoding="utf-8")
    print(f"Wrote {args.out} ({total_bytes} bytes)")

if __name__ == "__main__":
    main()