from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

RUNTIME_ID = "ubuntu-24.04|python=3.11.15|numpy=2.2.6|scipy=1.15.3|pytest=8.4.1|python-flint=0.9.0"


@dataclass(frozen=True)
class NodePlan:
    node_id: str
    shard: int
    command: str
    outdir: str
    fingerprint: str
    closure: tuple[str, ...]
    dynamic_import: bool
    baseline_match: bool
    action: str
    reason: str


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def load_manifest(path: Path) -> dict:
    data = _read_json(path)
    nodes = data.get("nodes", [])
    shard_count = int(data.get("shard_count", 0))
    if shard_count <= 0 or not nodes:
        raise ValueError("nonempty integration manifest and positive shard count required")
    ids, outdirs, commands = set(), set(), set()
    for node in nodes:
        node_id = str(node["id"])
        module = str(node["module"])
        command = str(node["command"])
        outdir = str(node["outdir"])
        shard = int(node["shard"])
        if node_id in ids or outdir in outdirs or command in commands:
            raise ValueError("integration node ids, outdirs and commands must be unique")
        ids.add(node_id); outdirs.add(outdir); commands.add(command)
        argv = shlex.split(command)
        if argv[:3] != ["python", "-m", module]:
            raise ValueError(f"node {node_id} command/module mismatch")
        if "--outdir" not in argv or argv[argv.index("--outdir") + 1] != outdir:
            raise ValueError(f"node {node_id} outdir provenance mismatch")
        if not (1 <= shard <= shard_count):
            raise ValueError(f"node {node_id} has invalid shard")
        if not outdir.startswith("results-integration/"):
            raise ValueError(f"node {node_id} must stay under results-integration")
    return data


def legacy_commands(workflow: Path) -> tuple[str, ...]:
    commands = []
    for line in workflow.read_text().splitlines():
        match = re.search(r"- run: (python -m src\..+)$", line)
        if match:
            commands.append(match.group(1).strip())
    if not commands:
        raise ValueError("legacy integration workflow exposes no module commands")
    return tuple(commands)


def verify_legacy_coverage(manifest: dict, legacy_workflow: Path) -> dict[str, int]:
    legacy = legacy_commands(legacy_workflow)
    current = tuple(str(n["command"]) for n in manifest["nodes"])
    missing = sorted(set(legacy) - set(current))
    if missing:
        raise AssertionError(f"v2 manifest dropped/changed legacy commands: {missing}")
    return {
        "legacy_command_count": len(legacy),
        "manifest_command_count": len(current),
        "new_v2_command_count": len(set(current) - set(legacy)),
    }


def _module_path(repo: Path, module: str) -> Path:
    return repo / (module.replace(".", "/") + ".py")


def _local_imports(path: Path) -> tuple[set[str], bool]:
    tree = ast.parse(path.read_text(), filename=str(path))
    imports: set[str] = set()
    dynamic = False
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("src."):
                    imports.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.module.startswith("src."):
                imports.add(node.module)
        elif isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id == "__import__":
                dynamic = True
            elif isinstance(func, ast.Attribute) and func.attr == "import_module":
                dynamic = True
    return imports, dynamic


def source_closure(repo: Path, module: str) -> tuple[tuple[str, ...], bool]:
    todo = [module]
    seen: set[str] = set()
    files: list[str] = []
    dynamic = False
    while todo:
        current = todo.pop()
        if current in seen:
            continue
        seen.add(current)
        path = _module_path(repo, current)
        if not path.exists():
            raise FileNotFoundError(f"manifest module is missing: {current}")
        rel = path.relative_to(repo).as_posix()
        files.append(rel)
        imports, has_dynamic = _local_imports(path)
        dynamic = dynamic or has_dynamic
        for dep in sorted(imports - seen):
            if _module_path(repo, dep).exists():
                todo.append(dep)
    return tuple(sorted(files)), dynamic


def output_tree_digest(root: Path) -> str:
    """Deterministic content digest of one certified node-output directory.

    Paths are hashed relative to the node root, followed by file bytes.  The
    digest intentionally ignores mtimes, ownership and filesystem enumeration
    order so an artifact can move between GitHub storage and the repository
    without changing identity.
    """
    base = Path(root)
    if not base.is_dir():
        raise FileNotFoundError(f"certified node-output directory missing: {base}")
    files = sorted(p for p in base.rglob('*') if p.is_file())
    if not files:
        raise ValueError(f"certified node-output directory is empty: {base}")
    h = hashlib.sha256()
    for path in files:
        rel = path.relative_to(base).as_posix()
        h.update(rel.encode() + b'\0' + path.read_bytes() + b'\0')
    return h.hexdigest()


def node_fingerprint(repo: Path, node: dict, closure: Iterable[str]) -> str:
    h = hashlib.sha256()
    h.update(("runtime\0" + RUNTIME_ID + "\0").encode())
    h.update(("command\0" + str(node["command"]) + "\0").encode())
    requirements = repo / "requirements.txt"
    h.update(b"requirements.txt\0" + requirements.read_bytes() + b"\0")
    for rel in sorted(closure):
        h.update(rel.encode() + b"\0" + (repo / rel).read_bytes() + b"\0")
    return h.hexdigest()


def make_baseline(repo: Path, manifest: dict, *, theorem_sha: str, record_sha: str, run_id: int, artifact_root: str) -> dict:
    if subprocess.run(["git", "-C", str(repo), "diff", "--quiet", theorem_sha, record_sha, "--", "src", "requirements.txt"]).returncode != 0:
        raise AssertionError("baseline theorem/record SHAs do not have identical executable src/requirements trees")
    rows = {}
    baseline_root = repo / artifact_root
    for node in manifest["nodes"]:
        rel = _relative_output(str(node["outdir"]))
        if not (baseline_root / rel).is_dir():
            continue
        closure, dynamic = source_closure(repo, str(node["module"]))
        if dynamic:
            raise AssertionError(f"baseline node {node['id']} uses dynamic imports and is not reusable")
        rows[str(node["id"])] = {
            "fingerprint": node_fingerprint(repo, node, closure),
            "output_digest": output_tree_digest(baseline_root / rel),
            "command": str(node["command"]),
            "outdir": str(node["outdir"]),
        }
    return {
        "schema": 1,
        "runtime_id": RUNTIME_ID,
        "theorem_sha": theorem_sha,
        "record_sha": record_sha,
        "integration_run_id": int(run_id),
        "artifact_root": artifact_root,
        "nodes": rows,
    }


def make_plan(repo: Path, manifest: dict, baseline: dict, *, force_full: bool = False) -> dict:
    if baseline.get("runtime_id") != RUNTIME_ID:
        force_full = True
    baseline_nodes = dict(baseline.get("nodes", {}))
    plans: list[NodePlan] = []
    for node in manifest["nodes"]:
        closure, dynamic = source_closure(repo, str(node["module"]))
        fp = node_fingerprint(repo, node, closure)
        base = baseline_nodes.get(str(node["id"]))
        baseline_match = bool(
            not dynamic
            and base is not None
            and isinstance(base.get("output_digest"), str)
            and len(base.get("output_digest")) == 64
            and base.get("command") == node["command"]
            and base.get("outdir") == node["outdir"]
            and base.get("fingerprint") == fp
        )
        if force_full:
            action, reason = "execute", "forced_full_sweep"
        elif dynamic:
            action, reason = "execute", "dynamic_import_fails_closed"
        elif base is None:
            action, reason = "execute", "new_integration_node"
        elif not isinstance(base.get("output_digest"), str) or len(base.get("output_digest", "")) != 64:
            action, reason = "execute", "certified_output_digest_missing"
        elif not baseline_match:
            action, reason = "execute", "transitive_source_runtime_or_command_changed"
        else:
            action, reason = "reuse", "exact_content_fingerprint_match"
        plans.append(NodePlan(str(node["id"]), int(node["shard"]), str(node["command"]), str(node["outdir"]), fp, closure, dynamic, baseline_match, action, reason))
    execute_shards = sorted({p.shard for p in plans if p.action == "execute"})
    return {
        "schema": 1,
        "runtime_id": RUNTIME_ID,
        "baseline_record_sha": baseline.get("record_sha"),
        "baseline_integration_run_id": baseline.get("integration_run_id"),
        "force_full": force_full,
        "node_count": len(plans),
        "execute_count": sum(p.action == "execute" for p in plans),
        "reuse_count": sum(p.action == "reuse" for p in plans),
        "execution_shards": execute_shards,
        "nodes": [p.__dict__ for p in plans],
    }


def _relative_output(outdir: str) -> Path:
    p = Path(outdir)
    if len(p.parts) < 2 or p.parts[0] != "results-integration":
        raise ValueError("integration output escaped root")
    return Path(*p.parts[1:])


def materialize_reusable(repo: Path, plan: dict, baseline: dict, result_root: Path) -> None:
    baseline_root = repo / str(baseline["artifact_root"])
    for node in plan["nodes"]:
        if node["action"] != "reuse":
            continue
        rel = _relative_output(node["outdir"])
        src = baseline_root / rel
        dst = result_root / rel
        if not src.is_dir():
            raise FileNotFoundError(f"certified baseline result missing for {node['node_id']}: {src}")
        expected = baseline["nodes"][node["node_id"]].get("output_digest")
        actual = output_tree_digest(src)
        if actual != expected:
            raise AssertionError(
                f"certified baseline output digest mismatch for {node['node_id']}: {actual} != {expected}"
            )
        if dst.exists():
            shutil.rmtree(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst)


def execute_shard(repo: Path, plan: dict, shard: int) -> dict:
    rows = []
    for node in plan["nodes"]:
        if int(node["shard"]) != int(shard) or node["action"] != "execute":
            continue
        argv = shlex.split(str(node["command"]))
        subprocess.run(argv, cwd=repo, check=True)
        rows.append({"node_id": node["node_id"], "command": node["command"], "fingerprint": node["fingerprint"]})
    return {"schema": 1, "shard": int(shard), "executed": rows}


def verify_result_set(repo: Path, manifest: dict, plan: dict, baseline: dict, result_root: Path, *, compare_reused: bool) -> dict:
    baseline_root = repo / str(baseline["artifact_root"])
    missing = []
    compared = 0
    for node in plan["nodes"]:
        rel = _relative_output(node["outdir"])
        current = result_root / rel
        if not current.is_dir():
            missing.append(node["node_id"])
            continue
        if compare_reused and node.get("baseline_match"):
            base = baseline_root / rel
            expected = baseline["nodes"][node["node_id"]].get("output_digest")
            base_digest = output_tree_digest(base)
            current_digest = output_tree_digest(current)
            if base_digest != expected:
                raise AssertionError(
                    f"stored certified output for {node['node_id']} failed its baseline digest: {base_digest} != {expected}"
                )
            if current_digest != expected:
                raise AssertionError(
                    f"integration output for {node['node_id']} differs from certified digest: {current_digest} != {expected}"
                )
            proc = subprocess.run(["diff", "-qr", str(base), str(current)], text=True, capture_output=True)
            if proc.returncode != 0:
                raise AssertionError(f"baseline-matched node {node['node_id']} differs from certified baseline:\n{proc.stdout}{proc.stderr}")
            compared += 1
    if missing:
        raise AssertionError(f"integration result set missing nodes: {missing}")
    return {"schema": 1, "node_count": len(plan["nodes"]), "baseline_matched_nodes_byte_compared": compared, "all_results_present": True}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("verify", "baseline", "plan", "materialize", "execute", "verify-results"))
    ap.add_argument("--repo", type=Path, default=Path("."))
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--legacy-workflow", type=Path)
    ap.add_argument("--baseline", type=Path)
    ap.add_argument("--plan", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--theorem-sha")
    ap.add_argument("--record-sha")
    ap.add_argument("--run-id", type=int)
    ap.add_argument("--artifact-root")
    ap.add_argument("--force-full", action="store_true")
    ap.add_argument("--shard", type=int)
    ap.add_argument("--result-root", type=Path, default=Path("results-integration"))
    args = ap.parse_args()
    repo = args.repo.resolve(); manifest = load_manifest(args.manifest)
    if args.command == "verify":
        if not args.legacy_workflow: raise SystemExit("--legacy-workflow required")
        payload = verify_legacy_coverage(manifest, args.legacy_workflow)
    elif args.command == "baseline":
        if not all((args.theorem_sha, args.record_sha, args.run_id, args.artifact_root, args.out)): raise SystemExit("baseline args missing")
        payload = make_baseline(repo, manifest, theorem_sha=args.theorem_sha, record_sha=args.record_sha, run_id=args.run_id, artifact_root=args.artifact_root)
    else:
        if not args.baseline: raise SystemExit("--baseline required")
        baseline = _read_json(args.baseline)
        if args.command == "plan":
            payload = make_plan(repo, manifest, baseline, force_full=args.force_full)
        else:
            if not args.plan: raise SystemExit("--plan required")
            plan = _read_json(args.plan)
            if args.command == "materialize":
                materialize_reusable(repo, plan, baseline, args.result_root); payload={"materialized": plan["reuse_count"]}
            elif args.command == "execute":
                if args.shard is None: raise SystemExit("--shard required")
                payload = execute_shard(repo, plan, args.shard)
            else:
                payload = verify_result_set(repo, manifest, plan, baseline, args.result_root, compare_reused=True)
    text=json.dumps(payload, indent=2, sort_keys=True)+"\n"
    if args.out: args.out.write_text(text)
    else: print(text, end="")

if __name__ == "__main__":
    main()
