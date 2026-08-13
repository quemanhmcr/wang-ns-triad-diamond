from __future__ import annotations

import json
from pathlib import Path

import pytest

from tools.physical_integration_plan import (
    RUNTIME_ID,
    load_manifest,
    make_plan,
    materialize_reusable,
    node_fingerprint,
    source_closure,
    verify_legacy_coverage,
)


def _write_repo(root: Path) -> None:
    (root / 'src').mkdir()
    (root / 'src/a.py').write_text('from src.b import VALUE\n')
    (root / 'src/b.py').write_text('VALUE = 1\n')
    (root / 'requirements.txt').write_text('numpy==2.2.6\n')


def _node() -> dict:
    return {
        'id': 'a',
        'module': 'src.a',
        'command': 'python -m src.a --outdir results-integration/a',
        'outdir': 'results-integration/a',
        'shard': 1,
    }


def test_transitive_source_change_invalidates_fingerprint(tmp_path: Path):
    _write_repo(tmp_path)
    node = _node()
    closure, dynamic = source_closure(tmp_path, 'src.a')
    assert dynamic is False
    assert closure == ('src/a.py', 'src/b.py')
    first = node_fingerprint(tmp_path, node, closure)
    (tmp_path / 'src/b.py').write_text('VALUE = 2\n')
    second = node_fingerprint(tmp_path, node, closure)
    assert first != second


def test_dynamic_import_fails_closed(tmp_path: Path):
    _write_repo(tmp_path)
    (tmp_path / 'src/a.py').write_text("import importlib\nimportlib.import_module('src.b')\n")
    _, dynamic = source_closure(tmp_path, 'src.a')
    assert dynamic is True


def test_exact_fingerprint_reuses_but_new_node_executes(tmp_path: Path):
    _write_repo(tmp_path)
    node = _node()
    manifest = {'shard_count': 1, 'nodes': [node]}
    closure, _ = source_closure(tmp_path, 'src.a')
    fp = node_fingerprint(tmp_path, node, closure)
    baseline = {'runtime_id': RUNTIME_ID, 'record_sha': 'base', 'integration_run_id': 1, 'artifact_root': 'baseline', 'nodes': {'a': {'fingerprint': fp, 'command': node['command'], 'outdir': node['outdir']}}}
    plan = make_plan(tmp_path, manifest, baseline)
    assert plan['reuse_count'] == 1 and plan['execute_count'] == 0
    manifest2 = {'shard_count': 1, 'nodes': [node, {'id':'new','module':'src.b','command':'python -m src.b --outdir results-integration/new','outdir':'results-integration/new','shard':1}]}
    plan2 = make_plan(tmp_path, manifest2, baseline)
    assert plan2['reuse_count'] == 1 and plan2['execute_count'] == 1
    assert next(n for n in plan2['nodes'] if n['node_id']=='new')['reason'] == 'new_integration_node'


def test_materialize_reusable_copies_only_certified_directory(tmp_path: Path):
    _write_repo(tmp_path)
    baseline_root = tmp_path / 'baseline/a'; baseline_root.mkdir(parents=True); (baseline_root/'summary.md').write_text('certified\n')
    plan = {'reuse_count':1,'nodes':[{'node_id':'a','action':'reuse','outdir':'results-integration/a'}]}
    baseline = {'artifact_root':'baseline'}
    result = tmp_path / 'results-integration'
    materialize_reusable(tmp_path, plan, baseline, result)
    assert (result/'a/summary.md').read_text() == 'certified\n'


def test_repository_manifest_never_drops_a_legacy_command():
    repo = Path(__file__).resolve().parents[1]
    manifest = load_manifest(repo / 'ci/physical_energy_integration_manifest.json')
    out = verify_legacy_coverage(manifest, repo / '.github/workflows/physical-energy-causal-integration.yml')
    assert out['legacy_command_count'] == 57
    assert out['manifest_command_count'] >= out['legacy_command_count']
