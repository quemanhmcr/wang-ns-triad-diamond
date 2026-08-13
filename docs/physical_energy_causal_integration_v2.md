# Physical-energy causal integration v2: exact-content reuse without weaker experiments

Status: certified engineering layer on exact implementation SHA `026b49337083433fbcf9764ab145ca5f7066f4a2`.

This layer changes **scheduling and provenance**, not Navier--Stokes/PDE mathematics.  It was introduced because the old physical-energy integration job reran every historical stress program sequentially whenever a new theorem was added, even when those programs and their entire executable dependency closures were unchanged.

## 1. What is never weakened

Every integration node retains its exact command line, sample count, seed, tolerance and output semantics.  The v2 manifest contains all `57` legacy module commands unchanged and adds the resolved-contact integration node as the first new node.  The full repository theorem suite and the exact `20,000`-trace master episode stress still run on every v2 invocation, including an invocation where all integration-node outputs are reusable.

A reusable node is therefore not a skipped theorem check.  It is a previously executed deterministic integration experiment whose executable identity and certified output identity are both proved unchanged, while the global test/master guards are rerun on the current SHA.

## 2. Content identity is fail-closed

For each node the planner hashes:

- the exact command, including sample count and seed;
- the pinned runtime identity: Ubuntu 24.04, Python 3.11.15, NumPy 2.2.6, SciPy 1.15.3, pytest 8.4.1 and python-flint 0.9.0;
- `requirements.txt`;
- the full transitive local `src.*` import closure and every byte of those source files.

AST-discovered dynamic import makes a node non-reusable.  The current integration closure has no runtime external-file reads and no dynamic-import node.  New nodes, changed commands, changed runtime/source closure, missing baseline metadata, or an explicit full-sweep request all execute rather than reuse.

The runtime pins are also regression-tested against the workflow and requirements file, so changing the runner contract without changing the fingerprint contract fails the full theorem suite.

## 3. Certified output identity is separate from source identity

Source identity alone is insufficient: a recorded artifact could be corrupted after certification.  Every baseline node therefore also stores a deterministic digest of its output directory.  Reuse first recomputes that digest on the stored certified tree and fails if it differs.  The aggregate stage then verifies the materialized/current node output against the same digest and performs a recursive byte comparison against the baseline tree.

Thus reuse requires both

`exact executable fingerprint` **and** `exact certified output digest`.

A test deliberately tampers with a stored baseline output and requires the materializer to reject it.

## 4. Five-way full-sweep sharding

The old exact-SHA integration `31655796159` took `30m54s`.  Profiling showed that its `57` module stress commands consumed about `1545` sequential step-seconds, while the full theorem suite was already a separate roughly five-minute block.  Longest-processing-time partitioning of the measured module runtimes into five shards gave predicted loads `309,310,310,308,308` seconds.  Six shards would save little wall time relative to the full-suite floor while spending an extra runner.

The first v2 run deliberately disabled reuse and executed **all 58 nodes** in those five shards.  Run `31658965931` completed successfully in `6m03s`.  Its global lane passed `912` tests and the master stress with worst margin `0.0`.  The aggregate stage found all 58 result directories and byte-compared all `57` baseline-matched legacy nodes successfully.  This is the differential certification that sharding changed only scheduling.

## 5. Exact-content incremental path

After the shadow-full result became the 58-node certified baseline, run `31659354320` planned

`58 reuse / 0 execute / 0 execution shards`.

It completed in `5m21s`; the full `912`-test suite and master stress still ran, and all `58/58` node outputs were byte-checked against the certified baseline.

The final integrity-hardened implementation, exact SHA `026b49337083433fbcf9764ab145ca5f7066f4a2`, ran as `31659688199`.  It completed in `3m28s`, passed `914` tests, retained master worst margin `0.0`, planned `58 reuse / 0 execute`, and verified all `58` outputs.  The additional tests cover baseline-output tampering and runtime-contract drift.

Relative to the measured `30m54s` legacy run, the shadow **full-compute** v2 route is about `5.1x` faster in wall time, and the certified exact-content route is about `8.9x` faster on this unchanged-node case.  No physical sample was removed to obtain either reduction.

## 6. How a new theorem enters

A theorem adds its own integration node to the manifest.  That new node has no baseline fingerprint, so it must execute.  Any historical node whose transitive executable closure changes also executes.  Unchanged historical nodes reuse their exact certified outputs.  The global theorem suite/master guard still runs.

This is deliberately analogous to the physical ontology of the project: a new event does not cause unrelated old physical events to happen again merely because the observer has advanced.  But the analogy is only engineering intuition; CI cache/reuse is **not** a Navier--Stokes causal law and creates no theorem owner, clock, stock or probability.

The workflow retains an explicit `force_full` mode.  Full sweeps remain available for periodic differential certification, baseline renewal, runtime/toolchain changes, or any situation in which provenance cannot prove safe reuse.
