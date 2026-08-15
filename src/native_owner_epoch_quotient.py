from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence

from src.material_native_joint_stop_projection import MaterialNativeJointStopProjection


STATUS = (
    "DRAFT_NATIVE_OWNER_EPOCH_QUOTIENT__"
    "MATERIAL_ONLY_OBSERVATIONS_DO_NOT_BREAK_PHYSICAL_EPOCHS__"
    "EXISTING_PURE_EPOCH_TELESCOPES_APPLY_TO_EVENT_SUBSEQUENCE"
)


class NativeEpochKind(str, Enum):
    HIGH_STRAIN = "high_strain"
    SIGNED_GOOD_GENERATED_HH = "signed_good_generated_hh"
    OTHER_NATIVE_EVENT = "other_native_event"
    MATERIAL_ONLY_OBSERVATION = "material_only_observation"


@dataclass(frozen=True)
class NativeEpochRecord:
    time: float
    kind: NativeEpochKind
    projection: MaterialNativeJointStopProjection | None = None
    witness: str = ""

    def __post_init__(self) -> None:
        if self.time < 0.0:
            raise ValueError("nonnegative physical time required")
        if self.kind is NativeEpochKind.MATERIAL_ONLY_OBSERVATION:
            if self.projection is None or not self.projection.no_causal_stop:
                raise ValueError("material-only record must be certified as no causal stop")
        elif self.projection is not None and self.projection.no_causal_stop:
            raise ValueError("a native event record cannot carry a no-causal-stop projection")


@dataclass(frozen=True)
class NativeOwnerEpochQuotient:
    event_indices: tuple[int, ...]
    event_kinds: tuple[str, ...]
    event_times: tuple[float, ...]
    high_strain_runs: tuple[tuple[int, int], ...]
    signed_good_hh_runs: tuple[tuple[int, int], ...]
    material_observations_removed: int
    material_epoch_breakers_created: int = 0
    event_order_changed: bool = False

    def __post_init__(self) -> None:
        if len(self.event_indices) != len(self.event_kinds) or len(self.event_indices) != len(self.event_times):
            raise ValueError("matching quotient event coordinates required")
        if any(b <= a for a, b in zip(self.event_indices, self.event_indices[1:])):
            raise ValueError("event subsequence indices must be strictly increasing")
        if any(t1 > t0 for t0, t1 in zip(self.event_times, self.event_times[1:])):
            raise ValueError("backward event path times must be nonincreasing")
        if self.material_observations_removed < 0:
            raise ValueError("nonnegative removed-observation count required")
        if self.material_epoch_breakers_created or self.event_order_changed:
            raise ValueError("material quotient cannot create epoch breakers or reorder physical events")
        for runs in (self.high_strain_runs, self.signed_good_hh_runs):
            for start, stop in runs:
                if not (0 <= start < stop <= len(self.event_kinds)):
                    raise ValueError("invalid maximal event-run coordinates")


def _maximal_runs(kinds: Sequence[str], target: str) -> tuple[tuple[int, int], ...]:
    out: list[tuple[int, int]] = []
    i = 0
    while i < len(kinds):
        if kinds[i] != target:
            i += 1
            continue
        j = i + 1
        while j < len(kinds) and kinds[j] == target:
            j += 1
        out.append((i, j))
        i = j
    return tuple(out)


def quotient_material_observations_from_epoch_path(
    records: Sequence[NativeEpochRecord],
) -> NativeOwnerEpochQuotient:
    """Delete only certified no-causal-stop material observations from a path.

    The operation is order-preserving and keeps every actual physical event at
    its original physical time. Therefore a material-only observation between
    two high-strain events cannot split their high-strain epoch; likewise for
    signed-good generated-HH events. Any independently witnessed different
    native event remains in the subsequence and is a genuine epoch breaker.
    """
    if not records:
        raise ValueError("nonempty backward observation path required")
    times = tuple(float(r.time) for r in records)
    if any(t1 > t0 for t0, t1 in zip(times, times[1:])):
        raise ValueError("records must be ordered backward in nonincreasing physical time")

    event_indices = tuple(i for i, r in enumerate(records) if r.kind is not NativeEpochKind.MATERIAL_ONLY_OBSERVATION)
    event_records = tuple(records[i] for i in event_indices)
    kinds = tuple(r.kind.value for r in event_records)
    event_times = tuple(float(r.time) for r in event_records)
    removed = len(records) - len(event_records)

    return NativeOwnerEpochQuotient(
        event_indices=event_indices,
        event_kinds=kinds,
        event_times=event_times,
        high_strain_runs=_maximal_runs(kinds, NativeEpochKind.HIGH_STRAIN.value),
        signed_good_hh_runs=_maximal_runs(kinds, NativeEpochKind.SIGNED_GOOD_GENERATED_HH.value),
        material_observations_removed=removed,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "event_subsequence": "remove only material observations already certified no_causal_stop; keep every genuine physical/internal event at the same physical time and in the same order",
        "high_strain": "material-only observations cannot split a maximal high-strain event run; the existing certified descending high-strain epoch telescope applies to the resulting event subsequence",
        "signed_good_hh": "material-only observations cannot split a maximal signed-good generated-HH event run; the existing certified physical-time telescope applies to the resulting event subsequence",
        "genuine_breakers": "any independently witnessed different native event remains in the subsequence and still breaks the pure epoch",
        "scope": "topological quotient only; it does not prove that the surviving mixed native-event word terminates and makes no Navier-Stokes global-regularity claim",
    }
