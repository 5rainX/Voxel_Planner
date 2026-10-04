"""Deterministic split by environment content hash."""

import random
from collections import defaultdict


def split_samples(samples, seed=42):
    groups = defaultdict(list)
    for sample in samples:
        groups[sample["provenance"]["map_id"]].append(sample["sample_id"])
    keys = sorted(groups)
    random.Random(seed).shuffle(keys)
    if len(keys) >= 3:
        n_val = max(1, round(len(keys) * 0.2))
        n_test = max(1, round(len(keys) * 0.2))
        if n_val + n_test >= len(keys):
            n_val = n_test = 1
        assignments = dict.fromkeys(keys[:n_val], "val")
        assignments.update(dict.fromkeys(keys[n_val:n_val + n_test], "test"))
        assignments.update(dict.fromkeys(keys[n_val + n_test:], "train"))
    else:
        assignments = {key: "train" for key in keys}
    result = {name: [] for name in ("train", "val", "test")}
    for key in sorted(groups):
        result[assignments[key]].extend(sorted(groups[key]))
    return result, "PROVISIONAL_SMALL_DATASET" if len(keys) < 30 else "GROUPED"
