"""Write the A06 → A06.1 comparison from real derived artifacts."""

import argparse
import json
from pathlib import Path


def correction_counts(root):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    train = set((root / "splits" / "train.txt").read_text(encoding="utf-8").splitlines())
    total = in_train = 0
    for record in manifest["samples"]:
        sample = json.loads((root / record["path"]).read_text(encoding="utf-8"))
        nonzero = any(abs(target["du"]) > 1e-4 or abs(target["dv"]) > 1e-4
                      for target in sample["targets"])
        total += nonzero
        in_train += nonzero and sample["sample_id"] in train
    return total, in_train


def compare(a06, a061):
    a06, a061 = Path(a06), Path(a061)
    before = json.loads((a06 / "stats.json").read_text(encoding="utf-8"))
    after = json.loads((a061 / "dataset" / "stats.json").read_text(encoding="utf-8"))
    reference = json.loads((a061 / "reference_geometry" / "stats.json").read_text(encoding="utf-8"))
    before_nonzero, before_train_nonzero = correction_counts(a06)
    after_nonzero, after_train_nonzero = correction_counts(a061 / "dataset")
    rows = [
        ("Source samples", before["total_source_samples"], after["total_source_samples"]),
        ("Geometry-valid", before["geometry_aligned_samples"], after["geometry_aligned_samples"]),
        ("Rejected", before["dataset_samples_rejected"], after["dataset_samples_rejected"]),
        ("Bend samples recovered", 0, reference["bend_recovered"]),
        ("du/dv-valid", before["du_dv_valid_samples"], after["du_dv_valid_samples"]),
        ("theta-valid samples", before["theta_valid_samples"], after["theta_valid_samples"]),
        ("theta-valid nodes", before["theta_valid_nodes"], after["theta_valid_nodes"]),
        ("SDF-valid", before["sdf_valid_samples"], after["sdf_valid_samples"]),
        ("Full schema-valid", before["fully_valid_training_samples"], after["fully_valid_training_samples"]),
        ("Nonzero correction samples", before_nonzero, after_nonzero),
        ("Train nonzero correction samples", before_train_nonzero, after_train_nonzero)]
    report = ["# A06 → A06.1 comparison", "", "| Metric | A06 | A06.1 |", "|---|---:|---:|"]
    report += [f"| {name} | {old} | {new} |" for name, old, new in rows]
    report += ["", f"Reference records: {reference['enriched_valid']} valid, {reference['unrecoverable']} unrecoverable.",
               f"Dense points exported: {reference['dense_points']}; orientation-valid points: {reference['orientation_valid_nodes']}.",
               f"Nonzero correction nodes: {after['nonzero_correction_nodes']}; nonzero theta nodes: {after['theta_nonzero_nodes']}.",
               f"A06.1 reconstruction error mean/max/P95 (voxel): {after['reconstruction_error_voxel']['mean']:.6g} / {after['reconstruction_error_voxel']['max']:.6g} / {after['reconstruction_error_voxel']['p95']:.6g}.",
               f"Split: {after['split_label_coverage']}; status PROVISIONAL_SMALL_DATASET.",
               "", "The external fine-planner bend curve and twist interior orientation were not available in the historical export.",
               "A07 data readiness: NO-GO; training lacks nonzero du/dv examples and bend coverage.", ""]
    output = a061 / "report.md"
    with output.open("x", encoding="utf-8") as stream:
        stream.write("\n".join(report))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a06", type=Path, required=True)
    parser.add_argument("--a06-1", type=Path, required=True)
    args = parser.parse_args()
    print(compare(args.a06, args.a06_1))


if __name__ == "__main__":
    main()
