"""Run collaborator -> Voxel_Planner -> collaborator fine planning in isolation.

FINE_ACCEPTED means accepted by the supplied fine planner, not independently
certified manufacturing ground truth. Only complete successful pairs enter
accepted.jsonl. All attempts, including timeouts and errors, enter results.jsonl.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def read_text(path):
    return path.read_text(encoding="utf-8", errors="replace")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def key_values(path):
    if not path.is_file():
        return {}
    return dict(line.split("=", 1) for line in read_text(path).splitlines() if "=" in line)


def validate_map(path):
    # Both supplied callers share this sparse TXT representation. The fine
    # caller cannot load VTK, and silently ignores several malformed records.
    with path.open(encoding="utf-8-sig") as stream:
        header = stream.readline().split()
        if len(header) != 4:
            raise ValueError("Map header must be 'map: nx ny nz'")
        shape = tuple(int(v) for v in header[1:])
        if any(v <= 0 or v > 10000 for v in shape):
            raise ValueError("Map dimensions must be in [1, 10000]")
        count = 0
        for number, line in enumerate(stream, 2):
            if not line.strip():
                continue
            tokens = line.split()
            if len(tokens) not in (3, 4):
                raise ValueError(f"Invalid map record on line {number}")
            point = tuple(int(v) for v in tokens[:3])
            if any(v < 0 or v >= size for v, size in zip(point, shape)):
                raise ValueError(f"Out-of-bounds voxel on line {number}")
            if len(tokens) == 4 and tokens[3] not in ("0", "1"):
                raise ValueError(f"Unsupported voxel flag on line {number}")
            count += 1
    return {"shape": shape, "listed_voxels": count}


def numeric_rows(path, columns):
    rows = []
    for line in read_text(path).splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row = [float(value) for value in line.split()]
        if len(row) != columns or not all(math.isfinite(v) for v in row):
            raise ValueError("Invalid exported numeric data: " + path.name)
        rows.append(row)
    return rows


def read_segments(path):
    segments = []
    vectors = {"Start", "End", "Start Normal", "End Normal", "Start Tangent", "End Tangent"}
    scalars = {"Angle", "Length", "Cost"}
    for line in read_text(path).splitlines():
        match = re.fullmatch(r"Segment (\d+):", line.strip())
        if match:
            if int(match[1]) != len(segments):
                raise ValueError("Nonsequential segment indices")
            segments.append({})
        elif segments and not line.lstrip().startswith("#") and ":" in line:
            key, value = (part.strip() for part in line.split(":", 1))
            if key == "Type":
                segments[-1]["type"] = value
            elif key in vectors:
                segments[-1][key.lower().replace(" ", "_")] = [float(v.strip()) for v in value.strip("()").split(",")]
            elif key in scalars:
                segments[-1][key.lower()] = float(value)
    required = {key.lower().replace(" ", "_") for key in vectors | scalars} | {"type"}
    for segment in segments:
        if not required.issubset(segment):
            raise ValueError("Incomplete fine segment export")
        for key, value in segment.items():
            if key == "type":
                continue
            values = value if isinstance(value, list) else [value]
            if not all(math.isfinite(v) for v in values) or (isinstance(value, list) and len(value) != 3):
                raise ValueError("Invalid fine geometry export")
    if not segments:
        raise ValueError("Empty fine path export")
    return segments


def collect_result(directory, returncode, timed_out=False):
    status_files = sorted(directory.glob("voxel_planner_release*_status.txt"))
    coarse = {p.name: key_values(p) for p in status_files}
    errors = [v.get("message", "") for v in coarse.values() if v.get("status") == "ERROR"]
    stdout = read_text(directory / "stdout.log")
    errors += [line for line in stdout.splitlines() if "Voxel_Planner_Release error:" in line]
    table = directory / "planning_candidates.csv"
    rows = []
    if table.exists():
        with table.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
    result = {"coarse_status": coarse, "coarse_errors": errors, "candidates": rows,
              "returncode": returncode, "accepted_count": 0}
    if timed_out:
        result["status"] = "PROCESS_TIMEOUT"
        return result, []
    # Original main can overwrite ERROR with NO_PATH; preserve request errors.
    if not rows and errors:
        result["status"] = "COARSE_ERROR"
        return result, []
    if returncode not in (0, 1):
        result["status"] = "PROCESS_ERROR"
        return result, []
    if not rows:
        root_status = coarse.get("voxel_planner_release_status.txt", {}).get("status")
        result["status"] = "COARSE_NO_PATH" if returncode == 1 and root_status == "NO_PATH" else "OUTPUT_ERROR"
        return result, []
    successful = [row for row in rows if row.get("fine_success") == "true"]
    if not successful:
        result["status"] = "FINE_TIMEOUT" if any("PLANNER_TIMEOUT" in row.get("message", "") for row in rows) else "FINE_NO_PATH"
        if returncode == 0:
            result["status"] = "OUTPUT_ERROR"
        return result, []
    if returncode != 0:
        result["status"] = "OUTPUT_ERROR"
        return result, []
    pairs = []
    for row in successful:
        candidate = int(row["candidate_id"])
        prefix = "voxel_planner_release_coarse_group_{:03d}_route_{:03d}".format(
            int(row["coarse_tangent_group"]), int(row["route_index"]))
        status = key_values(directory / (prefix + "_status.txt"))
        fine_status = key_values(directory / f"busbar_path_candidate_{candidate:03d}_status.txt")
        if status.get("status") != "OK" or fine_status.get("success") != "true":
            raise ValueError("Candidate status files do not confirm both planning stages")
        centerline = numeric_rows(directory / (prefix + "_guide.txt"), 4)
        poses = numeric_rows(directory / (prefix + "_centerline_poses.txt"), 7)
        segments = read_segments(directory / f"busbar_path_candidate_{candidate:03d}.txt")
        if not centerline or len(centerline) != int(row["coarse_voxels"]) or len(segments) != int(row["segments"]):
            raise ValueError("Exported path counts do not match candidate summary")
        if len(poses) != len(centerline):
            raise ValueError("Missing per-point coarse poses")
        if any(point[0] != i or pose[0] != i for i, (point, pose) in enumerate(zip(centerline, poses))):
            raise ValueError("Coarse point/pose indices do not match")
        pairs.append({"candidate_id": candidate, "metrics": row,
                      "coarse_centerline": [p[1:] for p in centerline],
                      "coarse_poses": [p[1:] for p in poses], "fine_segments": segments,
                      "validation_scope": "collaborator_fine_planner",
                      "independent_geometry_validation": False,
                      "global_optimality_proven": False})
    result["status"] = "FINE_ACCEPTED"
    result["accepted_count"] = len(pairs)
    return result, pairs


def make_command(executable, case, fine_timeout, coarse_limit, export_obj):
    arguments = case.get("args", [])
    if not isinstance(arguments, list) or not all(isinstance(v, str) for v in arguments):
        raise ValueError("case.args must be an array of CLI argument strings")
    forbidden = {"--coarse-only", "--no-release-initial-path", "--fine-on-processed-map",
                 "--fine-on-processed-map-with-guide", "--astar-jump", "--legacy-coarse", "--nc",
                 "--start-tangent-rotate-deg", "--end-tangent-rotate-deg", "--goal-tangent-rotate-deg"}
    if forbidden.intersection(arguments):
        raise ValueError("Options would change or bypass the required coarse-to-fine workflow")
    for flag in ("--width", "--thickness", "--start-normal", "--end-normal", "--start-tangent", "--end-tangent"):
        if flag not in arguments:
            raise ValueError("Specify " + flag + " explicitly in case.args")
    endpoints = case["start"] + case["goal"]
    if len(case["start"]) != 3 or len(case["goal"]) != 3 or any(type(v) is not int for v in endpoints):
        raise ValueError("start and goal must each contain three integers")
    command = [str(executable), "plan", "input.txt"] + [str(v) for v in endpoints] + arguments
    if "--voxel-planner-release" not in arguments:
        command.append("--voxel-planner-release")
    # Explicit portfolio prevents the collaborator's legacy fallback from
    # silently succeeding without any Voxel_Planner guide.
    if "--endpoint-tangent-candidate" not in arguments:
        command += ["--endpoint-tangent-candidate", "0", "0"]
    if "--max-time" not in arguments:
        command += ["--max-time", str(fine_timeout)]
    if "--coarse-max-expansions-per-stage" not in arguments:
        command += ["--coarse-max-expansions-per-stage", str(coarse_limit)]
    if export_obj:
        command += [flag for flag in ("--coarse-obj", "--path-obj") if flag not in arguments]
    return command


def run_case(executable, case, manifest_dir, output, options):
    directory = output / case["id"]
    directory.mkdir()
    write_json(directory / "request.json", case)
    started = time.monotonic()
    pairs = []
    error_stage = "INPUT_ERROR"
    try:
        source = (manifest_dir / case["map"]).resolve()
        shutil.copyfile(source, directory / "input.txt")
        # Reject a UTF-8 BOM: the original native loaders treat it differently.
        with (directory / "input.txt").open("rb") as stream:
            bom = stream.read(3) == b"\xef\xbb\xbf"
        if bom:
            raise ValueError("Use a TXT map encoded without a UTF-8 BOM")
        metadata = validate_map(directory / "input.txt")
        command = make_command(executable, case, options.fine_timeout, options.coarse_limit, options.export_obj)
        for point in (case["start"], case["goal"]):
            if any(v < 0 or v >= n for v, n in zip(point, metadata["shape"])):
                raise ValueError("Endpoint lies outside map bounds")
        write_json(directory / "command.json", command)
        env = dict(os.environ, OMP_NUM_THREADS=str(options.threads))
        timed_out = False
        error_stage = "PROCESS_ERROR"
        with (directory / "stdout.log").open("wb") as out, (directory / "stderr.log").open("wb") as err:
            try:
                process = subprocess.run(command, cwd=directory, stdout=out, stderr=err,
                                         env=env, timeout=options.timeout, check=False)
                returncode = process.returncode
            except subprocess.TimeoutExpired:
                timed_out, returncode = True, None
        error_stage = "OUTPUT_ERROR"
        result, pairs = collect_result(directory, returncode, timed_out)
        result["map"] = metadata
        result["map_sha256"] = file_hash(directory / "input.txt")
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        result = {"status": error_stage,
                  "error": str(error), "accepted_count": 0}
        pairs = []
    result.update(case_id=case["id"], elapsed_s=time.monotonic() - started, request=case)
    expected = case.get("expected")
    result["expected"] = expected
    if expected == "NO_SOLUTION":
        result["expectation_met"] = result["status"] in {"COARSE_NO_PATH", "FINE_NO_PATH"}
    else:
        result["expectation_met"] = result["status"] == expected if expected else None
    write_json(directory / "result.json", result)
    for pair in pairs:
        pair.update(case_id=case["id"], sample_id=f"{case['id']}:{pair['candidate_id']}",
                    request=case, map_sha256=result["map_sha256"], map_file=f"{case['id']}/input.txt",
                    artifacts=case["id"])
    return result, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=60, help="Whole-process wall-clock limit, seconds")
    parser.add_argument("--fine-timeout", type=float, default=10, help="Limit per fine-planning attempt, seconds")
    parser.add_argument("--coarse-limit", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--export-obj", action="store_true")
    options = parser.parse_args()
    sys.stdout.reconfigure(errors="replace")
    executable, manifest_path, output = options.exe.resolve(), options.manifest.resolve(), options.output.resolve()
    if not executable.is_file() or not manifest_path.is_file():
        parser.error("Executable and manifest must exist")
    if any(not math.isfinite(v) or v <= 0 for v in (options.timeout, options.fine_timeout, options.coarse_limit, options.threads)):
        parser.error("Limits and thread count must be finite and positive")
    if output.exists():
        parser.error("Use a new --output directory; old results are never reused")
    manifest = json.loads(read_text(manifest_path))
    cases = manifest.get("cases", [])
    ids = [case.get("id", "") for case in cases]
    if not cases or len(set(ids)) != len(ids) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in ids):
        parser.error("Cases need unique nonempty ASCII ids: letters, digits, underscore or dash")
    output.mkdir(parents=True)
    write_json(output / "manifest.json", manifest)
    build_provenance = executable.parent.parent / "provenance.json"
    write_json(output / "run.json", {"executable": str(executable), "executable_sha256": file_hash(executable),
               "manifest": str(manifest_path), "limits": vars(options) | {"exe": str(executable), "manifest": str(manifest_path), "output": str(output)},
               "build": json.loads(read_text(build_provenance)) if build_provenance.exists() else None})
    results = []
    with (output / "results.jsonl").open("w", encoding="utf-8") as attempts, (output / "accepted.jsonl").open("w", encoding="utf-8") as accepted:
        for case in cases:
            result, pairs = run_case(executable, case, manifest_path.parent, output, options)
            results.append(result)
            attempts.write(json.dumps(result, ensure_ascii=False, allow_nan=False) + "\n")
            attempts.flush()
            for pair in pairs:
                accepted.write(json.dumps(pair, ensure_ascii=False, allow_nan=False) + "\n")
            accepted.flush()
            print(f"{case['id']}: {result['status']} ({result['accepted_count']} accepted, {result['elapsed_s']:.2f}s)", flush=True)
    with (output / "summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        columns = ["case_id", "status", "accepted_count", "elapsed_s", "expected", "expectation_met"]
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    failed = [r for r in results if r["expectation_met"] is False or r["status"] in
              {"INPUT_ERROR", "OUTPUT_ERROR", "PROCESS_ERROR", "COARSE_ERROR", "PROCESS_TIMEOUT", "FINE_TIMEOUT"}]
    report = ["# Integration test results", "", "FINE_ACCEPTED = supplied fine planner accepted the generated path. No global optimality or independent CAD certification is claimed.", "",
              "| Case | Status | Accepted | Expected match | Seconds |", "|---|---|---:|---|---:|"]
    report += [f"| {r['case_id']} | {r['status']} | {r['accepted_count']} | {r['expectation_met']} | {r['elapsed_s']:.2f} |" for r in results]
    report += ["", f"Attempts: {len(results)}; accepted pairs: {sum(r['accepted_count'] for r in results)}; test/errors: {len(failed)}."]
    (output / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
