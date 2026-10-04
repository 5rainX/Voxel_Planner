"""Generate small synthetic integration fixtures; dimensions are in voxels."""
import argparse
import json
from pathlib import Path
import random


def create_cases(destination, count=0, seed=20260927):
    destination.mkdir(parents=True, exist_ok=False)
    rng = random.Random(seed)
    cases = []
    common = ["--width", "2", "--thickness", "1", "--clearance", "0",
              "--min-begin-length", "1", "--min-end-length", "1",
              "--flat-bend-factor", "1", "--vertical-bend-factor", "1",
              "--twist-factor", "1", "--bend-factor", "1",
              "--start-normal", "0", "0", "1", "--end-normal", "0", "0", "1",
              "--start-tangent", "1", "0", "0", "--end-tangent", "1", "0", "0",
              "--no-bend45", "--heuristic", "3"]

    def add(case_id, blocked=(), start=(8, 16, 12), goal=(40, 16, 12), expected=None, extra=()):
        map_name = case_id + ".txt"
        with (destination / map_name).open("w", encoding="ascii") as stream:
            stream.write("map: 48 32 24\n")
            for x, y, z in sorted(blocked):
                stream.write(f"{x} {y} {z} 0\n")
        item = {"id": case_id, "map": map_name, "start": list(start), "goal": list(goal),
                "args": common + list(extra), "family": case_id.split("_")[0]}
        if expected:
            item["expected"] = expected
        cases.append(item)

    add("open_straight", expected="FINE_ACCEPTED")
    add("open_waypoint", expected="FINE_ACCEPTED", extra=("--via", "24", "16", "12", "1", "0", "0", "1"))
    add("blocked_wall", ((24, y, z) for y in range(32) for z in range(24)), expected="COARSE_NO_PATH")
    # The coarse planner intentionally has endpoint obstacle immunity. The
    # complete pipeline must reject this case, regardless of rejection stage.
    add("blocked_start", [(8, 16, 12)], expected="NO_SOLUTION")
    add("single_obstacle", ((x, y, z) for x in range(22, 26) for y in range(13, 20) for z in range(24)))
    # Random cases deliberately have no expected verdict. A planning failure
    # is recorded as a failure to find a solution, not a proof of infeasibility.
    for index in range(count):
        obstacles = set()
        for _ in range(rng.randint(1, 3)):
            x0, y0 = rng.randint(18, 29), rng.randint(7, 22)
            for x in range(x0, x0 + rng.randint(2, 4)):
                for y in range(y0, y0 + rng.randint(2, 5)):
                    for z in range(24):
                        obstacles.add((x, y, z))
        y_start = rng.randint(10, 21)
        add(f"random_{index:04d}", obstacles, (8, y_start, 12), (40, y_start, 12))
    manifest = {"schema_version": 1, "units": "voxel", "synthetic": True,
                "seed": seed, "cases": cases}
    path = destination / "cases.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--random-cases", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260927)
    args = parser.parse_args()
    if args.random_cases < 0:
        parser.error("--random-cases must be nonnegative")
    print(create_cases(args.output, args.random_cases, args.seed))
