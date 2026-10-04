"""Build the real collaborator CLI against the supplied Voxel_Planner source."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    sys.stdout.reconfigure(errors="replace")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--busbar-source", type=Path, required=True)
    parser.add_argument("--voxel-source", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--generator", help="Optional CMake generator; uses CMake's default otherwise")
    args = parser.parse_args()
    build = args.build_dir.resolve()
    # A fresh directory preserves old build provenance and test evidence.
    if build.exists():
        parser.error("--build-dir must be a new directory")
    for source in (args.busbar_source, args.voxel_source):
        if not (source / "CMakeLists.txt").is_file():
            parser.error("Not a source directory: " + str(source))
    build.mkdir(parents=True)
    provenance = {"sources": {}, "generator": args.generator, "configuration": "Release"}
    for name, source in (("busbar", args.busbar_source), ("voxel", args.voxel_source)):
        source = source.resolve()
        snapshot = build / "sources" / name
        snapshot.mkdir(parents=True)
        shutil.copy2(source / "CMakeLists.txt", snapshot)
        for part in ("src", "include", "tests", "Module1_Main", "Module2_Morphology", "Module3_AStar"):
            if (source / part).is_dir():
                shutil.copytree(source / part, snapshot / part)
        hashes = {str(p.relative_to(snapshot)).replace("\\", "/"):
                  hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sorted(snapshot.rglob("*")) if p.is_file()}
        provenance["sources"][name] = {"original": str(source), "sha256": hashes}
    project = Path(__file__).resolve().parent
    provenance["integration_cmake_sha256"] = hashlib.sha256((project / "CMakeLists.txt").read_bytes()).hexdigest()
    configure = ["cmake", "-S", str(project), "-B", str(build),
                 "-DCMAKE_BUILD_TYPE=Release",
                 "-DBUSBAR_SOURCE_DIR=" + str(build / "sources/busbar"),
                 "-DVOXEL_PLANNER_SOURCE_DIR=" + str(build / "sources/voxel")]
    if args.generator:
        configure += ["-G", args.generator]
    compile_cmd = ["cmake", "--build", str(build), "--target", "hybrid_validation_cli",
                   "--config", "Release", "-j", str(args.jobs)]
    provenance["commands"] = [configure, compile_cmd]
    (build / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
    for name, command in (("configure", configure), ("build", compile_cmd)):
        print("Running " + name + "...", flush=True)
        with (build / (name + ".log")).open("wb") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            print((build / (name + ".log")).read_text(encoding="utf-8", errors="replace")[-8000:])
            return result.returncode
    executable = build / "bin" / ("hybrid_validation_cli.exe" if os.name == "nt" else "hybrid_validation_cli")
    provenance["toolchain"] = [line for line in (build / "CMakeCache.txt").read_text(encoding="utf-8").splitlines()
                               if line.startswith(("CMAKE_CXX_COMPILER:", "CMAKE_BUILD_TYPE:", "OpenMP_CXX_FLAGS:", "CMAKE_GENERATOR:"))]
    provenance["executable_sha256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
    (build / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
    print("Built: " + str(executable))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
