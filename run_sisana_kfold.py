#!/usr/bin/env python3
"""
Run SISANA preprocess + generate across Set_* folds using a template params.yml.

Per run, this script updates:
- preprocess.exp_file
- preprocess.outdir
- generate.exp
- generate.pandafilepath
- generate.lionessfilepath
- generate.compute
- generate.ncores (forced to 1 when compute=gpu)

It can either:
1) only write params.yml files, or
2) write params.yml files and run both:
   - sisana preprocess <params.yml>
   - sisana generate <params.yml>

Example usage (write params only):
python /storage/kuijjerarea/plos/net_reg_pipeline/run_sisana_kfold.py \
  --kfold-dir /storage/kuijjerarea/plos/kFold_splits_ProtCod \
  --template /storage/kuijjerarea/plos/sisana160/example_inputs/params.yml \
  --out-root /storage/kuijjerarea/plos/sisana_runs \
  --include-test \
  --compute cpu

Example usage (write params and run):
python /storage/kuijjerarea/plos/net_reg_pipeline/run_sisana_kfold.py \
  --kfold-dir /storage/kuijjerarea/plos/kFold_splits_ProtCod \
  --template /storage/kuijjerarea/plos/sisana160/example_inputs/params.yml \
  --out-root /storage/kuijjerarea/plos/sisana_runs \
  --include-test \
  --compute gpu \
  --run \
  --conda-env sisana_env
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import yaml  # pip install pyyaml


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: The parsed command-line arguments.
    """
    p = argparse.ArgumentParser(description="Run SISANA preprocess+generate across k-fold splits.")
    p.add_argument("--kfold-dir", required=True, help="Folder containing Set_* directories")
    p.add_argument("--template", required=True, help="Template params.yml")
    p.add_argument("--out-root", default="./sisana_runs", help="Output root for all runs")
    p.add_argument("--include-test", action="store_true", help="Also run test split")
    p.add_argument("--run", action="store_true", help="Execute sisana commands")
    p.add_argument("--conda-env", default=None, help="Conda environment name (required with --run)")
    p.add_argument("--sisana-cmd", default="sisana", help="SISANA executable name")
    p.add_argument(
        "--generate-exp-mode",
        choices=["preprocess_output", "raw_input"],
        default="preprocess_output",
        help="Set generate.exp to preprocess output path or raw input expression file",
    )
    p.add_argument(
        "--compute",
        choices=["cpu", "gpu"],
        default="cpu",
        help="Set generate.compute. If gpu, ncores is forced to 1.",
    )
    return p.parse_args()


def load_yaml(path: Path) -> dict:
    """
    Load data from a YAML file.

    Args:
        path (Path): The path to the YAML file to be loaded.

    Returns:
        dict: The data loaded from the YAML file.

    Raises:
        ValueError: If the loaded YAML contains anything other than a mapping/dictionary.
    """
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Template YAML must be a mapping/object: {path}")
    return data


def save_yaml(path: Path, data: dict) -> None:
    """
    Save a dictionary to a YAML file, creating parent directories if necessary.

    Args:
        path (Path): The destination path for the YAML file.
        data (dict): The dictionary to write out to the file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def run_cmd(cmd: list[str], log_path: Path) -> int:
    """
    Run a subprocess command and write its standard output and error to a log file.

    Args:
        cmd (list[str]): The command and its arguments to execute.
        log_path (Path): The file path where the output will be logged.

    Returns:
        int: The return code of the executed command.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


def expected_preprocessed_file(preprocess_exp_file: Path, preprocess_outdir: Path) -> Path:
    """
    Determine the expected file path of the SISANA preprocessed output.

    Args:
        preprocess_exp_file (Path): The initial expression file path.
        preprocess_outdir (Path): The directory where the preprocessed file is saved.

    Returns:
        Path: The full expected path of the preprocessed file.
    """
    # Expected SISANA preprocess output naming convention from  template
    # e.g. input expression.tsv -> expression_preprocessed.txt
    return preprocess_outdir / f"{preprocess_exp_file.stem}_preprocessed.txt"


def main() -> int:
    """
    The main execution workflow. 
    Iterates through k-fold splits, prepares SISANA parameters, and optionally runs the SISANA steps.

    Returns:
        int: Exit status code (0 for success, 1 for errors).
    """
    args = parse_args()

    kfold_dir = Path(args.kfold_dir).resolve()
    template_path = Path(args.template).resolve()
    out_root = Path(args.out_root).resolve()

    if not kfold_dir.is_dir():
        print(f"ERROR: missing kfold dir: {kfold_dir}", file=sys.stderr)
        return 1
    if not template_path.is_file():
        print(f"ERROR: missing template file: {template_path}", file=sys.stderr)
        return 1
    if args.run and not args.conda_env:
        print("ERROR: --conda-env is required when --run is set", file=sys.stderr)
        return 1

    template = load_yaml(template_path)
    out_root.mkdir(parents=True, exist_ok=True)

    set_dirs = sorted([d for d in kfold_dir.glob("Set_*") if d.is_dir()])
    if not set_dirs:
        print(f"ERROR: no Set_* directories found in {kfold_dir}", file=sys.stderr)
        return 1

    splits = [("train", "Train_Set_VSTnorm_Counts.tsv")]
    if args.include_test:
        splits.append(("test", "Test_Set_VSTnorm_Counts.tsv"))

    manifest_rows: list[list[str]] = []

    for set_dir in set_dirs:
        for split_name, exp_filename in splits:
            exp_file = set_dir / "preprocessing_output" / exp_filename

            run_dir = out_root / set_dir.name / split_name
            preprocess_outdir = run_dir / "preprocess"
            network_dir = run_dir / "network"
            params_path = run_dir / "params.yml"

            if not exp_file.exists():
                manifest_rows.append([set_dir.name, split_name, str(params_path), "missing_expression"])
                continue

            params = deepcopy(template)
            params.setdefault("preprocess", {})
            params.setdefault("generate", {})

            # preprocess
            params["preprocess"]["exp_file"] = str(exp_file)
            params["preprocess"]["outdir"] = str(preprocess_outdir)

            # generate
            if args.generate_exp_mode == "preprocess_output":
                gen_exp = expected_preprocessed_file(exp_file, preprocess_outdir)
            else:
                gen_exp = exp_file

            params["generate"]["exp"] = str(gen_exp)
            params["generate"]["pandafilepath"] = str(network_dir / "panda_network.txt")
            params["generate"]["lionessfilepath"] = str(network_dir / "lioness_networks.npy")
            params["generate"]["compute"] = args.compute

            # enforce ncores policy
            if args.compute == "gpu":
                params["generate"]["ncores"] = 1

            save_yaml(params_path, params)

            status = "params_written"
            if args.run:
                pre_log = run_dir / "sisana_preprocess.log"
                gen_log = run_dir / "sisana_generate.log"

                pre_cmd = ["conda", "run", "-n", args.conda_env, args.sisana_cmd, "preprocess", str(params_path)]
                pre_rc = run_cmd(pre_cmd, pre_log)

                if pre_rc == 0:
                    gen_cmd = ["conda", "run", "-n", args.conda_env, args.sisana_cmd, "generate", str(params_path)]
                    gen_rc = run_cmd(gen_cmd, gen_log)
                    status = "ok" if gen_rc == 0 else f"generate_failed({gen_rc})"
                else:
                    status = f"preprocess_failed({pre_rc})"

            manifest_rows.append([set_dir.name, split_name, str(params_path), status])

    manifest_path = out_root / "manifest.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["set_id", "split", "params_file", "status"])
        writer.writerows(manifest_rows)

    print(f"Done. Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())