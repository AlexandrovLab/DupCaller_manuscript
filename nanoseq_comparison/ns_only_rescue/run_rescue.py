#!/usr/bin/env python3
"""Run DupCaller --rescue over each sample's region file, with the exact
DupCaller that produced compare_20261117_dmgfix's <sample>_dmgfix runs (1.2.12
main 64dac19 + dmg uint8 fix, frozen DC4.call.rerun/code_v1212main_dmgfix_20261116/src
run under the DupCaller_test env with PYTHONPATH=src), same masks and
(no --skipCoveragePass: the 64dac19 code forbids it with --rescue), reusing that run's error profiles (-E) and round-1
mutation-rate tables (-mr). Copy of ns_only_rescue_investigation_20261010/run_rescue.py
with only the code path, _latest -> _dmgfix and the output dir changed; -nms is
left off here (the _dmgfix runs used it) so masked candidates are still reported
under --rescue."""
import os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

BENCH = "/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking"
OUT = os.path.join(BENCH, "ns_only_rescue_investigation_20261117_dmgfix")
REF = "/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/reference/hg38_gdc/dupcaller_index/GRCh38.d1.vd1.fa"
CODE = os.path.join(BENCH, "DC4.call.rerun/code_v1212main_dmgfix_20261116/src")
CONDA_SH = "/tscc/nfs/home/yuc211/miniforge3/etc/profile.d/conda.sh"
AA_SAMPLES = {"40mM-72h", "70mM-48h", "80mM-48h"}
THREADS_PER_SAMPLE = 4
MAX_CONCURRENT = 8


def bam_args(s):
    if s in AA_SAMPLES:
        return f"DC_mkdped/{s}_D.bam", f"DC_mkdped/{s}_UD.bam"
    return f"DC_mkdped/{s}.bam", f"DC_mkdped/{s}_n.bam"


def run_one(s):
    region_file = os.path.join(OUT, f"{s}_regions.bed.gz")
    bam, nbam = bam_args(s)
    out_dir = os.path.join(OUT, f"rescue_{s}")
    os.makedirs(out_dir, exist_ok=True)
    cmd = (
        f"source {CONDA_SH} && conda activate DupCaller_test && cd {BENCH} && "
        f"PYTHONPATH={CODE} python {CODE}/DupCaller.py call -b {bam} -n {nbam} -f {REF} "
        f"-m {BENCH}/snp_mask_dc_benchmark.bed.gz {BENCH}/noise_mask_dc_benchmark.bed.gz "
        f"-E {s}_dmgfix/ERROR/{s}_dmgfix -mr {s}_dmgfix/tmp/{s}_dmgfix "
        f"-R {region_file} --rescue -p {THREADS_PER_SAMPLE} "
        f"-o {out_dir}/rescue_{s}"
    )
    with open(os.path.join(OUT, f"{s}_rescue.log"), "w") as logfh:
        rc = subprocess.run(["bash", "-c", cmd], stdout=logfh, stderr=subprocess.STDOUT).returncode
    return s, rc


if __name__ == "__main__":
    samples = sys.argv[1:] or sorted(
        f[: -len("_regions.bed.gz")] for f in os.listdir(OUT) if f.endswith("_regions.bed.gz"))
    print(f"Running rescue for {len(samples)} samples: {samples}", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=MAX_CONCURRENT) as ex:
        for fut in as_completed([ex.submit(run_one, s) for s in samples]):
            s, rc = fut.result()
            print(f"[{'OK' if rc == 0 else f'FAILED rc={rc}'}] {s}", flush=True)
            results.append((s, rc))
    n_fail = sum(1 for _, rc in results if rc)
    print(f"\n{len(results) - n_fail}/{len(results)} succeeded")
    sys.exit(1 if n_fail else 0)
