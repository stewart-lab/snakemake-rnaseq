import os
import shutil
import shlex
import subprocess
import argparse
import hashlib
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('config_yaml', type=str)
args = parser.parse_args()

def main():
    config_yaml = os.path.normpath(args.config_yaml)
    config_dir = os.path.dirname(config_yaml)

    # check that config_yaml exists
    if not os.path.exists(config_yaml):
        print(f"Error: config file '{config_yaml}' does not exist")
        exit(1)

    # check that config_yaml ends with .yaml or .yml
    if not str(config_yaml).endswith(".yaml") and not str(config_yaml).endswith(".yml"):
        print(f"Error: config file '{config_yaml}' must end with .yaml or .yml")
        exit(1)

    # read the config yaml
    with open(config_yaml) as f:
        config = yaml.safe_load(f) or dict()
    cores = config.get("cores", 1)
    results_dir = config.get("results_dir")
    workflow_dir = config.get("workflow_dir")
    
    if not workflow_dir:
        print("Error: config must specify 'workflow_dir'")
        exit(1)
    if not results_dir:
        print("Error: config must specify 'results_dir'")
        exit(1)

    # get workflow + results as absolute paths
    workflow_dir = os.path.normpath(os.path.join(config_dir, workflow_dir))
    results_dir = os.path.normpath(os.path.join(config_dir, results_dir))

    if not os.path.exists(workflow_dir):
        os.makedirs(workflow_dir)

    # copy the workflow into workflow_dir: copy envs/rules/scripts from the
    # container, replacing any file whose contents differ (checksum mismatch)
    # and leaving unchanged files in place.
    subdirs = ['envs', 'rules', 'scripts']
    for subdir in subdirs:
        src_subdir = os.path.join("/app/workflow", subdir)
        dst_subdir = os.path.join(workflow_dir, subdir)
        os.makedirs(dst_subdir, exist_ok=True)

        for src_file in os.listdir(src_subdir):
            src_path = os.path.join(src_subdir, src_file)
            dst_path = os.path.join(dst_subdir, src_file)

            if not os.path.isfile(src_path):
                continue  # skip nested dirs

            if os.path.exists(dst_path):
                if file_checksum(src_path) == file_checksum(dst_path):
                    continue  # unchanged, leave it
                os.remove(dst_path)  # optional; copy2 overwrites anyway

            shutil.copy2(src_path, dst_path)


    # always refresh the Snakefile from the container
    if os.path.exists(os.path.join(workflow_dir, "Snakefile")):
        os.remove(os.path.join(workflow_dir, "Snakefile"))
    
    shutil.copy2(os.path.join("/app/workflow/Snakefile"), os.path.join(workflow_dir, "Snakefile"))
    
    # change working directory to workflow dir
    os.chdir(workflow_dir)

    # write DAG image (best-effort: a failure here shouldn't abort the run)
    if results_dir:
        dag_dir = os.path.join(results_dir, "000_meta")
        os.makedirs(dag_dir, exist_ok=True)
        dag_png = os.path.join(dag_dir, "dag.png")
        dag_cmd = (
            f"snakemake --rulegraph --configfile {shlex.quote(config_yaml)} | "
            f"dot -Tpng -Grankdir=TB -Gsplines=ortho > {shlex.quote(dag_png)}"
        )
        dag_result = subprocess.run(dag_cmd, shell=True)
        if dag_result.returncode != 0:
            print(f"Warning: failed to render DAG image (exit {dag_result.returncode})")

    # run snakemake in the workflow dir; propagate its exit code so failures are
    # visible to callers (CI, schedulers, the `docker run` exit status).
    result = subprocess.run([
        "snakemake",
        "--use-conda",
        "--configfile", config_yaml,
        "--cores", str(cores),
        "--forcerun", "render_quarto",
    ])
    if result.returncode != 0:
        print(f"Error: snakemake failed with exit code {result.returncode}")
        exit(result.returncode)



def file_checksum(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

if __name__ == "__main__":
    main()