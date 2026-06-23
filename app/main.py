import os
import shutil
import argparse
import hashlib
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('config_yaml', type=str)
args = parser.parse_args()

def main():
    config_yaml = args.config_yaml
    config_dir = os.path.dirname(config_yaml)
    project_dir = os.path.dirname(config_dir)
    workflow_dir = os.path.join(project_dir, 'workflow')

    # check that config_yaml exists
    if not os.path.exists(config_yaml):
        print(f"Error: config file '{config_yaml}' does not exist")
        exit(1)

    # check that config_yaml ends with .yaml or .yml
    if not str(config_yaml).endswith(".yaml") and not str(config_yaml).endswith(".yml"):
        print(f"Error: config file '{config_yaml}' must end with .yaml or .yml")
        exit(1)

    # read the number of cores from the config (default to 1 if unset)
    with open(config_yaml) as f:
        config = yaml.safe_load(f) or {}
    cores = config.get("cores", 1)

    if not os.path.exists(workflow_dir):
        os.makedirs(workflow_dir)

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


    # remove existing workflow/envs, workflow/rules, workflow/scripts, and workflow/Snakemake
    if os.path.exists(os.path.join(workflow_dir, "Snakefile")):
        os.remove(os.path.join(workflow_dir, "Snakefile"))
    
    shutil.copy2(os.path.join("/app/workflow/Snakefile"), os.path.join(workflow_dir, "Snakefile"))
    
    # change working directory to workflow dir
    os.chdir(workflow_dir)

    # write DAG image
    os.system(f'snakemake --rulegraph --configfile {config_yaml} | dot -Tpng > {config_dir}/dag.png')

    # run snakemake in the workflow dir
    os.system(f'snakemake --use-conda --configfile {config_yaml} --cores {cores} --forcerun render_quarto')



def file_checksum(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

if __name__ == "__main__":
    main()