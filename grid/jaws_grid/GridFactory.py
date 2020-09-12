import subprocess

class Slurm(object):
    def submit(self, params):
        """
        sbatch
        """
        cmd=["sbatch"]
        if "account" in params:
            cmd.extend(["-A", params['account']])
        if "constraint" in params:
            cmd.extend(["-C", params['constraint']])
        if "cpu" in params:
            cmd.extend(["-c", params['cpu']])
        if "exclusive" in params:
            cmd.append("--exclusive")
        if "name" in params:
            cmd.extend(["-N", params['name']])
        if "memory" in params:
            cmd.append(f"--mem={params['memory']}")
        elif "mem-per-cpu" in params:
            cmd.append(f"--mem-per-cpu={params['mem-per-cpu']}")
        if "output" in params:
            cmd.extend(["-o", params['output']])
        if "error" in params:
            cmd.extend(["-e", params['error']])
        if "qos" in params:
            cmd.extend(["-q", params['qos']])
        if "time" in params:
            cmd.extend(["-t", params['time']])
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except FileNotFoundError as error:
            print(f"FileNotFound: {error}")
        except subprocess.CalledProcessError as error:
            print(f"CalledProcessError: {error}")
        else:
            print(f"rc={result.returncode}")
            print(f"stdout:\n{result.stdout}")
            print(f"stderr:\n{result.stderr}")

    def cancel(self, job_id):
        """
        scancel
        """
        cmd=["scancel", str(job_id)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)

    def queue(self, job_id):
        """
        squeue
        """
        cmd=f"squeue -j {job_id} -o '%T'"
        # TODO
        

class GridFactory():
    def create_grid(self, typ):
        targetclass = typ.capitalize()
        return globals()[targetclass]()

