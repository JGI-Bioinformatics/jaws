class Slurm(object):
    def submit(self, params):
        """
        sbatch
        """
        options=["sbatch"]
        if "account" in params:
            options.append(f"-A {params['account']}")
        if "constraint" in params:
            options.append(f"-C {params['constraint']}")
        if "cpu" in params:
            options.append(f"-c {params['cpu']}")
        if "exclusive" in params:
            options.append("--exclusive")
        if "name" in params:
            options.append(f"-N {params['name']}")
        if "memory" in params:
            options.append(f"--mem={params['memory']}")
        elif "mem-per-cpu" in params:
            options.append(f"--mem-per-cpu={params['mem-per-cpu']}")
        if "output" in params:
            options.append(f"-o {params['output']}")
        if "error" in params:
            options.append(f"-e {params['error']}")
        if "qos" in params:
            options.append(f"-q {params['qos']}")
        if "time" in params:
            options.append(f"-t {params['time']}")
        cmd=' '.join(options)
        # TODO

    def cancel(self, job_id):
        """
        scancel
        """
        cmd=f"scancel {job_id}" 
        # TODO

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

