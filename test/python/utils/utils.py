import subprocess as sp

def run_snakemake(invocation_str):
        p = sp.Popen(invocation_str, shell = True, stdout = sp.PIPE, stderr = sp.PIPE)

        return p
