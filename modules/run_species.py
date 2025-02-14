from .configure import db_list, exe
import os,subprocess

fastANI = exe['fastANI']
reference_genome = db_list['species_ref']

def species_identify(query_genome,prefix,outdir_abpath,threads):
    name = prefix + '.ani.txt'
    outfile = os.path.join(outdir_abpath,name)
    ani_cmd = f"{fastANI} -q {query_genome} -r {reference_genome} -o {outfile} -t {threads}"
    subprocess.run(ani_cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.path.getsize(outfile)==0:
        return [False, '<75% ']
    else:
        with open(outfile, 'r') as ani_file:
            for line in ani_file:
                part = line.strip().split('\t')
                ani = part[2]
                if float(ani) >= 95:
                    return [True, ani]
                else:
                    return [False, ani]


