import subprocess,tempfile
from .configure import current_directory,exe

amrfinder = exe['amrfinder']

def run_amrfinder(input,threads):
    with tempfile.NamedTemporaryFile(mode='w+', delete=True,dir=current_directory) as temp_file:
        output_name = temp_file.name
        amrfinder_cmd = f'{amrfinder} -n  {input} --threads {threads}  -o {output_name}'
        process = subprocess.run(amrfinder_cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, check=True)
        data = []
        amr_class = set()
        if process.returncode == 0 :
            with open(output_name) as fin :
                fin.readline()
                for line in fin :
                    part = line.strip().split('\t')
                    #print(part)
                    if part[8] == 'AMR':
                        data.append(dict(
                            locus_tag='-',
                            contig=part[1],
                            start=int(part[2]),
                            end=int(part[3]),
                            direct=part[4],
                            identity = part[16],
                            coverage = part[15],
                            accession=part[18],
                            gene=part[5],
                            category=part[8] if part[8] == part[9] else '{0}:{1}'.format(part[8], part[9]),
                            subcategory=part[10] if part[10] == part[11] else '{0}:{1}'.format(part[10], part[11]),
                            function=part[6],
                        ))
                        amr_class.add(part[10])
    #print(data,len(amr_class),len(data))
    return data,len(amr_class),len(data)


# if __name__ == '__main__':
#      run_amrfinder()
#run_amrfinder('/home/guilai/software/staphalytics/run/HC0304799.fasta','8')
