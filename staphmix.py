'''python version > 3.6'''
#!/usr/bin/env python3

import os,click,sys
from datetime import datetime
from modules.configure import current_directory,db_folder,db_list
from modules.run_species import species_identify
from modules.run_scc import run_scc_blastn
from modules.read_seq import process_input_file
from modules.run_vfdb import runVFDB
from modules.run_amr import run_amrfinder
from modules.run_print import print_name
from modules.run_cgMLST import run_cgmlst
from modules.run_score import calculate_virulence,calculate_amr
from modules.write_result import write_cgMLST,write_scc,write_vf,write_amr,write_summary



@click.command()
@click.option('-q', '--query', required=True, help='Path to the input genome file, support fasta(.gz)/fastq(.gz) file format.')
@click.option('-s', '--species', is_flag=True, help='Flag to run species identification.')
@click.option('-c', '--cgmlst', is_flag=True, help='Flag to run cgMLST typing. The species must be Staphylococcus aureus.')
@click.option('-m', '--sccmec', is_flag=True, help='Flag to run SCCmec typing.')
@click.option('-v', '--vf', is_flag=True, help='Flag to run virulence genes comparison.')
@click.option('-a', '--amr', is_flag=True, help='Flag to run antibiotic resistance genes comparison.')
@click.option('-p', '--prefix', required=True, help='Prefix for output files.')
@click.option('-o', '--outdir', required=True, help='Output directory for results.')
@click.option('-n', '--n_thread', default=8, type=int, help='Number of threads to use.')
def main(query, species, cgmlst,vf, amr, sccmec, prefix, outdir, n_thread):
    print_name()
    print(f'Current location : {current_directory}' )
    print(f'Database directory : {db_folder}')
    outdir_abpath = os.path.join(current_directory,outdir)
    if not os.path.exists(outdir_abpath):
        os.makedirs(outdir_abpath, exist_ok=True)
        print(f"The working directory does not exist. Created directory: {outdir_abpath}")
    else:
        print(f"The working directory already exists: {outdir_abpath}")


    fasta = process_input_file(query,outdir_abpath)
    #fasta = '/Terra/guilai/StaphMix/tmpbzy6fhs6.fasta'
    #fasta = '/titan/guilai/GCF_000025145.1_ASM2514v1_genomic.fna'
    #fasta = '/home/guilai/GCF_000025145.1_ASM2514v1_genomic.fna'
    summary_list = [prefix]


    print("=========================================================================================")
    if species:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"\033[33m{current_time} Running ANI identification of Staphylococcus aureus species\033[0m")
        species_out = species_identify(fasta,prefix,outdir_abpath,n_thread)
        if species_out[0]:
            print(f"Staphylococcus aureus : \033[31mTrue\033[0m  | Reference strain GCF000025145 with ANI value : \033[31m{species_out[1]}\033[0m\n"
                  "=========================================================================================")

        else:
            print(f"Staphylococcus aureus : \033[31mFalse\033[0m | Reference strain GCF000025145 with ANI value : \033[31m{species_out[1]}\033[0m\n"
                  f"=========================================================================================")
            sys.exit()
    if cgmlst:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'\033[33m{current_time} Running cgMLST typing and hierarchical clustering\033[0m')
        show_data,alleles, hiercc = run_cgmlst(fasta,n_thread)
        write_cgMLST(outdir_abpath,prefix,alleles)
        summary_list.append(hiercc['Hierarchical_cluster_level - HC6'])
        summary_list.append(hiercc['Hierarchical_cluster_level - HC650'])
        print('\n'.join(show_data))
        print("=========================================================================================")
    else:
        summary_list.append('None')
        summary_list.append('None')


    if sccmec:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'\033[33m{current_time} Running Staphylococcal chromosomal cassettes mec (SCCmec) typing\033[0m')
        scc_db = db_list['scc_type_db']
        scc_type_print,outputs = run_scc_blastn(fasta,scc_db,n_thread)
        if scc_type_print:
            write_scc(outdir_abpath,prefix,outputs)
            summary_list.append(scc_type_print)
            print(f"Staphylococcal chromosomal cassettes mec (SCCmec) tpye : \033[31m{scc_type_print}\033[0m\n"
                  "=========================================================================================")
        else:
            summary_list.append('None')
            print("Staphylococcal chromosomal cassettes mec (SCCmec) tpye : \033[31mNone\033[0m\n"
                  "=========================================================================================")

    if vf:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'\033[33m{current_time} Running virulence genes comparison\033[0m')
        vf_db = db_list['vfdb']
        vf_out = runVFDB(fasta,vf_db,n_thread)
        if vf_out:
            write_vf(outdir_abpath,prefix,vf_out)
            vf_score_out = calculate_virulence(vf_out)
            summary_list.append(vf_score_out[0])
            summary_list.append(vf_score_out[1])
            print(f"Enterotoxin genes score: \033[31m{vf_score_out[0]}\033[0m\n"
                  f"Blood infection genes score: \033[31m{vf_score_out[1]}\033[0m\n"
                  "=========================================================================================")
        else:
            summary_list.extend(['0','0'])
            print("The strain did not predict any virulence genes\n"
                  "=========================================================================================")
    else:
        summary_list.append('None')
        summary_list.append('None')      
    
    if amr:
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'\033[33m{current_time} Running antibiotic resistance genes comparison\033[0m')
        amr_out,amr_class,amr_num= run_amrfinder(fasta,n_thread)
        if amr_out:
            write_amr(outdir_abpath,prefix,amr_out)
            mec,amr_score,special_gene = calculate_amr(amr_out)
            if mec:
                summary_list.extend(['MRSA',amr_score])
                if special_gene:
                    summary_list.extend([amr_class,amr_num,','.join(special_gene),])
                    print("The strain is \033[31mMRSA\033[0m\n"
                        f"Antibiotic resistance genes score: \033[31m{amr_score}\033[0m\n"
                        f"Antibiotic resistance genes class: \033[31m{amr_class}\033[0m\n"
                        f"Antibiotic resistance genes number: \033[31m{amr_num}\033[0m\n"
                        f"\x1b[31mWarning: The strain may have resistance to {','.join(special_gene)}\x1b[0m\n"
                        "=========================================================================================")
                else:
                    summary_list.extend([amr_class,amr_num,'None'])
                    print("This strain is \033[31mMRSA\033[0m\n"
                        f"Antibiotic resistance genes score: \033[31m{amr_score}\033[0m\n"
                        f"Antibiotic resistance genes class:  \033[31m{amr_class}\033[0m\n"
                        f"Antibiotic resistance genes number: \033[31m{amr_num}\033[0m\n"
                        "=========================================================================================")
            else:
                summary_list.extend(['MSSA',amr_score,amr_class,amr_num,'None'])
                print("This strain is \033[31mMSSA\033[0m\n"
                      f"Antibiotic resistance genes score: \033[31m{amr_score}\033[0m\n"
                      f"Antibiotic resistance genes class:  \033[31m{amr_class}\033[0m\n"
                      f"Antibiotic resistance genes number: \033[31m{amr_num}\033[0m\n"
                      "=========================================================================================")
        else:
            summary_list.extend(['MSSA','None','None','None','None'])
            print("The strain did not predict any resistance genes\n"
                   "=========================================================================================")
    else:
        summary_list.append('None')
        summary_list.append('None') 
        summary_list.append('None')
        summary_list.append('None') 
        summary_list.append('None')

    summary_row = [str(value) for value in summary_list]
    write_summary(outdir_abpath,prefix,summary_row)
    if os.path.exists(fasta):
        os.remove(fasta)







if __name__ == '__main__':
     main()




#main('/home/guilai/staph_aureus/staph_cgmlst/saureus_2927_fasta/seq/GCA_903813275.1_22276_7_164_genomic.fna.gz','spcies','cgmlst','sccmec','vf','amr','test','temp','6')
