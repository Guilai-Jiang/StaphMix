import os


def write_cgMLST(outdir_abpath,prefix,alleles) :
    output_name = prefix + '.alleles.profile'
    outfile = os.path.join(outdir_abpath,output_name)
    with open(outfile, 'wt') as outfile :
        outfile.write('#Query\t{0}\n'.format('\t'.join([g for g in sorted(alleles.keys())])))
        outfile.write('{0}\t{1}\n'.format(outdir_abpath, '\t'.join([allele.get('value_md5', '-') for g, allele in sorted(alleles.items())])))


def write_scc(outdir_abpath,prefix,scc_out):
    output_name = prefix + '.sccmec.txt'
    header = [
    "Query contig", "Match ID", "Identity (%)", "Alignment Length", "Gap", 
    "Strand", "Start Position (query)", "End Position (query)", "E-value", 
    "Score", "Start Position (subject)", "End Position (subject)", "Bit Score"]
    outfile = os.path.join(outdir_abpath,output_name)
    with open(outfile, 'wt') as outfile :
        outfile.write('\t'.join(header)+'\n')
        for value in scc_out.values():
            row = value[3][0]
            row_str = [str(item) for item in row]
            outfile.write('\t'.join(row_str)+'\n')

def write_vf(outdir_abpath,prefix,vfdb_out):
    output_name = prefix + '.vf.txt'
    header = vfdb_out[0].keys()
    outfile = os.path.join(outdir_abpath,output_name)
    with open(outfile, 'wt') as outfile :
        outfile.write('\t'.join(header)+'\n')
        for entry in vfdb_out:
            row = [str(entry[key]) for key in header]
            outfile.write('\t'.join(row)+'\n')

def write_amr(outdir_abpath,prefix,amr_out):
    output_name = prefix + '.amr.txt'
    header = amr_out[0].keys() 
    outfile = os.path.join(outdir_abpath, output_name)
    with open(outfile, 'wt') as outfile:
        outfile.write('\t'.join(header) + '\n')  
        for entry in amr_out:
            row = [str(entry[key]) for key in header]  
            outfile.write('\t'.join(row) + '\n')  

def write_summary(outdir_abpath,prefix,summary_list):
    output_name = prefix + '.summary.txt'
    header = ['Query genome','HC6','HC650','SCCmec tpye', 'Enterotoxin score',
              'Blood infection score','MRSA/MSSA','Antibiotic resistance core','Antibiotic resistance gene_class','Antibiotic resistance gene num','Special antibiotic resistance']
    outfile = os.path.join(outdir_abpath,output_name)
    with open(outfile, 'wt') as outfile :
        outfile.write('\t'.join(header)+'\n')
        outfile.write('\t'.join(summary_list)+'\n')











# def process_blastn(db,query,min_cov,min_ident):
#     headers = ['contig','accession_number','identity','coverage','mismatch','gap','start','end','bit_score','gene']
#     with open(os.path.join(configure.dirname,'blast.txt'),'wt') as out:
#         out.write('\t'.join(headers)+'\n')
#         for hit in run_blast(db,query,min_cov,min_ident):
#             contig = hit.contig_name
#             accession_number = str(hit.gene_id.split('|')[1])
#             identity = str(hit.pident)
#             coverage = str(hit.ref_cov * 100)
#             mismatch = str(hit.mismatch_number)
#             gap = str(hit.gaps)
#             start = str(hit.contig_start)
#             end = str(hit.contig_end)
#             bit_score = str(hit.score)
#             gene = str(hit.gene_id.split('|')[4])
#             contents = [contig,accession_number,identity,coverage,mismatch,gap,start,end,bit_score,gene]
#             print(contents)
#             out.write('\t'.join(contents)+'\n')





