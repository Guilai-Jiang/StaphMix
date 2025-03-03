import numpy as np, pandas as pd,os
import hashlib
from collections import OrderedDict
import warnings
from .configure import db_list, scheme_info

from .uberBLAST import uberBlast, readFastq, rc

refseqs = db_list['refseqs']
core_genes = db_list['core_genes']
repr = db_list['repr']
hcc = db_list['hcc']



warnings.simplefilter(action='ignore', category=FutureWarning)

def get_md5(value, dtype=str):
    if value != '':
        m = hashlib.md5(str(value).encode()).hexdigest()
        if dtype == str:
            return m
        else:
            return int(m, 16)
    else:
        pass
    
def delete(value):
    if  value < 0:
        return ''
    else:
        return value



def run_cgmlst(query, n_thread) :
    alleles = nomenclature(query, n_thread)
    repr_profile = pd.read_csv(repr, sep='\t',na_values=['None'],dtype=str)
    #print(repr_profile)
    repr_profile = repr_profile.set_index(repr_profile.columns[0])

    repr_hcc = pd.read_csv(hcc, sep='\t')
    repr_hcc = repr_hcc.set_index(repr_hcc.columns[0])
    #print(repr_hcc)
    profile = np.array([('-' if v.startswith('-') else v) for v in [v.get('value_md5', '-') for k, v in sorted(alleles.items())]])
    
    #print(repr_profile.values == profile)
    
    relshare = (np.sum((repr_profile.values == profile) & (profile != ''), 1)*profile.size+0.1)/(np.max(( np.sum((repr_profile.values != '') & (profile != '-'), 1), \
                                                                                np.sum((repr_profile.values != ''), 1)*0.97 ), 0)+0.1)

    max_idx = np.argmax(relshare)
    min_dist = int(profile.size-relshare[max_idx]+0.5)
    ref_repr = repr_profile.index[max_idx][:13]
    hc = repr_hcc.loc[ref_repr].values.astype(str).tolist()
    hc[:min_dist] = ['ND'] * min_dist
    hiercc = {'Reference_genome':ref_repr, 'Allelic_distance':min_dist}
    hiercc['Hierarchical_cluster_level - HC650'] = hc[650]
    hiercc['Hierarchical_cluster_level - HC6'] = hc[6]

    show_data = []
    for key, value in hiercc.items():
        show_data.append(f"{key} : {value}")    
    return show_data, alleles, hiercc

    
def nomenclature(query, n_thread=8) :
    core = {}
    with open(core_genes, 'rt') as fin :
        for line in fin :
            core[line.strip().split()[0]] = 1
            
    blastab = uberBlast(
        '-r {0} -q {1} -f --blastn --diamond --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 21,21 -m --merge_gap 300'.format(
            query, refseqs, scheme_info['min_iden']-0.05, 0.05, n_thread).split())
    merges = {}
    for b in blastab.T[16] :
        if len(b) > 4 :
            key = tuple(b[3:])
            merges[key] = b[:3]
    bsn = { b[15]:b[:15] for b in blastab }
    for ids, (score, iden, size) in merges.items() :
        bs = np.array([ bsn[i] for i in ids ], dtype=object)
        if np.unique(bs.T[1]).size == 1 :
            bs[0][2], bs[0][3], bs[0][7], bs[0][9], bs[0][11], bs[0][14] = iden, size, bs[-1][7], bs[-1][9], score, 'COMPLEX'
        for i in ids[1:] :
            bsn.pop(i)
        bsn[ids[0]] = bs[0]
    #print(bsn)
    blastab = np.array(list(bsn.values()))

    blastab = blastab[(blastab.T[2] >= scheme_info['min_iden']) & ((blastab.T[7]-blastab.T[6]+1) >= scheme_info['min_frag'] * blastab.T[12])]
    blastab.T[11] = blastab.T[11]*(blastab.T[7] - blastab.T[6]+1)/blastab.T[12]
    blastab = blastab[np.lexsort((blastab[:, 8], blastab[:, 1]))]
    for i, b0 in enumerate(blastab[:-1]) :
        if b0[0] == '' :
            continue
        s0, e0 = sorted(b0[8:10])
        todel = []
        for b1 in blastab[i+1:] :
            s1, e1 = sorted(b1[8:10])
            if b0[1] != b1[1] or e0 < s1 :
                break
            ovl = min(e0, e1) - max(s0, s1) + 1
            if ovl >= 0.5 * (e0-s0+1) or ovl >= 0.5 * (e1-s1+1) :
                sc0, sc1 = abs(b0[11]), abs(b1[11])
                g0, g1 = b0[0].rsplit('_', 1)[0], b1[0].rsplit('_', 1)[0]
                if b0[2] < b1[2]*scheme_info['max_iden'] or (b1[2] >= b0[2]*scheme_info['max_iden']
                    and (sc0 < sc1 or (sc0 == sc1 and b0[0] > b1[0]))) :
                    b0[11] = -sc0
                    if g0 == g1 or sc0 < sc1 * scheme_info['max_iden'] or b0[2] < b1[2]*scheme_info['max_iden'] :
                        b0[0] = ''
                        break
                else :
                    b1[11] = -sc1
                    if g0 == g1 or sc1 < sc0 * scheme_info['max_iden'] or b1[2] < b0[2]*scheme_info['max_iden'] :
                        todel.append(b1)
        if b0[0] and len(todel) :
            for b1 in todel :
                b1[0] = ''
    blastab = blastab[blastab.T[0] != '']
    blastab = blastab[np.lexsort([-blastab.T[11], [b.rsplit('_', 1)[0] for b in blastab.T[0]]])]
    alleles = OrderedDict()
    for bsn in blastab :
        gene = bsn[0].rsplit('_', 1)[0]
        if gene in alleles :
            if alleles[gene]['score']*scheme_info['max_iden'] > bsn[11] :
                continue
            alleles[gene]['coordinates'].append((bsn[1], bsn[8], bsn[9]))
            alleles[gene]['flag'] |= 32
            continue
        flag = 0
        if bsn[6] > 1 or bsn[7] < bsn[12] :
            flag = 64
            if bsn[14] != 'COMPLEX' :
                if bsn[6] > 1 :
                    bsn[14] = '{0}D{1}'.format(bsn[6]-1, bsn[14])
                if bsn[7] < bsn[12] :
                    bsn[14] = '{0}{1}D'.format(bsn[14], bsn[12]-bsn[7])
        alleles[gene] = {'gene_name': gene, 'CIGAR':bsn[0]+':'+bsn[14], 'reference':os.path.basename(query),
                         'identity': bsn[2], 'coordinates':[(bsn[1], bsn[8], bsn[9])], 'flag':flag, 'score':bsn[11]}

    seq, qual = readFastq(query)
    for gene, allele in sorted(alleles.items()) :
        if allele['flag'] & 96 == 96 :
            alleles.pop(gene)
            continue
        if allele['flag'] & 32 > 0 :
            allele['sequence'] = 'DUPLICATED'
        else :
            c, s, e = allele['coordinates'][0]
            ss = seq[c][s - 1:e] if s < e else rc(seq[c][e - 1:s])
            qs = (min(qual[c][s - 1:e] if s < e else qual[c][e-1:s])) if len(qual) else 0
            
            allele['sequence'] = ss
            if qs < 10 :
                allele['flag'] |= 2
            if allele['flag'] == 0 :
                allele['flag'] = 1

        allele['coordinates'] = ','.join(['{0}:{1}..{2}'.format(*c) for c in  allele['coordinates']])
        allele['value_md5'] = ('' if allele['flag'] < 16 else '-') + get_md5(allele['sequence'])
    #print(alleles)
    return {g:alleles.get(g, {"gene_name":g, "value_md5":"-"}) for g in core }


# if __name__ == '__main__':
#     main()

#cgmlst('/Terra/guilai/staph_aureus/staph_cgmlst/saureus_2927_fasta/seq/GCA_903812305.1_22276_6_90_genomic.fna.gz')
