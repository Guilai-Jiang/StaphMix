from shutil import which
import os,numpy as np

exe = dict(makeblastdb=which('makeblastdb'), 
           blastn=which('blastn'), 
           diamond=which('diamond'),
           amrfinder=which('amrfinder'),
           fastANI=which('fastANI'))

current_directory = os.getcwd()
dirname = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(dirname)
db_folder = os.path.join(parent_dir, 'db')

db_list = dict(refseqs = os.path.join(db_folder, 'MLSTdb.fasta'), 
          core_genes = os.path.join(db_folder, 'staph_aureus.cgMLST'), 
          repr = os.path.join(db_folder, 'new_profile.csv'),
          hcc = os.path.join(db_folder, 'new_Mtb.HierCC.HierCC'), 
          species_ref = os.path.join(db_folder, 'GCF_000025145.1_ASM2514v1_genomic.fna'), 
          vfdb = os.path.join(db_folder,'VFDB_setB_pro.fas'),
           scc_type_db = os.path.join(db_folder,'whole_cassette_SCCmec_database_EXTENDED_20171117.fasta'))

scheme_info = dict(min_iden = 0.65, min_frag=0.6, max_iden=0.95)

blosum62 = np.array(
    [4., -2., 0., -2., -1., -2., 0., -2., -1., 0., -1., -1., -1., -2., 0., -1., -1., -1., 1., 0., -4., 0.,
     -3., 0., -2., -1., 0., 0., 0., 0., 0., 0., -2., 4., -3., 4., 1., -3., -1., 0., -3., 0., 0., -4.,
     -3., 3., 0., -2., 0., -1., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., 0., -3.,
     9., -3., -4., -2., -3., -3., -1., 0., -3., -1., -1., -3., 0., -3., -3., -3., -1., -1., -4., -1., -2., -2.,
     -2., -3., 0., 0., 0., 0., 0., 0., -2., 4., -3., 6., 2., -3., -1., -1., -3., 0., -1., -4., -3., 1.,
     0., -1., 0., -2., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., -1., 1., -4., 2.,
     5., -3., -2., 0., -3., 0., 1., -3., -2., 0., 0., -1., 2., 0., 0., -1., -4., -2., -3., -1., -2., 4.,
     0., 0., 0., 0., 0., 0., -2., -3., -2., -3., -3., 6., -3., -1., 0., 0., -3., 0., 0., -3., 0., -4.,
     -3., -3., -2., -2., -4., -1., 1., -1., 3., -3., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -2., -3.,
     6., -2., -4., 0., -2., -4., -3., 0., 0., -2., -2., -2., 0., -2., -4., -3., -2., -1., -3., -2., 0., 0.,
     0., 0., 0., 0., -2., 0., -3., -1., 0., -1., -2., 8., -3., 0., -1., -3., -2., 1., 0., -2., 0., 0.,
     -1., -2., -4., -3., -2., -1., 2., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -3., 0., -4., -3.,
     4., 0., -3., 2., 1., -3., 0., -3., -3., -3., -2., -1., -4., 3., -3., -1., -1., -3., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., -3., -1., 1., -3., -2., -1., -3., 0.,
     5., -2., -1., 0., 0., -1., 1., 2., 0., -1., -4., -2., -3., -1., -2., 1., 0., 0., 0., 0., 0., 0.,
     -1., -4., -1., -4., -3., 0., -4., -3., 2., 0., -2., 4., 2., -3., 0., -3., -2., -2., -2., -1., -4., 1.,
     -2., -1., -1., -3., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -2., 0., -3., -2., 1., 0., -1., 2.,
     5., -2., 0., -2., 0., -1., -1., -1., -4., 1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., 3.,
     -3., 1., 0., -3., 0., 1., -3., 0., 0., -3., -2., 6., 0., -2., 0., 0., 1., 0., -4., -3., -4., -1.,
     -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., -2., -3., -1.,
     -1., -4., -2., -2., -3., 0., -1., -3., -2., -2., 0., 7., -1., -2., -1., -1., -4., -2., -4., -2., -3., -1.,
     0., 0., 0., 0., 0., 0., -1., 0., -3., 0., 2., -3., -2., 0., -3., 0., 1., -2., 0., 0., 0., -1.,
     5., 1., 0., -1., -4., -2., -2., -1., -1., 3., 0., 0., 0., 0., 0., 0., -1., -1., -3., -2., 0., -3.,
     -2., 0., -3., 0., 2., -2., -1., 0., 0., -2., 1., 5., -1., -1., -4., -3., -3., -1., -2., 0., 0., 0.,
     0., 0., 0., 0., 1., 0., -1., 0., 0., -2., 0., -1., -2., 0., 0., -2., -1., 1., 0., -1., 0., -1.,
     4., 1., -4., -2., -3., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., -1., -1., -1., -1., -2., -2., -2.,
     -1., 0., -1., -1., -1., 0., 0., -1., -1., -1., 1., 5., -4., 0., -2., 0., -2., -1., 0., 0., 0., 0.,
     0., 0., -4., -4., -4., -4., -4., -4., -4., -4., -4., 0., -4., -4., -4., -4., 0., -4., -4., -4., -4., -4.,
     1., -4., -4., -4., -4., -4., 0., 0., 0., 0., 0., 0., 0., -3., -1., -3., -2., -1., -3., -3., 3., 0.,
     -2., 1., 1., -3., 0., -2., -2., -3., -2., 0., -4., 4., -3., -1., -1., -2., 0., 0., 0., 0., 0., 0.,
     -3., -4., -2., -4., -3., 1., -2., -2., -3., 0., -3., -2., -1., -4., 0., -4., -2., -3., -3., -2., -4., -3.,
     11., -2., 2., -3., 0., 0., 0., 0., 0., 0., 0., -1., -2., -1., -1., -1., -1., -1., -1., 0., -1., -1.,
     -1., -1., 0., -2., -1., -1., 0., 0., -4., -1., -2., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., -3.,
     -2., -3., -2., 3., -3., 2., -1., 0., -2., -1., -1., -2., 0., -3., -1., -2., -2., -2., -4., -1., 2., -1.,
     7., -2., 0., 0., 0., 0., 0., 0., -1., 1., -3., 1., 4., -3., -2., 0., -3., 0., 1., -3., -1., 0.,
     0., -1., 3., 0., 0., -1., -4., -2., -3., -1., -2., 4., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])