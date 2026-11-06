# SPDX-License-Identifier: GPL-3.0-or-later
"""Command line interface for StaphMix."""

import os
import sys
from datetime import datetime

import click

from . import __version__, configure
from .configure import (check_tools, db_path, default_db_folder,
                        missing_db_files, set_db_folder)
from .read_seq import process_input_file
from .run_amr import run_amrfinder
from .run_cgMLST import run_cgmlst
from .run_print import print_name
from .run_scc import run_scc_blastn
from .run_score import calculate_amr, calculate_virulence
from .run_species import species_identify
from .run_vfdb import runVFDB
from .write_result import (write_amr, write_cgMLST, write_scc, write_summary,
                           write_vf)

SEPARATOR = '=' * 89


def _resolve_database(db):
    """Set the database directory, raising friendly CLI errors."""
    folder = db or default_db_folder()
    if folder is None:
        raise click.ClickException(
            'StaphMix database not found. Download it and pass --db PATH, or '
            'set the STAPHMIX_DB environment variable.')
    try:
        set_db_folder(folder)
    except ValueError as exc:
        raise click.ClickException(str(exc))
    missing = missing_db_files(folder)
    if missing:
        raise click.ClickException(
            'Database directory {0} is missing: {1}'.format(folder, ', '.join(missing)))


@click.command()
@click.option('-q', '--query', required=True,
              help='Input genome file (fasta/fastq, optionally gzipped).')
@click.option('-s', '--species', is_flag=True, help='Run species identification.')
@click.option('-c', '--cgmlst', is_flag=True,
              help='Run cgMLST typing. The species must be Staphylococcus aureus.')
@click.option('-m', '--sccmec', is_flag=True, help='Run SCCmec typing.')
@click.option('-v', '--vf', is_flag=True, help='Run virulence gene comparison.')
@click.option('-a', '--amr', is_flag=True, help='Run antimicrobial resistance gene comparison.')
@click.option('-p', '--prefix', required=True, help='Prefix for output files.')
@click.option('-o', '--outdir', required=True, help='Output directory for results.')
@click.option('-n', '--n_thread', default=8, type=int, help='Number of threads to use.')
@click.option('--db', default=None,
              help='Database directory. Defaults to $STAPHMIX_DB, the conda '
                   'share directory, or ./db next to the source checkout.')
@click.version_option(__version__, prog_name='staphmix')
def main(query, species, cgmlst, sccmec, vf, amr, prefix, outdir, n_thread, db):
    print_name()
    try:
        check_tools()
    except RuntimeError as exc:
        raise click.ClickException(str(exc))
    _resolve_database(db)

    print('Database directory : {0}'.format(configure.db_folder))
    outdir_abpath = os.path.abspath(outdir)
    if not os.path.exists(outdir_abpath):
        os.makedirs(outdir_abpath, exist_ok=True)
        print('The working directory does not exist. Created directory: {0}'.format(outdir_abpath))
    else:
        print('The working directory already exists: {0}'.format(outdir_abpath))

    fasta = process_input_file(query, outdir_abpath)
    summary_list = [prefix]

    print(SEPARATOR)
    try:
        if species:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('\033[33m{0} Running ANI identification of Staphylococcus aureus species\033[0m'.format(current_time))
            species_out = species_identify(fasta, prefix, outdir_abpath, n_thread)
            if species_out[0]:
                print('Staphylococcus aureus : \033[31mTrue\033[0m  | Reference strain GCF000025145 '
                      'with ANI value : \033[31m{0}\033[0m\n{1}'.format(species_out[1], SEPARATOR))
            else:
                print('Staphylococcus aureus : \033[31mFalse\033[0m | Reference strain GCF000025145 '
                      'with ANI value : \033[31m{0}\033[0m\n{1}'.format(species_out[1], SEPARATOR))
                sys.exit()

        if cgmlst:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('\033[33m{0} Running cgMLST typing and hierarchical clustering\033[0m'.format(current_time))
            show_data, alleles, hiercc = run_cgmlst(fasta, n_thread)
            write_cgMLST(outdir_abpath, prefix, alleles)
            summary_list.append(hiercc['Hierarchical_cluster_level - HC6'])
            summary_list.append(hiercc['Hierarchical_cluster_level - HC650'])
            print('\n'.join(show_data))
            print(SEPARATOR)
        else:
            summary_list.append('None')
            summary_list.append('None')

        if sccmec:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('\033[33m{0} Running Staphylococcal chromosomal cassettes mec (SCCmec) typing\033[0m'.format(current_time))
            scc_type_print, outputs = run_scc_blastn(fasta, db_path('scc_type_db'), n_thread)
            if scc_type_print:
                write_scc(outdir_abpath, prefix, outputs)
                summary_list.append(scc_type_print)
                print('Staphylococcal chromosomal cassettes mec (SCCmec) type : '
                      '\033[31m{0}\033[0m\n{1}'.format(scc_type_print, SEPARATOR))
            else:
                summary_list.append('None')
                print('Staphylococcal chromosomal cassettes mec (SCCmec) type : \033[31mNone\033[0m\n'
                      + SEPARATOR)

        if vf:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('\033[33m{0} Running virulence genes comparison\033[0m'.format(current_time))
            vf_out = runVFDB(fasta, db_path('vfdb'), n_thread)
            if vf_out:
                write_vf(outdir_abpath, prefix, vf_out)
                vf_score_out = calculate_virulence(vf_out)
                summary_list.append(vf_score_out[0])
                summary_list.append(vf_score_out[1])
                print('Enterotoxin genes score: \033[31m{0}\033[0m\n'
                      'Blood infection genes score: \033[31m{1}\033[0m\n{2}'.format(
                          vf_score_out[0], vf_score_out[1], SEPARATOR))
            else:
                summary_list.extend(['0', '0'])
                print('The strain did not predict any virulence genes\n' + SEPARATOR)
        else:
            summary_list.append('None')
            summary_list.append('None')

        if amr:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('\033[33m{0} Running antibiotic resistance genes comparison\033[0m'.format(current_time))
            amr_out, amr_class, amr_num = run_amrfinder(fasta, n_thread)
            if amr_out:
                write_amr(outdir_abpath, prefix, amr_out)
                mec, amr_score, special_gene = calculate_amr(amr_out)
                if mec:
                    summary_list.extend(['MRSA', amr_score])
                    if special_gene:
                        summary_list.extend([amr_class, amr_num, ','.join(sorted(special_gene))])
                        print('The strain is \033[31mMRSA\033[0m\n'
                              'Antibiotic resistance genes score: \033[31m{0}\033[0m\n'
                              'Antibiotic resistance genes class: \033[31m{1}\033[0m\n'
                              'Antibiotic resistance genes number: \033[31m{2}\033[0m\n'
                              '\x1b[31mWarning: The strain may have resistance to {3}\x1b[0m\n{4}'.format(
                                  amr_score, amr_class, amr_num,
                                  ','.join(sorted(special_gene)), SEPARATOR))
                    else:
                        summary_list.extend([amr_class, amr_num, 'None'])
                        print('This strain is \033[31mMRSA\033[0m\n'
                              'Antibiotic resistance genes score: \033[31m{0}\033[0m\n'
                              'Antibiotic resistance genes class:  \033[31m{1}\033[0m\n'
                              'Antibiotic resistance genes number: \033[31m{2}\033[0m\n{3}'.format(
                                  amr_score, amr_class, amr_num, SEPARATOR))
                else:
                    summary_list.extend(['MSSA', amr_score, amr_class, amr_num, 'None'])
                    print('This strain is \033[31mMSSA\033[0m\n'
                          'Antibiotic resistance genes score: \033[31m{0}\033[0m\n'
                          'Antibiotic resistance genes class:  \033[31m{1}\033[0m\n'
                          'Antibiotic resistance genes number: \033[31m{2}\033[0m\n{3}'.format(
                              amr_score, amr_class, amr_num, SEPARATOR))
            else:
                summary_list.extend(['MSSA', 'None', 'None', 'None', 'None'])
                print('The strain did not predict any resistance genes\n' + SEPARATOR)
        else:
            summary_list.extend(['None'] * 5)

        summary_row = [str(value) for value in summary_list]
        write_summary(outdir_abpath, prefix, summary_row)
    finally:
        if os.path.exists(fasta):
            os.remove(fasta)


if __name__ == '__main__':
    main()
