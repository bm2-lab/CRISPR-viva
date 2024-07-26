import argparse
from ..core import Elimination

def parse_elim(p_elim:argparse.ArgumentParser):
    elim_inp = p_elim.add_argument_group('Input options')
    elim_inp.add_argument('-s', '--strg-gtf', metavar='<path>', help='Stringtie result GTF file')
    elim_inp.add_argument('--fpkm-filter', metavar='<value>', type=float, default=1.0, help='FPKM lower bound for filtering GTF output of Stringtie(default: 1.0)')
    elim_inp.add_argument('--fixed-ref', metavar='<path>', help='Fixed reference FASTA file')
    elim_inp.add_argument('-c', '--chrom-name', metavar='<chrom_name>', help='Chromosome name')

    elim_out = p_elim.add_argument_group('Output options')
    elim_out.add_argument('-n', '--name', metavar='<name>', help='Output name')
    elim_out.add_argument('-o', '--output-dir', metavar='<dir>', help='Output directory')
    elim_out.add_argument('-x', '--csv', action='store_true', help='Generate csv report or not')

    elim_grna = p_elim.add_argument_group('gRNA options')
    elim_grna.add_argument('-vgene', '--virus-genes', metavar='<gene_name>', help='Virus gene list', nargs='*')
    elim_grna.add_argument('--cas-json', metavar='<path>', help='Cas setup')
    elim_grna.add_argument('-m', '--mismatch', metavar='<nt>', type=int, default=6, help='Number of mismatches(default: 6)')

    elim_bw = p_elim.add_argument_group('Bowtie2 arguments')
    elim_bw.add_argument('--bowtie-processes', metavar='<num_processes>', type=int, default=5, help='Number of processes for Bowtie2(default: 5)')
    elim_bw.add_argument('-b', '--bowtie-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for Bowtie2(default: 10)')


def run_elim(opts):
    elim = Elimination(
        strg_gtf=opts.strg_gtf,
        strg_gtf_fpkm_filter=opts.fpkm_filter,
        ref_fa_fixed=opts.fixed_ref,
        chrom_name=opts.chrom_name,
        output_dir=opts.output_dir,
        bl_csv=opts.csv,
        name=opts.name,
        bowtie_nprocesses=opts.bowtie_processes,
        bowtie_nthreads=opts.bowtie_threads,
        virus_gene_selected=opts.virus_genes,
        cas_definition_json=opts.cas_json,
        grna_mismatch=opts.mismatch
    )

    elim.run()