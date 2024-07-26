import argparse
from ..core import util


def parse_util(p_util:argparse.ArgumentParser):
    util_sra = p_util.add_argument_group('SRA options')
    util_sra.add_argument('-sra', '--sra-dir', metavar='<dir>', help='SRA directory(optional)')
    util_sra.add_argument('-fq', '--fq-dir', metavar='<dir>', help='Fastq directory')
    util_sra.add_argument('--sra-processes', metavar='<num_processes>', type=int, default=5, help='Number of processes for processing SRA files(default: 5)')

    util_gff = p_util.add_argument_group('GFF options, coverting GFF to GTF')
    util_gff.add_argument('--gff', metavar='<path>', help='GFF file path')
    
    util_csv = p_util.add_argument_group('csv report options')
    util_csv.add_argument('-xn', '--csv-name', metavar='<name>', help='output name')
    util_csv.add_argument('-xo', '--csv-outdir', metavar='<dir>', help='Output directory')
    util_csv.add_argument('--cas-json', metavar='<path>', help='Cas setup')
    util_csv.add_argument('-xr', '--bg-rna', metavar='<path>', help='Background RNA fa')
    util_csv.add_argument('-xe', '--xelim', action='store_true', help='Generate report for elimination process')


# strg_gtf, strg_gtf_fpkm_filter, ref_fa_fixed, output_dir, bowtie_nthreads, suffix
    util_build_rna_index = p_util.add_argument_group('Build background RNA index')
    util_build_rna_index.add_argument('-is', '--strg-gtf', metavar='<path>', help='Stringtie result GTF file')
    util_build_rna_index.add_argument('--fpkm-filter', metavar='<value>', type=float, default=1.0, help='FPKM lower bound for filtering GTF output of Stringtie(default: 1.0)')
    util_build_rna_index.add_argument('--fixed-ref', metavar='<path>', help='Fixed reference FASTA file')
    util_build_rna_index.add_argument('-io', '--idx-outdir', metavar='<dir>', help='Output directory')
    util_build_rna_index.add_argument('-b', '--bowtie-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for Bowtie2(default: 10)')
    util_build_rna_index.add_argument('-sx', '--suffix-name', metavar='<suffix>', help='Suffix')

    

def run_util(opts):
    sra_dir=opts.sra_dir
    fq_dir=opts.fq_dir
    sra_nprocesses=opts.sra_processes
    gff_path = opts.gff
    csv_name = opts.csv_name
    csv_out = opts.csv_outdir
    csv_cas = opts.cas_json
    csv_bg_rna = opts.bg_rna
    csv_elim = opts.xelim
    if sra_dir is not None:
        util.proc_sra(sra_dir, sra_nprocesses, fq_dir)
    if gff_path is not None:
        util.proc_gff(gff_path)
    if csv_name is not None:
        util.proc_csv(csv_name, csv_out, csv_cas, csv_bg_rna, csv_elim)
    
    strg_gtf = opts.strg_gtf
    strg_gtf_fpkm_filter = opts.fpkm_filter
    ref_fa_fixed = opts.fixed_ref
    bidx_out = opts.idx_outdir
    bowtie_nthreads = opts.bowtie_threads
    bidx_suffix = opts.suffix_name
    if bidx_out is not None:
        util.proc_build_background_rna_index(strg_gtf, strg_gtf_fpkm_filter, ref_fa_fixed, bidx_out, bowtie_nthreads, bidx_suffix)

