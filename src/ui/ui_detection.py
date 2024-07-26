import argparse
from ..core import Detection

def parse_det(p_det:argparse.ArgumentParser):
    det_inp = p_det.add_argument_group('Input options')
    det_inp.add_argument('-s', '--strg-gtf', metavar='<path>', help='Stringtie result GTF file')
    det_inp.add_argument('--fpkm-filter', metavar='<value>', type=float, default=1.0, help='FPKM lower bound for filtering GTF output of Stringtie(default: 1.0)')
    det_inp.add_argument('--fixed-ref', metavar='<path>', help='Fixed reference FASTA file')
    det_inp.add_argument('-v', '--viral-gtf', metavar='<path>', help='Viral GTF file')
    det_inp.add_argument('-r', '--viral-fa', metavar='<path>', help='Viral Fasta file')
    det_inp.add_argument('-c', '--chrom-name', metavar='<chrom_name>', help='Chromosome name')

    det_out = p_det.add_argument_group('Output options')
    det_out.add_argument('-n', '--name', metavar='<name>', help='Output name')
    det_out.add_argument('-o', '--output-dir', metavar='<dir>', help='Output directory')
    det_out.add_argument('-x', '--csv', action='store_true', help='Generate csv report or not')

    det_grna = p_det.add_argument_group('gRNA options')
    det_grna.add_argument('-vgene', '--virus-genes', metavar='<gene_name>', help='Virus gene list', nargs='*')
    det_grna.add_argument('--cas-json', metavar='<path>', help='Cas setup')
    det_grna.add_argument('-m', '--mismatch', metavar='<nt>', type=int, default=6, help='Number of mismatches(default: 6)')

    det_bw = p_det.add_argument_group('Bowtie2 arguments')
    det_bw.add_argument('--bowtie-processes', metavar='<num_processes>', type=int, default=5, help='Number of processes for Bowtie2(default: 5)')
    det_bw.add_argument('-b', '--bowtie-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for Bowtie2(default: 10)')
    det_bw.add_argument('-bi', '--idx-prefix', metavar='<path>', help='Built bowtie index prefix')
    det_bw.add_argument('-br', '--bg-rna', metavar='<path>', help='Background RNA fa')


def run_det(opts):
    det = Detection(
        strg_gtf=opts.strg_gtf,
        strg_gtf_fpkm_filter=opts.fpkm_filter,
        ref_fa_fixed=opts.fixed_ref,
        viral_gtf=opts.viral_gtf,
        viral_ref_fa=opts.viral_fa,
        chrom_name=opts.chrom_name,
        output_dir=opts.output_dir,
        bl_csv=opts.csv,
        name=opts.name,
        bowtie_nprocesses=opts.bowtie_processes,
        bowtie_nthreads=opts.bowtie_threads,
        bowtie_idx_prefix = opts.idx_prefix,
        bowtie_bg_rna = opts.bg_rna,
        virus_gene_selected=opts.virus_genes,
        cas_definition_json=opts.cas_json,
        grna_mismatch=opts.mismatch
    )

    det.run()