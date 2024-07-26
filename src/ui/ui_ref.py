import argparse
from ..core import Ref


def parse_ref(p_ref:argparse.ArgumentParser):
    ref_genome = p_ref.add_argument_group('Building human-virus hybrid genome')
    ref_genome.add_argument('-fa', '--ref-fa', metavar='<path>', help='Human genome fasta')
    ref_genome.add_argument('-vfa', '--virus-fa', metavar='<path>', help='Virus genome fasta')
    ref_genome.add_argument('-anno', '--ref-anno', metavar='<path>', help='Human genome annotation(gtf)')
    ref_genome.add_argument('-vanno', '--virus-anno', metavar='<path>', help='Virus genome annotation[gtf/gff]')

    ref_star = p_ref.add_argument_group('STAR arguments')
    ref_star.add_argument('--gen-staridx', metavar='<flag 0/1>', type=int, default=0, help='Built STAR index or not (default: 0 [not])')
    ref_star.add_argument('--staridx-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for building STAR index(default: 10)')
    ref_star.add_argument('--sjdb-overhang', metavar='<length>', type=int, default=100, help='Length of sjdb overhang(default: 100)')

    ref_out = p_ref.add_argument_group('Output options')
    ref_out.add_argument('-o', '--output-dir', metavar='<dir>', help='Output directory')


def run_ref(opts):
    ref = Ref(
        ori_ref_fa=opts.ref_fa,
        virus_fa=opts.virus_fa,
        ori_ref_anno=opts.ref_anno,
        virus_anno=opts.virus_anno,
        output_dir=opts.output_dir,
        staridx_flag=opts.gen_staridx,
        star_idx_nthreads=opts.staridx_threads,
        star_sjdb_overhang=opts.sjdb_overhang
    )
    ref.run()