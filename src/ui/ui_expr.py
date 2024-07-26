import argparse
from ..core import Expression, Processor

def parse_expr(p_expr:argparse.ArgumentParser):
    expr_inp = p_expr.add_argument_group('Input options')
    expr_inp.add_argument('-f', '--fq-dir', metavar='<dir>', help='Fastq directory')
    expr_inp.add_argument('-e', '--fq-ext', metavar='<file ext>', type=str, default='fastq.gz', help='File extension of fastq files(default: fastq.gz)')
    expr_inp.add_argument('-fa', '--ref-fa', metavar='<path>', help='Reference genome fasta')
    expr_inp.add_argument('-anno', '--ref-anno', metavar='<path>', help='Reference genome annotation(gtf)')

    expr_out = p_expr.add_argument_group('Output options')
    expr_out.add_argument('-n', '--name', metavar='<name>', help='Output name')
    expr_out.add_argument('-o', '--output-dir', metavar='<dir>', help='Output directory')

    expr_star = p_expr.add_argument_group('STAR arguments')
    expr_star.add_argument('--staridx-dir', metavar='<dir>', help='STAR index directory')
    expr_star.add_argument('--star-processes', metavar='<num_processes>', type=int, default=5, help='Number of processes for mapping fastq files(default: 5)')
    expr_star.add_argument('-t', '--star-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for mapping RNA-Seq using STAR(default: 10)')
    expr_star.add_argument('--star-bam', metavar='<path>', help='Merged BAM files(optional)')

    expr_strg = p_expr.add_argument_group('Stringtie arguments')
    expr_strg.add_argument('--strg-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for assembling RNA-Seq using Stringtie(default: 10)')
    expr_strg.add_argument('--strg-gtf', metavar='<path>', help='Stringtie result GTF file(optional)')

    expr_gatk = p_expr.add_argument_group('GATK arguments')
    expr_gatk.add_argument('--gatk-rgprocesses', metavar='<num_processes>', type=int, default=5, help='Number of processes for adding readgroup for bam files(default: 5)')
    expr_gatk.add_argument('--gatk-mdprocesses', metavar='<num_processes>', type=int, default=5, help='Number of processes for marking duplicates for bam files(default: 5)')
    expr_gatk.add_argument('--gatk-sncprocesses', metavar='<num_processes>', type=int, default=5, help='Number of processes for spliting NCigar for bam files(default: 5)')
    expr_gatk.add_argument('--gatk-hcprocesses', metavar='<num_processes>', type=int, default=5, help='Number of processes for calling haplotype for bam files(default: 5)')
    expr_gatk.add_argument('-k', '--gatk-threads', metavar='<num_threads>', type=int, default=10, help='Number of threads for calling haplotype for bam files(default: 10)')
    expr_gatk.add_argument('--gatk-vcfprocesses', metavar='<num_processes>', type=int, default=5, help='Number of processes for filtering vcf files(default: 5)')
    expr_gatk.add_argument('--gatk-vcf', metavar='<path>', help='GATK result VCF file(optional)')
    expr_gatk.add_argument('--fixed-ref', metavar='<path>', help='Fixed reference FASTA file(optional)')


def run_expr(opts):
    expression = Expression(
        fq_dir=opts.fq_dir,
        fq_ext=opts.fq_ext,
        ref_fa=opts.ref_fa,
        ref_anno=opts.ref_anno,
        output_dir=opts.output_dir,
        name=opts.name,
        star_idx_dir=opts.staridx_dir,
        star_nprocesses=opts.star_processes,
        star_nthreads=opts.star_threads,
        star_bam=opts.star_bam,
        strg_nthreads=opts.strg_threads,
        strg_gtf=opts.strg_gtf,
        gatk_rg_nprocesses=opts.gatk_rgprocesses,
        gatk_rg_md_nprocesses=opts.gatk_mdprocesses,
        gatk_rg_md_snc_nprocesses=opts.gatk_sncprocesses,
        gatk_haplotype_caller_nprocesses=opts.gatk_hcprocesses,
        gatk_haplotype_caller_nthreads=opts.gatk_threads,
        gatk_filter_vcf_nprocesses=opts.gatk_vcfprocesses,
        gatk_vcf=opts.gatk_vcf,
        ref_fa_fixed=opts.fixed_ref,
    )
    print('<Step 1 / 3>')
    expression.proc_mapping_rna_seq()
    print('<Step 2 / 3>')
    ps1 = Processor.exec_proc(expression.proc_assemble_rna_seq)
    print('<Step 3 / 3>')
    ps2 = Processor.exec_proc(expression.proc_gatk_pipeline)
    ps1.join()
    ps2.join()
