from .processor import Processor
import os
import glob

class Expression:
    def __init__(self,
                 fq_dir,
                 fq_ext,
                 ref_fa,
                 ref_anno,
                 output_dir,
                 name,
                 star_idx_dir,
                 star_nprocesses,
                 star_nthreads,
                 star_bam,
                 strg_nthreads,
                 strg_gtf,
                 gatk_rg_nprocesses,
                 gatk_rg_md_nprocesses,
                 gatk_rg_md_snc_nprocesses,
                 gatk_haplotype_caller_nprocesses,
                 gatk_haplotype_caller_nthreads,
                 gatk_filter_vcf_nprocesses,
                 gatk_vcf,
                 ref_fa_fixed):
        self.fq_dir = fq_dir
        self.fq_ext = fq_ext
        self.ref_fa = ref_fa
        self.ref_anno = ref_anno
        self.output_dir = output_dir
        self.star_idx_dir = star_idx_dir
        self.star_nprocesses = star_nprocesses
        self.star_nthreads = star_nthreads
        self.star_bam = star_bam
        self.strg_nthreads = strg_nthreads
        self.strg_gtf = strg_gtf
        self.gatk_rg_nprocesses = gatk_rg_nprocesses
        self.gatk_rg_md_nprocesses = gatk_rg_md_nprocesses
        self.gatk_rg_md_snc_nprocesses = gatk_rg_md_snc_nprocesses
        self.gatk_haplotype_caller_nprocesses = gatk_haplotype_caller_nprocesses
        self.gatk_haplotype_caller_nthreads = gatk_haplotype_caller_nthreads
        self.gatk_filter_vcf_nprocesses = gatk_filter_vcf_nprocesses
        self.gatk_vcf = gatk_vcf
        self.ref_fa_fixed = ref_fa_fixed
        self.name = name

        os.makedirs(self.output_dir, exist_ok=True)
        if len(glob.glob(f'{fq_dir}/*.{fq_ext}')) == 0:
            raise RuntimeError('FASTQ file cannot be found!')
        self.star_out_dir = f'{self.output_dir}/star_out'
        os.makedirs(self.star_out_dir, exist_ok=True)
        if self.star_bam is None:
            self.star_bam = f'{self.output_dir}/{self.name}_merged.bam'
            self.is_star_bam_processed = False
        else:
            self.is_star_bam_processed = True
        if strg_gtf is None:
            self.strg_gtf = f'{self.output_dir}/{self.name}_strg.gtf'
            self.is_strg_gtf_processed = False
        else:
            self.is_strg_gtf_processed = True
        self.gatk_out_dir = f'{self.output_dir}/gatk_out'
        os.makedirs(self.gatk_out_dir, exist_ok=True)
        if gatk_vcf is None:
            self.gatk_vcf = f'{self.output_dir}/{self.name}_merged.vcf'
            self.is_gatk_vcf_processed = False
        else:
            self.is_gatk_vcf_processed = True

        if ref_fa_fixed is None:
            ref_fa_name = os.path.splitext(os.path.basename(self.ref_fa))[0]
            self.ref_fa_fixed = f'{self.output_dir}/{ref_fa_name}_{self.name}.fa'
            self.is_ref_fa_fixed_processed = False
        else:
            self.is_ref_fa_fixed_processed = True
        
        

    def proc_mapping_rna_seq(self):
        bypassed = (self.is_star_bam_processed | self.is_strg_gtf_processed) & (self.is_gatk_vcf_processed |self.is_ref_fa_fixed_processed)
        if bypassed:
            print('No need to map RNA sequences')
            return
        commands = Processor.mapping_rna_seq(self.fq_dir, self.fq_ext, self.star_nthreads, self.star_idx_dir, self.star_out_dir)
        Processor.multiproc_command(self.star_nprocesses, commands)
    
    def proc_assemble_rna_seq(self):
        if self.is_strg_gtf_processed:
            print('No need to assemble RNA-Seq')
            return
        print('Merging bam files...')
        command = Processor.merge_star_bam(self.star_bam, self.star_out_dir)
        Processor.proc_command(command)

        print('Assembling merged bam...')
        command = Processor.assemble_rna_transcripts(self.ref_anno, self.star_bam, self.strg_gtf, self.strg_nthreads)
        Processor.proc_command(command)

    def proc_gatk_pipeline(self):
        bypassed = self.is_gatk_vcf_processed | self.is_ref_fa_fixed_processed
        if not bypassed:
            print('GATK adding readgroups...')
            commands = Processor.gatk_add_read_group(self.star_out_dir, self.gatk_out_dir, self.name)
            Processor.multiproc_command(self.gatk_rg_nprocesses, commands)
            print('GATK mark duplicates...')
            commands = Processor.gatk_mark_duplicates(self.gatk_out_dir)
            Processor.multiproc_command(self.gatk_rg_md_nprocesses, commands)
            print('GATK spliting NCigarReads...')
            commands = Processor.gatk_split_ncigar_reads(self.ref_fa, self.gatk_out_dir)
            Processor.multiproc_command(self.gatk_rg_md_snc_nprocesses, commands)
            print('GATK calling haplotype...')
            commands = Processor.gatk_haplotype_caller(self.ref_fa, self.gatk_out_dir, self.gatk_haplotype_caller_nthreads)
            Processor.multiproc_command(self.gatk_haplotype_caller_nprocesses, commands)
            print('GATK filtering vcf files...')
            commands = Processor.gatk_filter_vcf(self.gatk_out_dir)
            Processor.multiproc_command(self.gatk_filter_vcf_nprocesses, commands)
            print('Merging vcf files...')
            command = Processor.gatk_merge_vcf(self.gatk_out_dir, self.gatk_vcf)
            Processor.proc_command(command)
        if not self.is_ref_fa_fixed_processed:
            command = Processor.reset_ref_fa_with_gatk_vcf(self.ref_fa, self.gatk_vcf, self.ref_fa_fixed)
            Processor.proc_command(command)

