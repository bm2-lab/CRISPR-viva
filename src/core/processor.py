import enum
import subprocess
from multiprocessing import Pool, Process
import os
import glob
import shutil
from pyfaidx import Fasta
import numpy as np
import json
import re
import pandas as pd
import openpyxl as xls
import sys
import dask

class Processor:

    g_nucleotide_dict = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'W': 'AT',
        'S': 'CG',
        'K': 'GT',
        'M': 'AC',
        'R': 'AG',
        'Y': 'CT',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }

    g_default_cas_definition = {
        "Cas12a": {
            "pfs": "TTTV",
            "loc": "up",
            "min_len": 18,
            "max_len": 25
        },
        "Cas13d": {
            "pfs": None,
            "loc": None,
            "min_len": 22,
            "max_len": 30
        },
        "LwCas13a": {
            "pfs": "H",
            "loc": "down",
            "min_len": 22,
            "max_len": 30
        }
    }
    
    # g_cas_types = {
    #     'LwCas13a': ('H', 'down'),
    #     'Cas13': (None, None),
    #     'Cas12': ('TTTV', 'up')
    # }

    
    @classmethod
    def exec_proc(cls, f):
        '''This function will return process handle that must be joined outside manually!'''
        process = Process(target=f)
        process.start()
        return process

    @classmethod
    def proc_command(cls, command):
        if command is None:
            print('Pass')
            return 0
        print(f'Execuate command: {command}...')
        status = subprocess.call(command, executable='/bin/bash', shell=True)
        print(f'End command: {command}')
        return status

    @classmethod
    def multiproc_command(cls, num_processes, commands):
        if commands is None:
            print('Pass')
            return 0
        pool = Pool(processes=num_processes)
        results = []
        for command in commands:
            result = pool.apply_async(cls.proc_command, args=(command,))
            results.append(result)
        pool.close()
        pool.join()
        status_value = [r.get() for r in results]
        return status_value

    @classmethod
    def multiproc_exec(cls, num_processes, func, arg_lst):
        if func is None:
            return 0
        pool = Pool(processes=num_processes)
        results = []
        for args in arg_lst:
            result = pool.apply_async(func, args=args)
            results.append(result)
        pool.close()
        pool.join()
        results_value = [r.get() for r in results]
        return results_value
        


    @classmethod
    def dump_fq(cls, fq_dir, sra_dir):
        if sra_dir is None:
            return None
        os.makedirs(fq_dir, exist_ok=True)
        sra_files = glob.glob(f'{sra_dir}/*.sra')
        commands = [
            f'fastq-dump --split-3 {sra} -O {fq_dir} --gzip' for sra in sra_files]
        return commands

    @classmethod
    def concat_fasta(cls, ori_ref_fa, virus_fa, output_dir):
        ref_fa_name = os.path.splitext(os.path.basename(ori_ref_fa))[0]
        virus_fa_name = os.path.splitext(os.path.basename(virus_fa))[0]
        ref_fa = f'{output_dir}/{ref_fa_name}_{virus_fa_name}.fa'
        shutil.copy(ori_ref_fa, ref_fa)
        gn_virus_fa = (line.strip() for line in open(virus_fa) if line.strip() != '')
        with open(ref_fa, 'rb+') as ref_f:
            while True:
                ref_f.seek(-1, os.SEEK_END)
                char = next(ref_f)
                if char == b'\n':
                    ref_f.seek(-1, os.SEEK_END)
                    ref_f.truncate()
                else:
                    break
            ref_f.write(b'\n')
        with open(ref_fa, 'r+') as ref_f:
            ref_f.seek(0, os.SEEK_END)
            for line in gn_virus_fa:
                ref_f.write(f'{line}\n')
        commands = []
        command = f'samtools faidx {ref_fa}'
        commands.append(command)
        ref_dict =  os.path.splitext(ref_fa)[0] + '.dict'
        command = f'gatk CreateSequenceDictionary -R {ref_fa} -O {ref_dict}'
        commands.append(command)
        cls.multiproc_command(len(commands), commands)
        return ref_fa
    
    @classmethod
    def prepare_fasta(cls, ref_fa):
        commands = []
        command = f'samtools faidx {ref_fa}'
        commands.append(command)
        ref_dict =  os.path.splitext(ref_fa)[0] + '.dict'
        command = f'gatk CreateSequenceDictionary -R {ref_fa} -O {ref_dict}'
        commands.append(command)
        cls.multiproc_command(len(commands), commands)


    @classmethod
    def concat_gtf(cls, ori_ref_anno, virus_anno, output_dir):
        ref_anno_name = os.path.splitext(os.path.basename(ori_ref_anno))[0]
        virus_anno_name = os.path.splitext(os.path.basename(virus_anno))[0]
        ref_anno = f'{output_dir}/{ref_anno_name}_{virus_anno_name}.gtf'
        ext = os.path.splitext(virus_anno)[1]
        if ext == '.gtf':
            pass
        elif ext.startswith('.gff'):
            print(f'Converting {virus_anno} to gtf file.')
            filename = f'{os.path.splitext(os.path.basename(virus_anno))[0]}.gtf'
            filepath = f'{output_dir}/{filename}'
            command = f'gffread {virus_anno} -T -o {filepath}'
            cls.proc_command(command)
            virus_anno = filepath
        shutil.copy(ori_ref_anno, ref_anno)
        gn_virus_anno = (line.strip() for line in open(virus_anno) if ((line.strip() != '') and (not line.startswith('#'))))
        with open(ref_anno, 'rb+') as ref_f:
            while True:
                ref_f.seek(-1, os.SEEK_END)
                char = next(ref_f)
                if char == b'\n':
                    ref_f.seek(-1, os.SEEK_END)
                    ref_f.truncate()
                else:
                    break
            ref_f.write(b'\n')
        with open(ref_anno, 'r+') as ref_f:
            ref_f.seek(0, os.SEEK_END)
            for line in gn_virus_anno:
                ref_f.write(f'{line}\n')
        return ref_anno

    @classmethod
    def prepare_gtf(cls, ref_anno):
        ext = os.path.splitext(ref_anno)[1]
        if ext == '.gtf':
            pass
        elif ext.startswith('.gff'):
            print(f'Converting {ref_anno} to gtf file.')
            ref_gtf =  os.path.splitext(ref_anno)[0] + '.gtf'
            command = f'gffread {ref_anno} -T -o {ref_gtf}'
            cls.proc_command(command)

    @classmethod
    def check_ref_consistency(cls, ref_fa, ref_anno):
        fa = Fasta(ref_fa)
        fa_chroms = set(fa.keys())
        anno_chroms = {line.strip().split('\t')[0] for line in open(ref_anno) if ((line.strip() != '') and (not line.startswith('#')))}
        assert len(fa_chroms - anno_chroms) == 0, f'{fa_chroms-anno_chroms} do not have annotation.'
        print('Reference consistency: Yes.')
    
    @classmethod
    def build_star_index(cls, ref_fa, ref_anno, star_idx_dir, star_idx_nthreads, star_sjdb_overhang):
        command = f'STAR --runThreadN {star_idx_nthreads} --runMode genomeGenerate --genomeDir {star_idx_dir} --genomeFastaFiles {ref_fa} --sjdbGTFfile {ref_anno} --sjdbOverhang {star_sjdb_overhang}'
        return command
    
    @classmethod
    def mapping_rna_seq(cls, fq_dir, ext, star_nthreads, star_idx_dir, star_out_dir):
        fq_files = glob.glob(f'{fq_dir}/*.{ext}')
        assert len(fq_files) > 0, f'No {ext} file found.'
        pe_fq_files = glob.glob(f'{fq_dir}/*_[1-2].{ext}')
        se_fq_files = [fq for fq in fq_files if fq not in pe_fq_files]
        se_fq_files_prefix = [os.path.basename(
            se_fq)[:-(1+len(ext))] for se_fq in se_fq_files]
        pe_fq_files_prefix = list(
            {os.path.basename(pe_fq)[:-(3+len(ext))] for pe_fq in pe_fq_files})
        se_commands = []
        for se in se_fq_files_prefix:
            se_path = f'{fq_dir}/{se}.{ext}'
            command = f'STAR --runThreadN {star_nthreads} --genomeDir {star_idx_dir} --readFilesIn {se_path} --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {star_out_dir}/{se} --readFilesCommand zcat'
            se_commands.append(command)
        pe_commands = []
        for pe in pe_fq_files_prefix:
            pe1_path = f'{fq_dir}/{pe}_1.{ext}'
            pe2_path = f'{fq_dir}/{pe}_2.{ext}'
            command = f'STAR --runThreadN {star_nthreads} --genomeDir {star_idx_dir} --readFilesIn {pe1_path} {pe2_path} --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {star_out_dir}/{pe} --readFilesCommand zcat'
            pe_commands.append(command)
        commands = se_commands + pe_commands
        return commands

    @classmethod
    def merge_star_bam(cls, star_bam, star_out_dir):
        bams = glob.glob(f'{star_out_dir}/*.bam')
        bam_str = ' '.join(bams)
        command = f'samtools merge {star_bam} {bam_str}'
        return command
    
    @classmethod
    def assemble_rna_transcripts(cls, ref_anno, star_bam, strg_gtf, strg_nthreads):
        command = f'stringtie -p {strg_nthreads} -G {ref_anno} -o {strg_gtf} {star_bam}'
        return command

    @classmethod
    def gatk_add_read_group(cls, star_out_dir, gatk_out_dir, name):
        star_bam_files = glob.glob(
            f'{star_out_dir}/*.sortedByCoord.out.bam')
        star_bam_file_names = [os.path.basename(bam).split(
            '.sortedByCoord.out.')[0] for bam in star_bam_files]
        gatk_bam_rgs = [
            f'{gatk_out_dir}/{nm}_gatk_rg.bam' for nm in star_bam_file_names]
        commands = [f'gatk AddOrReplaceReadGroups -I {star_bam} -O {gatk_bam_rg} -LB lib1 -PL illumina -PU unit1 -SM {name}' for (
            star_bam, gatk_bam_rg) in zip(star_bam_files, gatk_bam_rgs)]
        return commands

    @classmethod
    def gatk_mark_duplicates(cls, gatk_out_dir):
        gatk_bam_rgs = glob.glob(f'{gatk_out_dir}/*_gatk_rg.bam')
        gatk_bam_file_names = [os.path.basename(bam).split('_gatk_rg.')[0] for bam in gatk_bam_rgs]
        gatk_bam_rg_mds = [f'{gatk_out_dir}/{nm}_gatk_rg_md.bam' for nm in gatk_bam_file_names]
        gatk_bam_rg_md_metrics = [f'{gatk_out_dir}/{nm}_gatk_rg_md_metrics.txt' for nm in gatk_bam_file_names]
        commands = [f'gatk MarkDuplicates -I {gatk_bam_rg} -O {gatk_bam_rg_md} -M {gatk_bam_rg_md_metric}' for (
            gatk_bam_rg, gatk_bam_rg_md, gatk_bam_rg_md_metric) in zip(gatk_bam_rgs, gatk_bam_rg_mds, gatk_bam_rg_md_metrics)]
        return commands

    @classmethod
    def gatk_split_ncigar_reads(cls, ref_fa, gatk_out_dir):
        gatk_bam_rg_mds = glob.glob(f'{gatk_out_dir}/*_gatk_rg_md.bam')
        gatk_bam_file_names = [os.path.basename(bam).split(
            '_gatk_rg_md.')[0] for bam in gatk_bam_rg_mds]
        gatk_bam_rg_md_sncs = [
            f'{gatk_out_dir}/{nm}_gatk_rg_md_snc.bam' for nm in gatk_bam_file_names]
        commands = [f'gatk SplitNCigarReads -R {ref_fa} -I {gatk_bam_rg_md} -O {gatk_bam_rg_md_snc}' for (
            gatk_bam_rg_md, gatk_bam_rg_md_snc) in zip(gatk_bam_rg_mds, gatk_bam_rg_md_sncs)]
        return commands

    @classmethod
    def gatk_haplotype_caller(cls, ref_fa, gatk_out_dir, gatk_haplotype_caller_nthreads):
        gatk_bam_rg_md_sncs = glob.glob(
            f'{gatk_out_dir}/*_gatk_rg_md_snc.bam')
        gatk_bam_file_names = [os.path.basename(bam).split('_gatk_rg_md_snc.')[
            0] for bam in gatk_bam_rg_md_sncs]
        gatk_unfiltered_vcfs = [
            f'{gatk_out_dir}/{nm}_gatk_unfiltered.vcf' for nm in gatk_bam_file_names]
        commands = [f'gatk HaplotypeCaller --native-pair-hmm-threads {gatk_haplotype_caller_nthreads} -R {ref_fa} -I {gatk_bam_rg_md_snc} -O {gatk_unfiltered_vcf}' for (
            gatk_bam_rg_md_snc, gatk_unfiltered_vcf) in zip(gatk_bam_rg_md_sncs, gatk_unfiltered_vcfs)]
        return commands

    @classmethod
    def gatk_filter_vcf(cls, gatk_out_dir):
        gatk_unfiltered_vcfs = glob.glob(
            f'{gatk_out_dir}/*_gatk_unfiltered.vcf')
        gatk_vcf_file_names = [os.path.basename(vcf).split('_gatk_unfiltered.')[
            0] for vcf in gatk_unfiltered_vcfs]
        gatk_vcf_prefices = [
            f'{gatk_out_dir}/{nm}_gatk' for nm in gatk_vcf_file_names]
        commands = [f'vcftools --vcf {gatk_unfiltered_vcf} --remove-indels --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out {gatk_vcf_prefix}' for (
            gatk_unfiltered_vcf, gatk_vcf_prefix) in zip(gatk_unfiltered_vcfs, gatk_vcf_prefices)]
        return commands

    @classmethod
    def gatk_merge_vcf(cls, gatk_out_dir, gatk_vcf):
        gatk_vcfs = glob.glob(f'{gatk_out_dir}/*_gatk.recode.vcf')
        gatk_vcfs_str = ' '.join([f'-I {vcf}' for vcf in gatk_vcfs])
        command = f'gatk MergeVcfs {gatk_vcfs_str} -O {gatk_vcf}'
        return command

    @classmethod
    def reset_ref_fa_with_gatk_vcf(cls, ref_fa, gatk_vcf, ref_fa_fixed):
        fa_ref = Fasta(ref_fa, sequence_always_upper=True, as_raw=True)
        vcf_gn = (line.strip().split('\t') for line in open(
            gatk_vcf) if not (line.startswith('#') or line.strip() == ''))
        vcf_dt = {}
        for lst in vcf_gn:
            if lst[0] in vcf_dt:
                vcf_dt[lst[0]].append(
                    (int(lst[1]), lst[3].upper(), lst[4].upper()))
            else:
                vcf_dt[lst[0]] = [
                    (int(lst[1]), lst[3].upper(), lst[4].upper())]
        chrs = list(fa_ref.keys())
        with open(ref_fa_fixed, 'w') as of:
            for i, chrom in enumerate(chrs):
                print(f'Processing {chrom} {i+1}/{len(chrs)}')
                ref_seq = list(fa_ref[chrom][:])
                assert len(ref_seq) == len(fa_ref[chrom])
                if chrom in vcf_dt:
                    vcf_lst = vcf_dt[chrom]
                    for vl in vcf_lst:
                        coord, ref_nt, alt_nt = vl
                        fa_nt = fa_ref[chrom][coord-1]
                        assert alt_nt in 'NACGT'
                        assert ref_nt == fa_nt
                        assert ref_nt != alt_nt
                        ref_seq[coord-1] = alt_nt
                seq = np.asarray(
                    ref_seq, dtype='|S1').tostring().decode('utf-8')
                assert len(seq) == len(fa_ref[chrom])
                of.write(f'>{chrom}\n')
                of.write(f'{seq}\n')
                of.flush()
        command = f'samtools faidx {ref_fa_fixed}'
        return command

    @classmethod
    def analyze_strg_gtf(cls, strg_gtf, ref_fa_fixed, strg_gtf_json, strg_gtf_fpkm_filter, name):
        def _analyze_attributes(lst):
            lt = lst[8].split(';')
            st = [t.strip().replace('"', '') for t in lt if t.strip() != '']
            dt = {}
            for t in st:
                k, v = t.split(' ')
                dt[k] = v
            if 'exon_number' in dt:
                dt['exon_number'] = int(dt['exon_number'])
            for k in ['cov', 'FPKM', 'TPM']:
                if k in dt:
                    dt[k] = float(dt[k])
            return dt
        gn_gtf = (line.strip().split('\t') for line in open(strg_gtf) if (
            (not line.startswith('#')) and 'ref_gene_id' in line))
        fa_ref = Fasta(ref_fa_fixed, sequence_always_upper=True)
        trans_dts = []
        idx = -1
        global_counter = 0
        trans_dt = None
        for lst in gn_gtf:
            adt = _analyze_attributes(lst)
            if lst[2] == 'transcript':
                idx += 1
                if idx > 0:
                    trans_dts.append(trans_dt)
                trans_dt = {}
                trans_dt['transcript_id'] = adt['transcript_id']
                trans_dt['chr'] = lst[0]
                trans_dt['start'] = int(lst[3])
                trans_dt['end'] = int(lst[4])
                trans_dt['strand'] = lst[6]
                trans_dt['gene_id'] = adt['gene_id']
                trans_dt['reference_id'] = adt['reference_id']
                trans_dt['ref_gene_id'] = adt['ref_gene_id']
                if 'ref_gene_name' in adt:
                    trans_dt['ref_gene_name'] = adt['ref_gene_name']
                trans_dt['FPKM'] = adt['FPKM']
                trans_dt['TPM'] = adt['TPM']
                trans_dt['cov'] = adt['cov']
                trans_dt['exon'] = []
            elif lst[2] == 'exon':
                assert adt['transcript_id'] == trans_dt['transcript_id'], lst
                exon_dt = {}
                global_counter += 1
                exon_dt['id'] = f'{name}.{global_counter}'
                exon_dt['exon_number'] = adt['exon_number']
                exon_dt['start'] = int(lst[3])
                exon_dt['end'] = int(lst[4])
                exon_dt['strand'] = lst[6]
                rc = False
                if lst[6] == '-':
                    rc = True
                seq = fa_ref.get_seq(
                    trans_dt['chr'], exon_dt['start'], exon_dt['end'], rc=rc)
                exon_dt['seq'] = str(seq)
                exon_dt['cov'] = adt['cov']
                trans_dt['exon'].append(exon_dt)
        else:
            trans_dts.append(trans_dt)

        global_counter = 0
        filtered_trans_dts = []
        for trans_dt in trans_dts:
            if trans_dt['FPKM'] < strg_gtf_fpkm_filter:
                continue
            for exon_dt in trans_dt['exon']:
                global_counter += 1
                exon_dt['id'] = f'{name}.{global_counter}'
            filtered_trans_dts.append(trans_dt)
        json.dump(filtered_trans_dts, open(
            strg_gtf_json, 'w', encoding='utf-8'), indent=3)

    @classmethod
    def analyze_viral_gtf(cls, gtf, viral_ref_fa, gtf_json, name):
        def _analyze_attributes(lst):
            lt = lst[8].split(';')
            st = [t.strip().replace('"', '') for t in lt if t.strip() != '']
            dt = {}
            for t in st:
                k, v = t.split(' ')
                dt[k] = v
            return dt
        gn_gtf = (line.strip().split('\t') for line in open(gtf) if (
            (not line.startswith('#')) and 'transcript_id' in line))
        fa_ref = Fasta(viral_ref_fa, sequence_always_upper=True)
        trans_dts = []
        global_counter = 0
        trans_dt = None
        for lst in gn_gtf:
            adt = _analyze_attributes(lst)
            if lst[2] == 'transcript':
                trans_dt = {}
                trans_dt['transcript_id'] = adt['transcript_id']
                trans_dt['chr'] = lst[0]
                trans_dt['start'] = int(lst[3])
                trans_dt['end'] = int(lst[4])
                trans_dt['strand'] = lst[6]
                trans_dt['gene_id'] = adt['gene_id']
                trans_dt['reference_id'] = adt['gene_id']
                trans_dt['ref_gene_id'] = adt['gene_id']
                if 'gene_name' in adt:
                    trans_dt['ref_gene_name'] = adt['gene_name']
                else:
                    if adt['gene_id'].startswith('gene'):
                        trans_dt['ref_gene_name'] = adt['gene_id']
                    else:
                        continue
                trans_dt['exon'] = []
                exon_dt = {}
                global_counter += 1
                exon_dt['id'] = f'{name}_viral.{global_counter}'
                exon_dt['start'] = int(lst[3])
                exon_dt['end'] = int(lst[4])
                exon_dt['strand'] = lst[6]
                rc = False
                if lst[6] == '-':
                    rc = True
                seq = fa_ref.get_seq(
                    trans_dt['chr'], exon_dt['start'], exon_dt['end'], rc=rc)
                exon_dt['seq'] = str(seq)
                trans_dt['exon'].append(exon_dt)
                trans_dts.append(trans_dt)
        json.dump(trans_dts, open(
            gtf_json, 'w', encoding='utf-8'), indent=3)

    @classmethod
    def analyze_viral_gtf_ori(cls, gtf, viral_ref_fa, gtf_json, name):
        def _analyze_attributes(lst):
            lt = lst[8].split(';')
            st = [t.strip().replace('"', '') for t in lt if t.strip() != '']
            dt = {}
            for t in st:
                k, v = t.split(' ')
                dt[k] = v
            return dt
        gn_gtf = (line.strip().split('\t') for line in open(gtf) if (
            (not line.startswith('#')) and 'transcript_id' in line))
        fa_ref = Fasta(viral_ref_fa, sequence_always_upper=True)
        trans_dts = []
        global_counter = 0
        trans_dt = None
        for lst in gn_gtf:
            adt = _analyze_attributes(lst)
            if lst[2] == 'transcript':
                trans_dt = {}
                trans_dt['transcript_id'] = adt['transcript_id']
                trans_dt['chr'] = lst[0]
                trans_dt['start'] = int(lst[3])
                trans_dt['end'] = int(lst[4])
                trans_dt['strand'] = lst[6]
                trans_dt['gene_id'] = adt['gene_id']
                trans_dt['reference_id'] = adt['gene_id']
                trans_dt['ref_gene_id'] = adt['gene_id']
                if 'gene_name' in adt:
                    trans_dt['ref_gene_name'] = adt['gene_name']
                else:
                    continue
                trans_dt['exon'] = []
                exon_dt = {}
                global_counter += 1
                exon_dt['id'] = f'{name}_viral.{global_counter}'
                exon_dt['start'] = int(lst[3])
                exon_dt['end'] = int(lst[4])
                exon_dt['strand'] = lst[6]
                rc = False
                if lst[6] == '-':
                    rc = True
                seq = fa_ref.get_seq(
                    trans_dt['chr'], exon_dt['start'], exon_dt['end'], rc=rc)
                exon_dt['seq'] = str(seq)
                trans_dt['exon'].append(exon_dt)
                trans_dts.append(trans_dt)
        json.dump(trans_dts, open(
            gtf_json, 'w', encoding='utf-8'), indent=3)

    @classmethod
    def build_background_rna_index(cls, strg_gtf_json, bg_rna_fa, bg_rna_idx_prefix, bowtie_nthreads):
        dts = json.load(open(strg_gtf_json, encoding='utf-8'))
        with open(bg_rna_fa, 'w') as of:
            for dt in dts:
                gene_name = dt['reference_id'] # transcript name, should be updated!!!
                chrom = dt['chr']
                if 'ref_gene_name' in dt:
                    gname = dt['ref_gene_name']
                    gene_name = f'{gene_name}|{gname}'
                for edt in dt['exon']:
                    exon_id = edt['id']
                    seq_name = f'>{exon_id}|{chrom}|{gene_name}'
                    seq = edt['seq']
                    of.write(f'{seq_name}\n{seq}\n')
        command = f'bowtie2-build {bg_rna_fa} {bg_rna_idx_prefix}'
        cls.proc_command(command)

    @classmethod
    def PFS_check(cls, pfs_pattern, pfs):
        if len(pfs_pattern) ==0:
            return True
        re_pfs = re.compile(''.join([f'[{cls.g_nucleotide_dict[p]}]' for p in pfs_pattern]))
        pfs_test = re_pfs.match(pfs)
        if pfs_test:
            return True
        else:
            return False
            

    @classmethod
    def generate_virus_grna(cls, strg_gtf_json, virus_gene_selected, bowtie_prefix, cas_dt, chr_name):
        dts = json.load(open(strg_gtf_json, encoding='utf-8'))
        dts = [dt for dt in dts if dt['chr'] == chr_name]
        assert len(dts) != 0, f'{chr_name} not found. Viral genome need to be checked!'
        for cas, l in bowtie_prefix:
            print(f'Generating guide RNA of {l}nt for {cas}...')
            pfs_pattern, pfs_loc = cas_dt[cas]
            if pfs_pattern is None:
                pfs_pattern = []
            pfs_len = len(pfs_pattern)
            with open(f'{bowtie_prefix[(cas, l)]}.fa', 'w') as of:
                for gene_name in virus_gene_selected:
                    exon_ids = []
                    seqs = []
                    for dt in dts:
                        if 'ref_gene_name' not in dt:
                            continue
                        if dt['ref_gene_name'] == gene_name:
                            for edt in dt['exon']:
                                seqs.append(edt['seq'])
                                exon_ids.append(edt['id'])
                    grna_lst = []
                    for eid, seq in zip(exon_ids, seqs):
                        num_grna = len(seq) - l - pfs_len + 1
                        if num_grna < 1:
                            continue
                        for i in range(num_grna):
                            grna_with_pfs = seq[i:i+l+pfs_len]
                            if pfs_loc == 'up':
                                pfs = grna_with_pfs[:pfs_len]
                                grna = grna_with_pfs[pfs_len:]
                            elif pfs_loc == 'down':
                                grna = grna_with_pfs[:-pfs_len]
                                pfs = grna_with_pfs[-pfs_len:]
                            else:
                                grna = grna_with_pfs
                                pfs = None
                            if cls.PFS_check(pfs_pattern, pfs):
                                grna_lst.append((eid, grna))
                                #grna_lst.append((eid, grna_with_pfs+'NN'))  ## test only
                   
                    grna_lst = list(filter(lambda g: 'TTT' not in g[1], grna_lst))
                    for i, (eid, grna) in enumerate(grna_lst, start=1):
                        seqname = f'>{eid}|{chr_name}|{gene_name}|{l}_{i}|{cas}'
                        of.write(f'{seqname}\n{grna}\n')

    @classmethod
    @dask.delayed
    def _analyze_sam(cls, i, cas, l, bowtie_prefix, cas_dt, bg_rna_fa):
        print(f'Processing #{i+1} of {len(bowtie_prefix)} sam files...')
        sam = f'{bowtie_prefix[(cas, l)]}.sam'
        json_file = f'{bowtie_prefix[(cas, l)]}.json'
        fa = Fasta(f'{bowtie_prefix[(cas, l)]}.fa',
                sequence_always_upper=True, as_raw=True)
        bg_fa = Fasta(bg_rna_fa, sequence_always_upper=True)
        sam_dt = {}
        sam_gn = (line.strip().split('\t') for line in open(sam) if line.strip() != '')
        pfs_pattern, pfs_loc = cas_dt[cas]
        if pfs_pattern is None:
            pfs_pattern = []
        for lst in sam_gn:
            gid = lst[0]
            others = lst[1:]
            if gid in sam_dt:
                sam_dt[gid].append(others)
            else:
                sam_dt[gid] = [others]
        js_lst = []
        for gid in sam_dt:
            gid_lst = gid.split('|')
            chrom = gid_lst[1]
            seq = fa[gid][:]
            gene = gid_lst[2]
            length = int(gid_lst[3].split('_')[0])
            gdt = {
                'gid': gid,
                'chr': chrom,
                'seq': seq,
                'length': length,
                'gene': gene,
                'profile': []
            }
            for lst in sam_dt[gid]:
                flag = int(lst[0])
                if flag not in (0, 16, 256, 272):
                    continue
                rid = lst[1]
                rchrom = rid.split('|')[1]
                rlen = len(bg_fa[rid])
                beg = int(lst[2])
                end = beg + length - 1
                rc = False
                strand = '+'
                if flag in (16, 272):
                    rc = True
                    strand = '-'
                    if cas != 'Cas12':
                        continue
                if len(pfs_pattern) == 0:
                    pfs = 'None'
                    # wbeg = beg #test only
                    # wend = end #test only
                else:
                    pfs_beg = 0
                    pfs_end = rlen + 1
                    if ((pfs_loc == 'up') and (strand == '+')) or ((pfs_loc == 'down') and (strand == '-')):
                        pfs_beg = beg - len(pfs_pattern)
                        pfs_end = beg - 1
                        # wbeg = pfs_beg #test only
                        # wend = end #test only
                    elif ((pfs_loc == 'up') and (strand == '-')) or ((pfs_loc == 'down') and (strand == '+')):
                        pfs_beg = end + 1
                        pfs_end = end + len(pfs_pattern)
                        # wbeg = beg #test only
                        # wend = pfs_end #test only
                    if (pfs_beg < 1) or (pfs_end > rlen):
                        continue
                    pfs = bg_fa.get_seq(rid, pfs_beg, pfs_end, rc=rc).seq
                
                if not cls.PFS_check(pfs_pattern, pfs):
                    continue
                # wseq = bg_fa.get_seq(rid, wbeg, wend, rc=rc).seq #test only
                mm = [int(l.split(':')[-1]) for l in lst if l.startswith('XM')][0]
                rseq = bg_fa.get_seq(rid, beg, end, rc=rc).seq
                refgene = rid.split('|')[-2]
                rgene = rid.split('|')[-1]
                target = False
                if rchrom == chrom:
                    if rseq == seq:
                        if rgene == gene:
                            if mm == 0:
                                target = True
                # for experiment
                forward_beg = max(1, beg - 100)
                forward_seq = bg_fa.get_seq(rid, forward_beg, end, rc=rc).seq
                backward_end = min(len(bg_fa[rid]), end+100)
                backward_seq = bg_fa.get_seq(rid, beg, backward_end, rc=rc).seq

                rdt = {
                    'rid': rid,
                    'chr': rchrom,
                    'beg': beg,
                    'end': end,
                    'strand': strand,
                    'seq': rseq,
                    'refgene': refgene,
                    'gene': rgene,
                    'mismatch': mm,
                    'pfs': pfs,
                    'target': target,
                    'forward':forward_seq,
                    'backward':backward_seq,
                    #'wseq': wseq
                }
                gdt['profile'].append(rdt)
            js_lst.append(gdt)
        json.dump(js_lst, open(json_file, 'w', encoding='utf-8'), indent=3)

    @classmethod
    def _analyze_sam_spk(cls, arg_lst):
        i, cas, l, bowtie_prefix, cas_dt, bg_rna_fa = arg_lst
        print(f'Processing #{i+1} of {len(bowtie_prefix)} sam files...')
        sam = f'{bowtie_prefix[(cas, l)]}.sam'
        json_file = f'{bowtie_prefix[(cas, l)]}.json'
        fa = Fasta(f'{bowtie_prefix[(cas, l)]}.fa',
                sequence_always_upper=True, as_raw=True)
        bg_fa = Fasta(bg_rna_fa, sequence_always_upper=True)
        sam_dt = {}
        sam_gn = (line.strip().split('\t') for line in open(sam) if line.strip() != '')
        pfs_pattern, pfs_loc = cas_dt[cas]
        if pfs_pattern is None:
            pfs_pattern = []
        for lst in sam_gn:
            gid = lst[0]
            others = lst[1:]
            if gid in sam_dt:
                sam_dt[gid].append(others)
            else:
                sam_dt[gid] = [others]
        js_lst = []
        for gid in sam_dt:
            gid_lst = gid.split('|')
            chrom = gid_lst[1]
            seq = fa[gid][:]
            gene = gid_lst[2]
            length = int(gid_lst[3].split('_')[0])
            gdt = {
                'gid': gid,
                'chr': chrom,
                'seq': seq,
                'length': length,
                'gene': gene,
                'profile': []
            }
            for lst in sam_dt[gid]:
                flag = int(lst[0])
                if flag not in (0, 16, 256, 272):
                    continue
                rid = lst[1]
                rchrom = rid.split('|')[1]
                rlen = len(bg_fa[rid])
                beg = int(lst[2])
                end = beg + length - 1
                rc = False
                strand = '+'
                if flag in (16, 272):
                    rc = True
                    strand = '-'
                    if cas != 'Cas12':
                        continue
                if len(pfs_pattern) == 0:
                    pfs = 'None'
                    # wbeg = beg #test only
                    # wend = end #test only
                else:
                    pfs_beg = 0
                    pfs_end = rlen + 1
                    if ((pfs_loc == 'up') and (strand == '+')) or ((pfs_loc == 'down') and (strand == '-')):
                        pfs_beg = beg - len(pfs_pattern)
                        pfs_end = beg - 1
                        # wbeg = pfs_beg #test only
                        # wend = end #test only
                    elif ((pfs_loc == 'up') and (strand == '-')) or ((pfs_loc == 'down') and (strand == '+')):
                        pfs_beg = end + 1
                        pfs_end = end + len(pfs_pattern)
                        # wbeg = beg #test only
                        # wend = pfs_end #test only
                    if (pfs_beg < 1) or (pfs_end > rlen):
                        continue
                    pfs = bg_fa.get_seq(rid, pfs_beg, pfs_end, rc=rc).seq
                
                if not cls.PFS_check(pfs_pattern, pfs):
                    continue
                # wseq = bg_fa.get_seq(rid, wbeg, wend, rc=rc).seq #test only
                mm = [int(l.split(':')[-1]) for l in lst if l.startswith('XM')][0]
                rseq = bg_fa.get_seq(rid, beg, end, rc=rc).seq
                refgene = rid.split('|')[-2]
                rgene = rid.split('|')[-1]
                target = False
                if rchrom == chrom:
                    if rseq == seq:
                        if rgene == gene:
                            if mm == 0:
                                target = True
                rdt = {
                    'rid': rid,
                    'chr': rchrom,
                    'beg': beg,
                    'end': end,
                    'strand': strand,
                    'seq': rseq,
                    'refgene': refgene,
                    'gene': rgene,
                    'mismatch': mm,
                    'pfs': pfs,
                    'target': target
                    #'wseq': wseq
                }
                gdt['profile'].append(rdt)
            js_lst.append(gdt)
        json.dump(js_lst, open(json_file, 'w', encoding='utf-8'), indent=3)

    @classmethod
    def mapping_virus_grna(cls, bowtie_prefix, bowtie_out_dir, bg_rna_idx_prefix, grna_mismatch, bowtie_nthreads,  bowtie_nprocesses):
        commands = []
        for cas, l in bowtie_prefix:
            penalty = -(1 + 6 * grna_mismatch)
            fa = f'{bowtie_prefix[(cas, l)]}.fa'
            sam = f'{bowtie_prefix[(cas, l)]}.sam'
            log = f'{bowtie_out_dir}/{os.path.basename(bowtie_prefix[(cas, l)])}_log.txt'
            command = f'bowtie2 --threads {bowtie_nthreads} --score-min C,{penalty},0 -a -f --no-hd -L 10 --gbar 31 -x {bg_rna_idx_prefix} {fa} -S {sam} 2>{log}'
            commands.append(command)
        cls.multiproc_command(bowtie_nprocesses, commands)

    @classmethod
    def gen_ava_dataframe_report(cls, grnas):
        chrs = []
        begs = []
        ends = []
        strands = []
        seqs = []
        refgenes = []
        genes = []
        mismatches = []
        pfs = []
        targets = []
        lengths = []
        for gdt in grnas:
            for rdt in gdt['profile']:
                lengths.append(gdt['length'])
                chrs.append(rdt['chr'])
                begs.append(rdt['beg'])
                ends.append(rdt['end'])
                strands.append(rdt['strand'])
                seqs.append(rdt['seq'])
                refgenes.append(rdt['refgene'])
                genes.append(rdt['gene'])
                mismatches.append(rdt['mismatch'])
                pfs.append(rdt['pfs'])
                targets.append(rdt['target'])

        grna_df_dt = {
            'chr': chrs,
            'beg': begs,
            'end': ends,
            'strand': strands,
            'seq': seqs,
            'refgene': refgenes,
            'gene': genes,
            'pfs': pfs,
            'length': lengths,
        }
        grna_df = pd.DataFrame(grna_df_dt)
        grna_df = grna_df.drop_duplicates().sort_values(by=['gene', 'length'])
        # grna_df_gdt = dict.fromkeys(grna_df['gene name'])
        # for g in grna_df_gdt:
        #     grna_df_gdt[g] = grna_df[grna_df['gene name'] == g]
        #     grna_df_gdt[g] = grna_df_gdt[g].drop_duplicates().sort_values(by=['length', 'beg'])
        return grna_df

    @classmethod
    def gen_ava_dataframe_report_det(cls, grnas):
        chrs = []
        seqs = []
        genes = []
        pfs = []
        lengths = []
        for gdt in grnas:
            chrs.append(gdt['chr'])
            seqs.append(gdt['seq'])
            genes.append(gdt['gene'])
            pfs.append(gdt['pfs'])
            lengths.append(gdt['length'])
        grna_df_dt = {
            'chr': chrs,
            'seq': seqs,
            'gene': genes,
            'pfs': pfs,
            'length': lengths
        }
        grna_df = pd.DataFrame(grna_df_dt)
        grna_df = grna_df.drop_duplicates().sort_values(by=['gene', 'length'])
        return grna_df


    @classmethod
    def gen_unava_dataframe_report(cls, grnas):
        gids = []
        gchrs = []
        gseqs = []
        ggenes = []
        chrs = []
        begs = []
        ends = []
        strands = []
        seqs = []
        forward_seqs = []
        backward_seqs = []
        refgenes = []
        genes = []
        mismatches = []
        pfs = []
        targets = []
        lengths = []
        for i, gdt in enumerate(grnas, start=1):
            chr_name = gdt['gid'].split('|')[1]
            cas = gdt['gid'].split('|')[-1]
            gid = f'{chr_name}_{i}_{cas}'
            for rdt in gdt['profile']:
                gids.append(gid)
                gchrs.append(gdt['chr'])
                gseqs.append(gdt['seq'])
                ggenes.append(gdt['gene'])
                lengths.append(gdt['length'])
                chrs.append(rdt['chr'])
                begs.append(rdt['beg'])
                ends.append(rdt['end'])
                strands.append(rdt['strand'])
                seqs.append(rdt['seq'])
                forward_seqs.append(rdt['forward'])
                backward_seqs.append(rdt['backward'])
                refgenes.append(rdt['refgene'])
                genes.append(rdt['gene'])
                mismatches.append(rdt['mismatch'])
                pfs.append(rdt['pfs'])
                targets.append(rdt['target'])

        grna_df_dt = {
            'gid': gids,
            'gchr': gchrs,
            'gseq': gseqs,
            'guide gene name': ggenes,
            'chr': chrs,
            'beg': begs,
            'end': ends,
            'strand': strands,
            'seq': seqs,
            'forward seq': forward_seqs,
            'backward seq': backward_seqs,
            'targeting refgene': refgenes,
            'targeting gene name': genes,
            '#mismatch': mismatches,
            'pfs': pfs,
            'length': lengths,
            'targeting status': targets
        }
        grna_df = pd.DataFrame(grna_df_dt)
        grna_df = grna_df.drop_duplicates().sort_values(by=['gid', 'guide gene name', 'length', 'chr', 'beg', '#mismatch'])
        return grna_df

    @classmethod
    def df2xls(cls, dfdt, xls_name):
        nwb = xls.Workbook()
        nwb.remove(nwb.active)
        for k in dfdt:
            ws = nwb.create_sheet(k)
            cols = list(dfdt[k].columns)
            ws.append(cols)
            for line in zip(*[dfdt[k][c] for c in cols]):
                ws.append(line)
        nwb.save(xls_name)



    @classmethod
    def elimination_report(cls, bowtie_prefix, cas_dt, rep_output_dir, name, bg_rna_fa):
        args = [(i, cas, l, bowtie_prefix, cas_dt, bg_rna_fa) for (i, (cas, l)) in enumerate(bowtie_prefix)]
        sam_analysis_results = []
        for arg in args:
            sam_analysis_results.append(cls._analyze_sam(*arg))
        dask.compute(sam_analysis_results)

        ava_grna_dt = {}
        unava_grna_dt = {}
        for cas, l in bowtie_prefix:
            if cas not in ava_grna_dt:
                ava_grna_dt[cas] = []
            if cas not in unava_grna_dt:
                unava_grna_dt[cas] = []
            json_file = f'{bowtie_prefix[(cas, l)]}.json'
            js_lst = json.load(open(json_file, encoding='utf-8'))
            available_grnas = []
            unavailable_grnas = []
            available_json = f'{bowtie_prefix[(cas, l)]}_available.json'
            unavailable_json = f'{bowtie_prefix[(cas, l)]}_unavailable.json'
            for gdt in js_lst:
                if len(gdt['profile']) == 1:
                    if gdt['profile'][0]['target'] == True:
                        available_grnas.append(gdt)
                    else:
                        unavailable_grnas.append(gdt)
                else:
                    unavailable_grnas.append(gdt)
            ava_grna_dt[cas].extend(available_grnas)
            unava_grna_dt[cas].extend(unavailable_grnas)
            json.dump(available_grnas, open(available_json, 'w', encoding='utf-8'), indent=3)
            json.dump(unavailable_grnas, open(unavailable_json, 'w', encoding='utf-8'), indent=3)

            stats_f = open(f'{bowtie_prefix[(cas, l)]}_stats.txt', 'w')
            num_ava = len(available_grnas)
            num_una = len(unavailable_grnas)
            num_total = num_ava + num_una
            stats_f.write(f'{num_total} guide RNA candidates in total; of these:\n')
            stats_f.write(f'\t{num_ava} ({round(num_ava / num_total * 100, 4)}%) guide RNA candidates available\n')
            stats_f.write(f'\t{num_una} ({round(num_una / num_total * 100, 4)}%) guide RNA candidates unavailable')
            stats_f.close()

        for i, cas in enumerate(ava_grna_dt):
            print(f'Processing #{i+1} of {len(ava_grna_dt)} Cas processed...')
            if len(ava_grna_dt[cas]) > 0:
                ava_grna_df = cls.gen_ava_dataframe_report(ava_grna_dt[cas])
                ava_grna_df.to_csv(f'{rep_output_dir}/{name}_{cas}_available.csv', index=None)
            
            if len(unava_grna_dt[cas]) > 0:
                unava_grna_df = cls.gen_unava_dataframe_report(unava_grna_dt[cas])
                unava_grna_df.to_csv(f'{rep_output_dir}/{name}_{cas}_unavailable.csv', index=None)

        # for cas in ava_grna_dt:
        #     if len(ava_grna_dt[cas]) > 0:
        #         ava_grna_df_gdt = cls.gen_ava_dataframe_report(ava_grna_dt[cas])
                
        #         cls.df2xls(ava_grna_df_gdt, f'{rep_output_dir}/{name}_{cas}_available.xlsx')
        #     if len(unava_grna_dt[cas]) > 0:
        #         unava_grna_df_gdt = cls.gen_unava_dataframe_report(unava_grna_dt[cas])
        #         cls.df2xls(unava_grna_df_gdt, f'{rep_output_dir}/{name}_{cas}_unavailable.xlsx')


    @classmethod
    def detection_report_all(cls, bowtie_prefix, cas_dt, rep_output_dir, name, bg_rna_fa):
        
        args = [(i, cas, l, bowtie_prefix, cas_dt, bg_rna_fa) for (i, (cas, l)) in enumerate(bowtie_prefix)]
        sam_analysis_results = []
        for arg in args:
            sam_analysis_results.append(cls._analyze_sam(*arg))
        dask.compute(sam_analysis_results)

        ava_grna_dt = {}
        unava_grna_dt = {}
        for cas, l in bowtie_prefix:
            if cas not in ava_grna_dt:
                ava_grna_dt[cas] = []
            if cas not in unava_grna_dt:
                unava_grna_dt[cas] = []
            json_file = f'{bowtie_prefix[(cas, l)]}.json'
            fa = Fasta(f'{bowtie_prefix[(cas, l)]}.fa', sequence_always_upper=True)
            available_json = f'{bowtie_prefix[(cas, l)]}_available.json'
            unavailable_json = f'{bowtie_prefix[(cas, l)]}_unavailable.json'
            available_grnas = []
            unavailable_grnas = []
            js_lst = json.load(open(json_file, encoding='utf-8'))
            jsdt = {gdt['gid']:gdt for gdt in js_lst}
            fadt = {c.name:[c[:].seq, True] for c in fa}
            gid_set = set(jsdt.keys())
            whole_set = set(fadt.keys())
            assert len(gid_set - whole_set) == 0
            for gid in jsdt:
                if len(jsdt[gid]['profile']) > 0:
                    fadt[gid][1] = False
            available_grnas = []
            unavailable_grnas = []
            for gid in fadt:
                if fadt[gid][1] == True:
                    gid_lst = gid.split('|')
                    chrom = gid_lst[1]
                    seq = fadt[gid][0]
                    gene = gid_lst[2]
                    length = int(gid_lst[3].split('_')[0])
                    gdt = {
                        'gid': gid,
                        'chr': chrom,
                        'seq': seq,
                        'gene': gene,
                        'pfs': cas_dt[cas][0],
                        'length': length
                    }
                    available_grnas.append(gdt)
                else:
                    unavailable_grnas.append(jsdt[gid])
            json.dump(available_grnas, open(available_json, 'w', encoding='utf-8'), indent=3)
            json.dump(unavailable_grnas, open(unavailable_json, 'w', encoding='utf-8'), indent=3)
            ava_grna_dt[cas].extend(available_grnas)
            unava_grna_dt[cas].extend(unavailable_grnas)

            stats_f = open(f'{bowtie_prefix[(cas, l)]}_stats.txt', 'w')
            num_ava = len(available_grnas)
            num_una = len(unavailable_grnas)
            num_total = num_ava + num_una
            stats_f.write(f'{num_total} guide RNA candidates in total; of these:\n')
            stats_f.write(f'\t{num_ava} ({round(num_ava / num_total * 100, 4)}%) guide RNA candidates available\n')
            stats_f.write(f'\t{num_una} ({round(num_una / num_total * 100, 4)}%) guide RNA candidates unavailable')
            stats_f.close()
        for i, cas in enumerate(ava_grna_dt):
            print(f'Processing #{i+1} of {len(ava_grna_dt)} Cas processed...')
            if len(ava_grna_dt[cas]) > 0:
                ava_grna_df = cls.gen_ava_dataframe_report_det(ava_grna_dt[cas])
                ava_grna_df.to_csv(f'{rep_output_dir}/{name}_{cas}_available.csv', index=None)
            if len(unava_grna_dt[cas]) > 0:
                unava_grna_df = cls.gen_unava_dataframe_report(unava_grna_dt[cas])
                unava_grna_df.to_csv(f'{rep_output_dir}/{name}_{cas}_unavailable.csv', index=None)
                

    @classmethod
    def detection_report(cls, bowtie_prefix, cas_dt, rep_output_dir, name, bg_rna_fa):
        
        args = [(i, cas, l, bowtie_prefix, cas_dt, bg_rna_fa) for (i, (cas, l)) in enumerate(bowtie_prefix)]
        sam_analysis_results = []
        for arg in args:
            sam_analysis_results.append(cls._analyze_sam(*arg))
        dask.compute(sam_analysis_results)

        ava_grna_dt = {}
        for cas, l in bowtie_prefix:
            if cas not in ava_grna_dt:
                ava_grna_dt[cas] = []
            json_file = f'{bowtie_prefix[(cas, l)]}.json'
            fa = Fasta(f'{bowtie_prefix[(cas, l)]}.fa', sequence_always_upper=True)
            available_json = f'{bowtie_prefix[(cas, l)]}_available.json'
            available_grnas = []
            js_lst = json.load(open(json_file, encoding='utf-8'))
            jsdt = {gdt['gid']:gdt for gdt in js_lst}
            fadt = {c.name:[c[:].seq, True] for c in fa}
            gid_set = set(jsdt.keys())
            whole_set = set(fadt.keys())
            assert len(gid_set - whole_set) == 0
            for gid in jsdt:
                if len(jsdt[gid]['profile']) > 0:
                    fadt[gid][1] = False
            available_grnas = []
            for gid in fadt:
                if fadt[gid][1] == True:
                    gid_lst = gid.split('|')
                    chrom = gid_lst[1]
                    seq = fadt[gid][0]
                    gene = gid_lst[2]
                    length = int(gid_lst[3].split('_')[0])
                    gdt = {
                        'gid': gid,
                        'chr': chrom,
                        'seq': seq,
                        'gene': gene,
                        'pfs': cas_dt[cas][0],
                        'length': length
                    }
                    available_grnas.append(gdt)
            json.dump(available_grnas, open(available_json, 'w', encoding='utf-8'), indent=3)
            ava_grna_dt[cas].extend(available_grnas)

            stats_f = open(f'{bowtie_prefix[(cas, l)]}_stats.txt', 'w')
            num_ava = len(available_grnas)
            stats_f.write(f'{num_ava} guide RNA candidates in total;\n')
            stats_f.close()
        for i, cas in enumerate(ava_grna_dt):
            print(f'Processing #{i+1} of {len(ava_grna_dt)} Cas processed...')
            if len(ava_grna_dt[cas]) > 0:
                ava_grna_df = cls.gen_ava_dataframe_report_det(ava_grna_dt[cas])
                ava_grna_df.to_csv(f'{rep_output_dir}/{name}_{cas}_available.csv', index=None)
                

            

    