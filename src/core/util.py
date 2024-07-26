from .processor import Processor
import json
import os
import subprocess
from pyfaidx import Fasta

def proc_sra(sra_dir, sra_nprocesses, fq_dir):
    commands = Processor.dump_fq(fq_dir, sra_dir)
    Processor.multiproc_command(sra_nprocesses, commands)


def proc_gff(gff_path):
    Processor.prepare_gtf(gff_path)


def proc_csv(name, output_dir, cas_definition_json, bg_rna_fa, bl_elim):
    if cas_definition_json is None:
        cas_definition = Processor.g_default_cas_definition
    else:
        cas_definition = json.load(open(cas_definition_json, encoding='utf-8'))
    os.makedirs(output_dir, exist_ok=True)
    bowtie_out_dir = f'{output_dir}/bowtie_out'
    os.makedirs(bowtie_out_dir, exist_ok=True)
    cas_dt = {cas:(cas_definition[cas]['pfs'], cas_definition[cas]['loc']) for cas in cas_definition}
    bowtie_prefix = {}
    for cas in cas_definition:
        for l in range(cas_definition[cas]['min_len'], cas_definition[cas]['max_len']+1):
            subprocess.call(f'mkdir -p {bowtie_out_dir}/{cas}/{name}_{cas}_{l}', executable='/bin/bash', shell=True)
            #os.makedirs(f'{bowtie_out_dir}/{cas}/{name}_{cas}_{l}', exist_ok=True)
            bowtie_prefix[(cas, l)] = f'{bowtie_out_dir}/{cas}/{name}_{cas}_{l}/{name}_{cas}_{l}'
            assert os.path.exists(f'{bowtie_prefix[(cas, l)]}.fa') is True, f'{bowtie_prefix[(cas, l)]}.fa'
            _ = Fasta(f'{bowtie_prefix[(cas, l)]}.fa')
            #print(os.path.exists(f'{bowtie_prefix[(cas, l)]}.fa'))
            #subprocess.call(f'samtools faidx {bowtie_prefix[(cas, l)]}.fa')
    
    if bl_elim:
        Processor.elimination_report(bowtie_prefix, cas_dt, output_dir, name, bg_rna_fa)
    else:
        Processor.detection_report(bowtie_prefix, cas_dt, output_dir, name, bg_rna_fa)
    

def proc_build_background_rna_index(strg_gtf, strg_gtf_fpkm_filter, ref_fa_fixed, output_dir, bowtie_nthreads, suffix):
    os.makedirs(output_dir, exist_ok=True)
    bowtie_out_dir = f'{output_dir}/bowtie_out_{suffix}'
    os.makedirs(bowtie_out_dir, exist_ok=True)
    strg_gtf_json = f'{bowtie_out_dir}/{os.path.splitext(os.path.basename(strg_gtf))[0]}.json'
    name = f'{os.path.splitext(os.path.basename(ref_fa_fixed))[0]}'
    bg_rna_fa = f'{bowtie_out_dir}/{name}_bg_rna.fa'
    
    bg_rna_idx_dir = f'{bowtie_out_dir}/bowtie2_idx'
    os.makedirs(bg_rna_idx_dir, exist_ok=True)
    bg_rna_idx_prefix = f'{bg_rna_idx_dir}/{name}_bg_rna'

    Processor.analyze_strg_gtf(strg_gtf, ref_fa_fixed, strg_gtf_json, strg_gtf_fpkm_filter, name)

    Processor.build_background_rna_index(strg_gtf_json, bg_rna_fa, bg_rna_idx_prefix, bowtie_nthreads)