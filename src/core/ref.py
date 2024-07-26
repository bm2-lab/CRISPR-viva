from .processor import Processor
import os

class Ref:
    def __init__(self, ori_ref_fa, virus_fa, ori_ref_anno, virus_anno, output_dir, staridx_flag, star_idx_nthreads, star_sjdb_overhang):
        self.ref_fa = ori_ref_fa
        self.virus_fa = virus_fa
        self.ref_anno = ori_ref_anno
        self.virus_anno = virus_anno
        self.output_dir = output_dir
        self.staridx_flag = staridx_flag
        self.star_idx_nthreads = star_idx_nthreads
        self.star_sjdb_overhang = star_sjdb_overhang
        os.makedirs(self.output_dir, exist_ok=True)
    
    def run(self):
        if (self.virus_fa is not None) and (self.virus_anno is not None):
            self.ref_fa = Processor.concat_fasta(self.ref_fa, self.virus_fa, self.output_dir)
            self.ref_anno = Processor.concat_gtf(self.ref_anno, self.virus_anno, self.output_dir)
        else:
            Processor.prepare_fasta(self.ref_fa)
            Processor.prepare_gtf(self.ref_anno)
        Processor.check_ref_consistency(self.ref_fa, self.ref_anno)

        if self.staridx_flag == 1:
            ref_fa_name = os.path.splitext(os.path.basename(self.ref_fa))[0]
            star_idx_dir = f'{self.output_dir}/star_idx/{ref_fa_name}'
            os.makedirs(star_idx_dir, exist_ok=True)
            command = Processor.build_star_index(self.ref_fa, self.ref_anno, star_idx_dir, self.star_idx_nthreads, self.star_sjdb_overhang)
            Processor.proc_command(command)

        




    