from pathlib import Path
from subprocess import Popen, PIPE, STDOUT


class VCF():
    def __init__(self, filepath, tmp_dir):
        self.filepath = filepath
        self.tmp_dir = Path(tmp_dir)
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

    def bgzip(self, n_threads=1):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.gz',
            '.vcf.bgz',
        )
        cmd = f'gzip -dc {input_filepath} | bgzip -@ {n_threads} > {output_filepath}'

        execute(cmd)

        return Vcf(output_filepath, self.tmp_dir)

    def subset_samples(self, sample_file):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-sample.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'   --samples-file {sample_file}'
               f'   --force-samples '
               f'   -O z'
               f'   -o {output_filepath}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def fix_header(self, genome_index_filepath):

        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-header.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools reheader'
               f'    --fai {genome_index_filepath}'
               f'    -o {output_filepath}'
               f'    {input_filepath}'
               f'    &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def normalize(self, genome_filepath, atomize=False):

        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-norm.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        if atomize:

            cmd = (''
                   f'bcftools norm'
                   f'      -a'
                   f'      -f {genome_filepath}'
                   f'      -c w'
                   f'      -O z'
                   f'      -o {output_filepath}'
                   f'      {input_filepath}'
                   f'      &> {log_filepath}'
                   '')
        else:
            cmd = (''
                   f'bcftools norm'
                   f'      -f {genome_filepath}'
                   f'      -c w'
                   f'      -O z'
                   f'      -o {output_filepath}'
                   f'      {input_filepath}'
                   f'      &> {log_filepath}'
                   '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def drop_info(self):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-info.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x INFO'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def fill_af(self):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-af.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools +fill-tags'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      -- -t AC,AN,AF,NS'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def drop_gt(self):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-site.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -G'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir)

    def index(self):
        input_filepath = self.filepath
        log_filepath = self.tmp_dir / f'{input_filepath.name}-index.log'

        cmd = (''
               f'bcftools index {input_filepath'
               f'      &> {log_filepath}'
               '')

        execute(cmd)


def execute(cmd):
    with Popen(cmd, shell=True, text=True, stdout=PIPE,
               stderr=STDOUT) as proc:
        for line in proc.stdout:
            print(line, end='')

        if proc.returnstatus:
            raise Exception(cmd)
