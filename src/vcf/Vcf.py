from pathlib import Path
from subprocess import Popen, PIPE, STDOUT


class Vcf():
    def __init__(self, filepath, workspace_dir):
        self.filepath = filepath
        self.workspace_dir = Path(workspace_dir)
        self.workspace_dir.mkdir(parents=True, exist_ok=True)

    def bgzip(self, n_threads=1):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.gz',
            '.vcf.bgz',
        )
        cmd = (''
               f'gzip -dc {input_filepath}'
               f' | bgzip -@ {n_threads} '
               f' > {output_filepath}'
               '')

        execute(cmd)

        return Vcf(output_filepath, self.workspace_dir)

    def fix_header(self, genome_index_filepath):

        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-header.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools reheader'
               f'    --fai {genome_index_filepath}'
               f'    -o {output_filepath}'
               f'    {input_filepath}'
               f'    &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def drop_info(self):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-info.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x INFO'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def subset_samples(self, sample_file):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-sample.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

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
        return Vcf(output_filepath, self.workspace_dir)

    def normalize(self, genome_filepath, atomize=False):

        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-norm.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

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
        return Vcf(output_filepath, self.workspace_dir)

    def fill_af(self):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-af.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools +fill-tags'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      -- -t AC,AN,AF,NS'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def drop_gt(self):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-site.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -G'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def sort(self):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-sort.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'
        tmp_dir = self.workspace_dir / 'tmp_sort'

        cmd = (''
               f'bcftools sort'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      --temp-dir {tmp_dir}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def index(self):
        input_filepath = self.filepath
        log_filepath = self.workspace_dir / f'{input_filepath.name}.csi.log'

        cmd = (''
               f'bcftools index {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)

    def subset_variants(self, chrom, start, stop):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{chrom}_{start}_{stop}.vcf.bgz',
        )
        log_filepath = self.workspace_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -r "{chrom}:{start}-{stop}"'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.workspace_dir)

    def to_tsv(self, fields):
        input_filepath = self.filepath
        output_filepath = self.workspace_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '.tsv',
        )

        header = '\t'.join(fields)
        header_pattern = '\t'.join([f'%{x}' for x in fields]) + '\n'
        # f'      -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\t%INFO/AC\t%INFO/NS\n"'

        cmd = (''
               f'echo "{header}" > {output_filepath};\n'
               f'bcftools query'
               f'      -f "{header_pattern}"'
               f'      {input_filepath}'
               f'      >> {output_filepath};'
               f'gzip {output_filepath}'
               '')

        execute(cmd, debug=True)
        return output_filepath

    @staticmethod
    def concat(vcf_list, output_filepath, tmp_dir):
        output_filepath = Path(output_filepath)
        tmp_dir = Path(tmp_dir)
        log_filepath = tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools concat'
               f'      --allow-overlaps'
               f'      --file-list {vcf_list}'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, tmp_dir)


def execute(cmd, debug=False):
    if debug:
        print(cmd)

    with Popen(cmd, shell=True, text=True, stdout=PIPE,
               stderr=STDOUT) as proc:
        for line in proc.stdout:
            print(line, end='')

        if proc.returncode:
            raise Exception(cmd)
