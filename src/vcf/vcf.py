from pathlib import Path
from subprocess import Popen, PIPE, STDOUT
import shutil
import gzip


class Vcf():
    def __init__(self, filepath, tmp_dir, n_threads=1):

        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        self.filepath = bgzip(filepath, tmp_dir, n_threads)
        index(self.filepath, tmp_dir)
        self.n_threads = n_threads
        self.tmp_dir = tmp_dir

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
               f'    -@ {self.n_threads}'
               f'    {input_filepath}'
               f'    &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

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
               f'      -@ {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

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
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

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
                   f'      -@ {self.n_threads}'
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
                   f'      -@ {self.n_threads}'
                   f'      {input_filepath}'
                   f'      &> {log_filepath}'
                   '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

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
               f'      -@ {self.n_threads}'
               f'      {input_filepath}'
               f'      -- -t AC,AN,AF,NS'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

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
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def sort(self):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-sort.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'
        tmp_dir = self.tmp_dir / 'tmp_sort'

        cmd = (''
               f'bcftools sort'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      --temp-dir {tmp_dir}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def subset_variants(self, start, stop, chrom=None):
        def get_chrom(self):
            with gzip.open(self.filepath, 'rt') as fd:
                for line in fd:
                    if line.startswith('#'):
                        continue
                    return line.split('\t', maxsplit=1)[0]

        if not chrom:
            chrom = get_chrom(self)

        print(chrom)

        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{chrom}_{start}_{stop}.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -r "{chrom}:{start}-{stop}"'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def annotate(self, annotations_vcf, columns):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-annot.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      --annotations {annotations_vcf.filepath}'
               f'      --columns {columns}'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      -@ {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def to_tsv(self, query_str):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '.tsv',
        )

        column_header = query_str.strip('[]\n').replace(' ', '\t').replace(
            '%', '')

        output_gz_filepath = output_filepath.with_suffix('.tsv.gz')

        if output_gz_filepath.exists():
            output_gz_filepath.unlink()

        cmd = (''
               f'echo "{column_header}" > {output_filepath};'
               f'bcftools query'
               f'      -f "{query_str}"'
               f'      {input_filepath}'
               f'      >> {output_filepath};'
               f'gzip {output_filepath}'
               '')

        execute(cmd)
        return output_gz_filepath

    def __str__(self):
        return self.filepath.__str__()

    def __len__(self):
        input_filepath = self.filepath

        cmd = f'bcftools index -n {input_filepath}'

        stdout, stderr = execute(cmd, pipe=True)
        return int(stdout.strip())


def concat(vcfs, output_filepath, tmp_dir, n_threads=1):

    vcf_list_file = tmp_dir / 'vcfs.tsv'

    with vcf_list_file.open('wt') as fd:
        for vcf in vcfs:
            print(vcf, file=fd)

    output_filepath = Path(output_filepath)
    tmp_dir = Path(tmp_dir)
    log_filepath = tmp_dir / f'{output_filepath.name}.log'

    cmd = (''
           f'bcftools concat'
           f'      --allow-overlaps'
           f'      --file-list {vcf_list_file}'
           f'      -O z'
           f'      -o {output_filepath}'
           f'      -@ {n_threads}'
           f'      &> {log_filepath}'
           '')

    execute(cmd)
    return Vcf(output_filepath, tmp_dir)


def bgzip(input_filepath, tmp_dir, n_threads=1):
    if input_filepath.name.endswith('.vcf'):
        output_filepath = tmp_dir / input_filepath.name.replace(
            '.vcf', '.vcf.bgz')
        cmd = (''
               f' bgzip -c -@ {n_threads} {input_filepath}'
               f' > {output_filepath}'
               '')
        execute(cmd)

    elif input_filepath.name.endswith('.vcf.gz'):
        output_filepath = tmp_dir / input_filepath.name.replace(
            '.vcf.gz', '.vcf.bgz')
        cmd = (''
               f'gzip -dc {input_filepath}'
               f' | bgzip -@ {n_threads} '
               f' > {output_filepath}'
               '')
        execute(cmd)
    elif input_filepath.name.endswith('.vcf.bgz'):
        output_filepath = tmp_dir / input_filepath.name
        shutil.copy2(input_filepath, output_filepath)
    else:
        raise Exception(input_filepath)
    return output_filepath


def index(input_filepath, tmp_dir):
    log_filepath = tmp_dir / f'{input_filepath.name}.csi.log'

    cmd = (''
           f'bcftools index'
           f'      {input_filepath}'
           f'      &> {log_filepath}'
           '')

    execute(cmd)


def execute(cmd, pipe=False, debug=False):
    if debug:
        print(cmd)

    with Popen(cmd, shell=True, text=True, stdout=PIPE,
               stderr=STDOUT) as proc:
        if pipe:
            stdout_data, stderr_data = proc.communicate()
        else:
            for line in proc.stdout:
                print(line, end='')

        if proc.returncode:
            raise Exception(cmd)

        if pipe:
            return stdout_data, stderr_data