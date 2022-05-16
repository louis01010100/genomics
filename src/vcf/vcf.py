from pathlib import Path
import pandas as pd
from subprocess import Popen, PIPE, STDOUT
import shutil
import gzip


class Vcf():
    def __init__(self, filepath, tmp_dir, n_threads=1):

        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        self.filepath = filepath
        self.n_threads = n_threads
        self.tmp_dir = tmp_dir

    def to_df(self):
        with gzip.open(self.filepath, 'rt') as fd:
            for line in fd:
                if line.startswith('##'):
                    continue
                if line.startswith('#'):
                    names = line.strip()[1:].split('\t')
                    break
                assert False

            df = pd.read_csv(fd, names=names, sep='\t', dtype='str')
        return df

    def bgzip(self):

        if self.filepath.name.endswith('.vcf'):
            output_filepath = self.tmp_dir / self.filepath.name.replace(
                '.vcf', '.vcf.bgz')
            log_filepath = self.tmp_dir / f'{output_filepath.name}.log'
            cmd = (''
                   f' bgzip -c -@ {self.n_threads} {self.filepath}'
                   f' > {output_filepath}'
                   f' 2> {log_filepath}'
                   '')
            execute(cmd)

        elif self.filepath.name.endswith('.vcf.gz'):
            output_filepath = self.tmp_dir / self.filepath.name.replace(
                '.vcf.gz', '.vcf.bgz')
            log_filepath = self.tmp_dir / f'{output_filepath.name}.log'
            cmd = (''
                   f'gzip -dc {self.filepath}'
                   f' | bgzip -@ {self.n_threads} '
                   f' > {output_filepath}'
                   f' 2> {log_filepath}'
                   '')
            execute(cmd)
        elif self.filepath.name.endswith('.vcf.bgz'):
            output_filepath = self.tmp_dir / self.filepath.name

            if not output_filepath.exists():
                shutil.copy2(self.filepath, output_filepath)
        else:
            raise Exception(self.filepath)

        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def index(self):
        log_filepath = self.tmp_dir / f'{self.filepath.name}.csi.log'

        cmd = (''
               f'bcftools index'
               f'      {self.filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)

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
               f'      --threads {self.n_threads}'
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
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def uppercase(self):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-uppercase.vcf',
        )
        with gzip(input_filepath,
                  'rt') as ifd, output_filepath.open('wt') as ofd:
            for line in ifd:
                line = line.strip()
                if line.startswith('#'):
                    print(line, file=ofd)
                    continue
                tokens = line.split('\t')
                tokens[3] = tokens[3].upper()
                tokens[4] = tokens[4].upper()

                line = '\t'.join(tokens)

                print(line, file=ofd)

        self.bgzip()

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
                   f'      --threads {self.n_threads}'
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
                   f'      --threads {self.n_threads}'
                   f'      {input_filepath}'
                   f'      &> {log_filepath}'
                   '')

        execute(cmd)
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def fill_tags(self, tags=['AC', 'AN', 'AF', 'NS']):

        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-tag.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        tags_str = ','.join(tags)

        cmd = (''
               f'bcftools +fill-tags'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      -- -t {tags_str}'
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
               f'      --threads {self.n_threads}'
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

    def subset_variants(self, chrom, start, stop):

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
               f'      --threads {self.n_threads}'
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
               f'      --threads {self.n_threads}'
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
        input_filepath = f'{self.filepath}.csi'

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
           f'      --threads {n_threads}'
           f'      &> {log_filepath}'
           '')

    execute(cmd)
    return Vcf(output_filepath, tmp_dir)


def execute(cmd, pipe=False, debug=False):
    if debug:
        print(cmd)
    if pipe:
        with Popen(cmd, shell=True, text=True, stdout=PIPE,
                   stderr=PIPE) as proc:
            stdout_data, stderr_data = proc.communicate()

            return stdout_data, stderr_data
    else:

        with Popen(cmd, shell=True, text=True, stdout=PIPE,
                   stderr=STDOUT) as proc:
            for line in proc.stdout:
                print(line, end='')

            if proc.returncode:
                raise Exception(cmd)
