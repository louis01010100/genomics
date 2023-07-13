import gzip
import random
import shutil
import string
from io import StringIO
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import pandas as pd
from Bio import bgzf

from .genomic_region import GenomicRegion
from .utils import is_gzipped, vcf2dict

COLUMN_IDX_MAP = {
    'CHROM': 0,
    'POS': 1,
    'ID': 2,
    'REF': 3,
    'ALT': 4,
    'QUAL': 5,
    'FILTER': 6,
    'INFO': 7,
}


class Vcf():

    def __init__(self, filepath, tmp_dir, n_threads=1):

        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        self.filepath = Path(filepath)
        self.n_threads = n_threads
        self.tmp_dir = Path(tmp_dir)

    @property
    def meta(self):
        if not getattr(self, '_meta', None):
            self._meta = _load_meta(self.filepath)
        return self._meta

    @property
    def header(self):
        if not getattr(self, '_header', None):
            self._header = _load_header(self.filepath)
        return self._header

    def delete(self):
        filepath = self.filepath
        log_filepath = self.tmp_dir / f'{filepath}.log'
        index_filepath = self.tmp_dir / f'{filepath}.csi'
        index_log_filepath = self.tmp_dir / f'{filepath}.csi.log'

        filepath.unlink(missing_ok=True)
        log_filepath.unlink(missing_ok=True)
        index_filepath.unlink(missing_ok=True)
        index_log_filepath.unlink(missing_ok=True)

    def rename(self, stem, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{stem}.vcf.bgz',
        )

        shutil.move(input_filepath, output_file)

        if delete_src:
            self.delete()

        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def bgzip(self, delete_src=False):

        if self.filepath.name.endswith('.vcf'):
            output_file = self.tmp_dir / self.filepath.name.replace(
                '.vcf', '.vcf.bgz')
            log_filepath = self.tmp_dir / f'{output_file.name}.log'
            cmd = (''
                   f' bgzip -c -@ {self.n_threads} {self.filepath}'
                   f' > {output_file}'
                   f' 2> {log_filepath}'
                   '')
            execute(cmd)

        elif self.filepath.name.endswith('.gz'):
            output_file = self.tmp_dir / self.filepath.name.replace(
                '.vcf', '').replace('.gz', '.vcf.bgz')
            log_filepath = self.tmp_dir / f'{output_file.name}.log'
            cmd = (''
                   f'gzip -dc {self.filepath}'
                   f' | bgzip -@ {self.n_threads} '
                   f' > {output_file}'
                   f' 2> {log_filepath}'
                   '')
            execute(cmd)

        elif self.filepath.name.endswith('.bgz'):
            output_file = self.tmp_dir / self.filepath.name.replace(
                '.vcf', '').replace('.bgz', '.vcf.bgz')

            if not output_file.samefile(self.filepath):
                shutil.copy2(self.filepath, output_file)

        else:
            raise Exception(self.filepath)

        if delete_src:
            self.delete()

        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def index(self):
        log_filepath = self.tmp_dir / f'{self.filepath.name}.csi.log'

        if (self.filepath.parent / f'{self.filepath.name}.csi').exists():
            return self

        cmd = (''
               f'bcftools index'
               f'      {self.filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return self

    def fix_header(self,
                   genome_index_filepath: Path,
                   delete_src: bool = False):

        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-header.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools reheader'
               f'    --fai {genome_index_filepath}'
               f'    -o {output_file}'
               f'    {input_filepath}'
               f'    &> {log_filepath}'
               '')
        execute(cmd)

        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def copy_to(self, output_file):
        output_file = Path(output_file)

        shutil.copy2(self.filepath, output_file)
        shutil.copy2(self.filepath.with_suffix('.bgz.csi'),
                     output_file.with_suffix('.bgz.csi'))

    def move_to(self, output_file):
        output_file = Path(output_file)

        shutil.move(self.filepath, output_file)
        index_file = self.filepath.with_suffix('.bgz.csi')
        if index_file.exists():
            shutil.move(index_file, output_file.with_suffix('.bgz.csi'))

    @property
    def samples(self) -> list[str]:

        cmd = (''
               f'bcftools query -l {self.filepath}'
               '')
        try:
            stdout, stderr = execute(cmd, pipe=True)
        except Exception as e:
            print(stderr)
            raise e

        return stdout.strip().split('\n')

    def view(self, chrom, start, stop=None):

        if not stop:
            stop = start
        cmd = (''
               f'bcftools view'
               f'      -H'
               f'      -r "{chrom}:{start}-{stop}"'
               f'      {self.filepath}'
               '')
        stdout, stderr = execute(cmd, pipe=True)
        return stdout.strip().split('\n')

    def drop_id(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-id.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x ID'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)

        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def include(self, criteria, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-ic.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -i {criteria}'
               f'      -O z'
               f'      -o {output_file}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def exclude(self, criteria, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-ex.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -e {criteria}'
               f'      -O z'
               f'      -o {output_file}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def filter(self, criteria: str = 'PASS', delete_src: bool = False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-filter.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -f "{criteria}"'
               f'      -O z'
               f'      -o {output_file}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def drop_format(self, fields, delete_src=False):
        expression = ','.join([f'FORMAT/{x}' for x in fields])
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-format.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x {expression}'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    @property
    def index_file(self):
        self.index()
        return self.filepath.with_suffix('.bgz.csi')

    def list_contig_names(self):
        self.index()
        cmd = (''
               f'bcftools index'
               f'       --stats {self.index_file}'
               '')

        stdout, stderr = execute(cmd, pipe=True)
        data = StringIO(stdout)

        df = pd.read_csv(data,
                         names=['contig_name', 'contig_size', 'n_records'],
                         sep='\t')
        return sorted(list(set(df['contig_name'])))

    def rename_chroms(self, chrom_map_file, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-rename_chroms.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      --rename-chrs {chrom_map_file}'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def drop_info(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-info.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x INFO'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def trim_alts(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-trim_alt.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   -a '
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def subset(self, expression, op_name, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{op_name}.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'   --include "{expression}"'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)

        if delete_src:
            self.delete()

        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def subset_biallelics(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-bi.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'   -M 2'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def subset_variants_by_id(self, ids: set, delete_src: bool = False):

        input_filepath = self.filepath

        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-variants.vcf.bgz',
        )

        with bgzf.open(input_filepath,
                       'rt') as ifd, bgzf.open(output_file, 'wt') as ofd:
            for line in ifd:
                if line.startswith('#'):
                    ofd.write(line)
                    continue

                id_ = line.split('\t', 4)[2]

                if id_ in ids:
                    ofd.write(line)

        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def subset_variants(self, coordinates_df: pd.DataFrame, delete_src=False):

        self.index()

        input_filepath = self.filepath

        coordinates_file = self.tmp_dir / ''.join(
            random.choices(string.ascii_letters, k=10))

        coordinates_df.to_csv(
            coordinates_file,
            header=False,
            index=False,
            sep='\t',
        )

        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-variants.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'      --regions-file {coordinates_file}'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    # def subset_chrom(self, chrom: str, delete_src: bool =False):
    #
    #     self.index()
    #     input_filepath = self.filepath
    #     output_file = self.tmp_dir / self.filepath.name.replace(
    #         '.vcf.bgz',
    #         f'-{chrom}.vcf.bgz',
    #     )
    #     log_filepath = self.tmp_dir / f'{output_file.name}.log'
    #
    #     cmd = (''
    #            f'bcftools view'
    #            f'      -r "{chrom}"'
    #            f'      -O z'
    #            f'      -o {output_file}'
    #            f'      --threads {self.n_threads}'
    #            f'      {input_filepath}'
    #            f'      &> {log_filepath}'
    #            '')
    #
    #     execute(cmd)
    #     if delete_src:
    #         self.delete()
    #     return Vcf(output_file, self.tmp_dir, self.n_threads)

    def include_chroms(self, chroms: set, delete_src: bool = False):
        input_filepath = self.filepath
        label = '_'.join(chroms)
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{label}.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        chroms = ','.join(chroms)

        cmd = (''
               f'bcftools view'
               f'   -t "{chroms}"'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def exclude_chroms(self, chroms: list, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-ex_chroms.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        ex_chroms = ','.join(chroms)

        cmd = (''
               f'bcftools view'
               f'   -t "^{ex_chroms}"'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def split_by_sample(self, delete_src=False):
        input_filepath = self.filepath
        output_dirpath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz', '-samples')

        output_dirpath.mkdir(exist_ok=True)
        log_filepath = output_dirpath / f'{input_filepath.name}.log'

        cmd = (''
               f'bcftools +split'
               f'   -O z'
               f'   -o {output_dirpath}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()

        bag = []

        for filepath in output_dirpath.glob('*.vcf.gz'):
            old_filepath = filepath
            filepath = old_filepath.with_suffix('.bgz')
            shutil.move(old_filepath, filepath)
            samplename = filepath.name.split('.')[0]

            tmp_dir = self.tmp_dir / samplename

            vcf = Vcf(filepath, tmp_dir)
            vcf.index()

            bag.append(vcf)

        return bag

    def subset_samples(self, samples: set, delete_src=False):
        self.index()
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-samples.vcf.bgz',
        )

        bag = []

        for sample in samples:
            bag.append({'sample_name': sample})

        samples_file = self.tmp_dir / 'samples.tsv'

        pd.DataFrame.from_records(bag).to_csv(
            samples_file,
            header=False,
            index=False,
            sep='\t',
        )

        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'   --samples-file {samples_file}'
               f'   --force-samples'
               f'   -O z'
               f'   -o {output_file}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def uppercase(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-uppercase.vcf.bgz',
        )
        with bgzf.open(input_filepath,
                       'rt') as ifd, bgzf.open(output_file, 'wt') as ofd:
            for line in ifd:
                if line.startswith('#'):
                    ofd.write(line)
                    continue
                tokens = line.split('\t')
                tokens[3] = tokens[3].upper()
                tokens[4] = tokens[4].upper()

                line = '\t'.join(tokens)

                ofd.write(line)

        if delete_src:
            self.delete()

        vcf = Vcf(output_file, self.tmp_dir, self.n_threads)

        return vcf

    def normalize(self,
                  genome_file: Path,
                  atomize=False,
                  split_multiallelics=False,
                  delete_src=False):

        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-norm.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools norm'
               f'      -f {genome_file}'
               f'      -c s'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               '')

        if atomize:
            cmd += '      -a'
        if split_multiallelics:
            cmd += '      -m -'

        cmd += f'      {input_filepath}'
        cmd += f'      &> {log_filepath}'

        execute(cmd)

        if delete_src:
            self.delete()

        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def fill_tags(self, tags=['AC', 'AN', 'AF', 'NS'], delete_src=False):

        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-tag.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        tags_str = ','.join(tags)

        cmd = (''
               f'bcftools +fill-tags'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      -- -t {tags_str}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def drop_gt(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-site.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'

        cmd = (''
               f'bcftools view'
               f'      -G'
               f'      -O z'
               f'      -o {output_file}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def sort(self, delete_src=False):
        input_filepath = self.filepath
        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-sort.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_file.name}.log'
        tmp_dir = self.tmp_dir / 'tmp_sort'

        cmd = (''
               f'bcftools sort'
               f'      -O z'
               f'      -o {output_file}'
               f'      --temp-dir {tmp_dir}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_file, self.tmp_dir, self.n_threads)

    @staticmethod
    def _fix_duplicates(input_filepath, output_file, columns):

        with gzip.open(input_filepath,
                       'rt') as ifd, output_file.open('wt') as ofd:

            ref_record = {
                'chrom': None,
                'pos': None,
                'ref': None,
                'alt_map': {},
            }

            for line in ifd:
                if line.startswith('#'):
                    ofd.write(line)
                    continue

                line = line.strip()

                items = line.split('\t', 5)
                chrom = items[0]
                pos = items[1]
                # id_ = items[2]
                ref = items[3]
                alt = items[4]

                if ref_record['chrom'] == chrom and \
                        ref_record['pos'] == pos and \
                        ref_record['ref'] == ref:
                    if alt in ref_record['alt_map']:
                        ref_line = ref_record['alt_map'][alt]
                        line = _new_vcf_record(line, ref_line, columns)
                    else:
                        ref_record['alt_map'][alt] = line
                        for a in alt.split(','):
                            ref_record['alt_map'][a] = line

                else:
                    alt_map = {alt: line}
                    for a in alt.split(','):
                        alt_map[a] = line

                    ref_record = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt_map': alt_map,
                    }
                ofd.write(line)
                ofd.write('\n')

    def annotate(self, annotations_vcf, columns: set, delete_src=False):
        self.index()
        input_filepath = self.filepath
        tmp1_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-annot_raw.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{tmp1_filepath.name}.log'

        columns_str = ','.join(columns)

        cmd = (''
               f'bcftools annotate'
               f'      --annotations {annotations_vcf.filepath}'
               f'      --columns {columns_str}'
               f'      -O z'
               f'      -o {tmp1_filepath}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)

        tmp2_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz', '-annot.vcf')

        Vcf._fix_duplicates(
            tmp1_filepath,
            tmp2_filepath,
            columns,
        )

        output_file = self.tmp_dir / tmp2_filepath.name.replace(
            '.vcf', '.vcf.bgz')

        log_filepath = self.tmp_dir / f'{tmp2_filepath.name}.log'
        cmd = (''
               f' bgzip -c -@ {self.n_threads} {tmp2_filepath}'
               f' > {output_file}'
               f' 2> {log_filepath}'
               '')

        execute(cmd)

        if delete_src:
            self.delete()
            tmp1_filepath.unlink()
            tmp2_filepath.unlink()

        return Vcf(output_file, self.tmp_dir, self.n_threads)

    def to_df(self,
              format_: str = None,
              site_only: bool = False,
              delete_src=False) -> pd.DataFrame:
        if not format_:
            if site_only:
                format_ = '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n'
            else:
                format_ = '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%SAMPLE\t%GT\t%TGT\n]'

        output_file = self.to_tsv(format_)

        if delete_src:
            self.delete()

        return pd.read_csv(output_file, header=0, sep='\t', dtype='str')

    def to_tsv(self, format_: str = None, site_only: bool = False) -> Path:
        if not format_:
            if site_only:
                format_ = '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\n'
            else:
                format_ = '[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%SAMPLE\t%GT\t%TGT\n]'

        input_filepath = self.filepath

        output_file = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '.tsv',
        )

        cnames_str = format_.strip('[]\n').replace(' ', '\t').replace('%', '')

        cnames = '\t'.join([x.lower() for x in cnames_str.split('\t')])

        format_ = repr(format_)[1:-1]    # turn \ into \\

        output_gz_filepath = output_file.with_suffix('.tsv.gz')

        if output_gz_filepath.exists():
            output_gz_filepath.unlink()

        cmd = (''
               f'echo "{cnames}" > {output_file};'
               f'bcftools query'
               f'      -f "{format_}"'
               f'      {input_filepath}'
               f'      >> {output_file};'
               f'gzip {output_file}'
               '')

        execute(cmd)

        return output_gz_filepath

    def __str__(self):
        return self.filepath.__str__()

    def __len__(self):
        input_filepath = Path(f'{self.filepath}.csi')
        if not input_filepath.exists():
            self.index()

        cmd = f'bcftools index -n {input_filepath}'

        stdout, stderr = execute(cmd, pipe=True)
        return int(stdout.strip())


def _fetch_seq(genome_file: Path, region: GenomicRegion) -> str:

    # samtools faidx hs38DH.fa 'chr12:1000000-1000010'

    chrom = region.chrom
    start = region.start
    stop = region.stop

    cmd = (''
           f'samtools faidx'
           f'     {genome_file}'
           f'     "{chrom}:{start}-{stop}"'
           '')

    stdout, stderr = execute(cmd, pipe=True)

    for line in stdout.split('\n'):
        if line.startswith('>'):
            continue
        return line.strip()


def _load_header(vcf):

    def fetch_header(fd):
        header = None
        for line in fd:
            line = line.strip()
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line
                break
            else:
                break
        return header

    if is_gzipped(vcf):
        with gzip.open(vcf, 'rt') as fd:
            return fetch_header(fd)
    else:
        with open(vcf, 'rt') as fd:
            return fetch_header(fd)


def _load_meta(vcf):

    def fetch_meta(fd):
        meta = []
        for line in fd:
            line = line.strip()
            if line.startswith('##'):
                meta.append(line)
                continue
            break
        return '\n'.join(meta)

    if is_gzipped(vcf):
        with gzip.open(vcf, 'rt') as fd:
            return fetch_meta(fd)
    else:
        with open(vcf, 'rt') as fd:
            return fetch_meta(fd)


def _new_vcf_record(current_line, ref_line, columns):
    current_record = current_line.split('\t')
    ref_record = ref_line.split('\t')

    for column in columns:
        if column in COLUMN_IDX_MAP:
            idx = COLUMN_IDX_MAP[column]
            current_record[idx] = ref_record[idx]

    tags = [x[5:] for x in columns if x.startswith('INFO/')]

    if len(tags) == 0:
        return '\t'.join(current_record)

    current_record_alt = current_record[COLUMN_IDX_MAP['ALT']]
    ref_record_alt = ref_record[COLUMN_IDX_MAP['ALT']]

    info_idx = COLUMN_IDX_MAP['INFO']

    if current_record_alt == ref_record_alt:
        current_record[info_idx] = ref_record[info_idx]
        return '\t'.join(current_record)

    current_info = current_record[COLUMN_IDX_MAP['INFO']]
    ref_info = ref_record[COLUMN_IDX_MAP['INFO']]

    current_record[info_idx] = _new_info(
        current_record_alt,
        current_info,
        ref_record_alt,
        ref_info,
        tags,
    )

    return '\t'.join(current_record)


def _get_prefix_suffix(new_ref, start, stop, ref, pos):

    allele_prefix = new_ref[0:pos - start]
    allele_suffix = new_ref[pos + len(ref) - start:]

    return allele_prefix, allele_suffix


def _new_info(
    current_alt,
    current_info,
    reference_alt,
    reference_info,
    tags,
):
    ralleles = reference_alt.split(',')
    calleles = current_alt.split(',')

    rtags = {}

    for tag_item in reference_info.split(';'):
        tag_name, tag_value_str = tag_item.split('=')
        tag_values = tag_value_str.split(',')
        if len(tag_values) > 1:
            rtags[tag_name] = {}
            for a, v in zip(ralleles, tag_values):
                rtags[tag_name][a] = v
        else:
            rtags[tag_name] = tag_values

    result = {}
    if current_info != '.':
        for tag_item in current_info.split(';'):
            tag_name, tag_value = tag_item.split('=')
            result[tag_name] = tag_value

    for tag_item in tags:
        if type(rtags[tag_name]) == str:

            result[tag_name] = rtags[tag_name]

        elif type(rtags[tag_name]) == dict:
            bag = []

            allele2value = rtags[tag_name]

            for allele in calleles:
                if allele in allele2value:
                    bag.append(allele2value[allele])
                else:
                    bag.append('.')
            result[tag_name] = ','.join(bag)
    sorted_keys = sorted(result.keys())

    bag = []

    for k in sorted_keys:
        bag.append(f'{k}={result[k]}')
    return ';'.join(bag)


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

            proc.wait()

            if proc.returncode:
                raise Exception(cmd)


def sync_alleles(vcf_file_x: Path, vcf_file_y: Path, output_dir) -> tuple:
    output_dir.mkdir(parents=True, exist_ok=True)

    vcf_file_x = Vcf(vcf_file_x, output_dir).bgzip().filepath
    vcf_file_y = Vcf(vcf_file_y, output_dir).bgzip().filepath

    modification = vcf2dict(vcf_file_x, vcf_file_y)

    vcf_file_x_modified = output_dir / vcf_file_x.name.replace(
        '.vcf.bgz', '-sync.vcf')

    vcf_file_y_modified = output_dir / vcf_file_y.name.replace(
        '.vcf.bgz', '-sync.vcf')

    _update_vcf_file(vcf_file_x, modification, vcf_file_x_modified)
    _update_vcf_file(vcf_file_y, modification, vcf_file_y_modified)

    vcf_file_x = Vcf(vcf_file_x_modified, output_dir).bgzip().index().filepath
    vcf_file_y = Vcf(vcf_file_y_modified, output_dir).bgzip().index().filepath

    return vcf_file_x, vcf_file_y


def _update_vcf_file(input_vcf_file, modification, output_vcf_file):

    with gzip.open(input_vcf_file,
                   'rt') as ifh, output_vcf_file.open('wt') as ofh:
        for line in ifh:
            if line.startswith('#'):
                ofh.write(line)
                continue

            chrom, pos, id_, ref, alt, rest = line.strip().split('\t', 5)
            if (chrom, pos) not in modification:
                ofh.write(line)
                continue

            updated_allele_pair = modification[(chrom,
                                                pos)].get_updated_allele_pair(
                                                    ref, alt)
            new_ref = updated_allele_pair['ref']
            new_alt = updated_allele_pair['alt']

            ofh.write('\t'.join([chrom, pos, id_, new_ref, new_alt, rest]))
            ofh.write('\n')


def merge(vcf_files: list,
          output_file: Path,
          tmp_dir: Path,
          flag: str = 'none',
          n_threads: int = 1) -> Vcf:
    tmp_dir.mkdir(parents=True, exist_ok=True)

    vcfs_file = tmp_dir / ''.join(random.choices(string.ascii_letters, k=10))

    with vcfs_file.open('wt') as fd:
        for vcf_file in vcf_files:
            Vcf(vcf_file, tmp_dir).bgzip().index()
            fd.write(str(vcf_file) + '\n')

    output_file = Path(output_file)

    tmp_filename = output_file.name.replace(
        '.vcf.bgz',
        '-merge.vcf.bgz',
    )

    tmp_dir = Path(tmp_dir)
    log_filepath = tmp_dir / f'{tmp_filename}.log'
    tmp_filepath = tmp_dir / f'{tmp_filename}'

    cmd = (''
           f'bcftools merge'
           f'      --merge {flag}'
           f'      --file-list {vcfs_file}'
           f'      -O z'
           f'      -o {tmp_filepath}'
           f'      --threads {n_threads}'
           f'      &> {log_filepath}'
           '')

    execute(cmd)
    result = Vcf(tmp_filepath, tmp_dir).sort(delete_src=True).index()

    result.move_to(output_file)

    vcfs_file.unlink()

    return Vcf(output_file, tmp_dir)


def concat(vcf_files: list,
           output_file: Path,
           tmp_dir: Path,
           n_threads: int = 1) -> Vcf:
    tmp_dir.mkdir(parents=True, exist_ok=True)

    vcfs_file = tmp_dir / 'vcfs.tsv'

    with vcfs_file.open('wt') as fd:
        for vcf_file in vcf_files:
            Vcf(vcf_file, tmp_dir).bgzip().index()
            print(vcf_file, file=fd)

    output_file = Path(output_file)

    tmp_filename = output_file.name.replace(
        '.vcf.bgz',
        '-concat.vcf.bgz',
    )

    tmp_dir = Path(tmp_dir)
    log_filepath = tmp_dir / f'{tmp_filename}.log'
    tmp_filepath = tmp_dir / f'{tmp_filename}'

    cmd = (''
           f'bcftools concat'
           f'      --allow-overlaps'
           f'      --file-list {vcfs_file}'
           f'      -O z'
           f'      -o {tmp_filepath}'
           f'      --threads {n_threads}'
           f'      &> {log_filepath}'
           '')

    execute(cmd, debug=True)
    result = Vcf(tmp_filepath, tmp_dir).sort(delete_src=True).index()

    result.copy_to(output_file)

    return Vcf(output_file, tmp_dir)
