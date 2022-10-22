from pathlib import Path
import pandas as pd
from subprocess import Popen, PIPE, STDOUT
import shutil
import gzip

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

    def to_df(self):

        with gzip.open(self.filepath, 'rt') as fd:
            for line in fd:

                if line.startswith('##'):
                    continue
                if line.startswith('#'):
                    cnames = line.strip()[1:].split('\t')

                    break
                assert False

            df = pd.read_csv(fd, names=cnames, sep='\t', dtype='str')
        return df

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
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{stem}.vcf.bgz',
        )

        shutil.move(input_filepath, output_filepath)

        if delete_src:
            self.delete()

        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def bgzip(self, delete_src=False):

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

        if delete_src:
            self.delete()

        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def index(self):
        log_filepath = self.tmp_dir / f'{self.filepath.name}.csi.log'

        cmd = (''
               f'bcftools index'
               f'      {self.filepath}'
               f'      &> {log_filepath}'
               '')

        execute(cmd)
        return self

    def fix_header(self, genome_index_filepath, delete_src=False):

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

        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def copy_to(self, output_file):
        output_file = Path(output_file)

        shutil.copy2(self.filepath, output_file)
        shutil.copy2(self.filepath.with_suffix('.bgz.csi'),
                     output_file.with_suffix('.bgz.csi'))

    def drop_id(self, delete_src=False):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-id.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools annotate'
               f'      -x ID'
               f'      -O z'
               f'      -o {output_filepath}'
               f'      --threads {self.n_threads}'
               f'      {input_filepath}'
               f'      &> {log_filepath}'
               '')
        execute(cmd)

        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def drop_info(self, delete_src=False):
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
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def trim_alts(self, delete_src):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-trim_alt.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'   -O z'
               f'   -o {output_filepath}'
               f'   --threads {self.n_threads}'
               f'   -a '
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def subset(self, expression, op_name, delete_src=False):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            f'-{op_name}.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'   --include "{expression}"'
               f'   -O z'
               f'   -o {output_filepath}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)

        if delete_src:
            self.delete()

        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def subset_biallelics(self, delete_src=False):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-bi.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools view'
               f'   -M 2'
               f'   -O z'
               f'   -o {output_filepath}'
               f'   --threads {self.n_threads}'
               f'   {input_filepath}'
               f'   &> {log_filepath}'
               '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def subset_variants(self, chrom, start=None, stop=None, delete_src=False):

        input_filepath = self.filepath
        if start and stop:

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
        else:
            output_filepath = self.tmp_dir / self.filepath.name.replace(
                '.vcf.bgz',
                f'-{chrom}.vcf.bgz',
            )
            log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

            cmd = (''
                   f'bcftools view'
                   f'      -r "{chrom}"'
                   f'      -O z'
                   f'      -o {output_filepath}'
                   f'      --threads {self.n_threads}'
                   f'      {input_filepath}'
                   f'      &> {log_filepath}'
                   '')

        execute(cmd)
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def subset_samples(self, sample_file, delete_src=False):
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
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def uppercase(self, delete_src=False):
        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-uppercase.vcf',
        )
        with gzip.open(input_filepath,
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

        if delete_src:
            self.delete()

        vcf = Vcf(output_filepath, self.tmp_dir, self.n_threads)

        return vcf.bgzip(delete_src)

    def normalize(self,
                  genome_filepath,
                  atomize=False,
                  split_multiallelics=False,
                  delete_src=False):

        input_filepath = self.filepath
        output_filepath = self.tmp_dir / self.filepath.name.replace(
            '.vcf.bgz',
            '-norm.vcf.bgz',
        )
        log_filepath = self.tmp_dir / f'{output_filepath.name}.log'

        cmd = (''
               f'bcftools norm'
               f'      -f {genome_filepath}'
               f'      -c s'
               f'      -O z'
               f'      -o {output_filepath}'
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

        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def fill_tags(self, tags=['AC', 'AN', 'AF', 'NS'], delete_src=False):

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
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def drop_gt(self, delete_src=False):
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
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    def sort(self, delete_src=False):
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
        if delete_src:
            self.delete()
        return Vcf(output_filepath, self.tmp_dir, self.n_threads)

    @staticmethod
    def _fix_duplicates(input_filepath, output_filepath, columns):

        with gzip.open(input_filepath,
                       'rt') as ifd, output_filepath.open('wt') as ofd:

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

    def annotate(self, annotations_vcf, *columns, delete_src=False):
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

        output_filepath = self.tmp_dir / tmp2_filepath.name.replace(
            '.vcf', '.vcf.bgz')

        log_filepath = self.tmp_dir / f'{tmp2_filepath.name}.log'
        cmd = (''
               f' bgzip -c -@ {self.n_threads} {tmp2_filepath}'
               f' > {output_filepath}'
               f' 2> {log_filepath}'
               '')

        execute(cmd)

        if delete_src:
            self.delete()
            tmp1_filepath.unlink()
            tmp2_filepath.unlink()

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
        input_filepath = Path(f'{self.filepath}.csi')
        if not input_filepath.exists():
            self.index()

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

    if _is_gzipped(vcf):
        with gzip.open(vcf, 'rt') as fd:
            return fetch_header(fd)
    else:
        with open(vcf, 'rt') as fd:
            return fetch_header(fd)


def _is_gzipped(filepath):
    with open(filepath, 'rb') as fd:
        magic_number = fd.read(2)
        if magic_number == b'\x1f\x8b':
            return True
        return False


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

    if _is_gzipped(vcf):
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

            if proc.returncode:
                raise Exception(cmd)
