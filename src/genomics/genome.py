from .gregion import GenomicRegion
import re
from pathlib import Path

# HG38_X_PAR_1 = GenomicRegion('chrX', 10000, 2781479)
# HG38_X_PAR_2 = GenomicRegion('chrX', 2781479, 156030895)
# HG38_Y_PAR_1 = GenomicRegion('chrY', 10000, 2781479)
# HG38_Y_PAR_2 = GenomicRegion('chrY', 56887902, 57217415)
#
# HG19_X_PAR_1 = GenomicRegion('chrX', 60001, 2699520)
# HG19_X_PAR_2 = GenomicRegion('chrX', 154931044, 155260560)
# HG19_Y_PAR_1 = GenomicRegion('chrY', 10001, 2649520)
# HG19_Y_PAR_2 = GenomicRegion('chrY', 59034050, 59363566)


class Genome():

    def __init__(self, genome_fh):
        chroms = dict()
        id_ = None
        seq = None

        for line in genome_fh:
            line = line.strip()
            if line.startswith('>'):
                if id_:
                    chroms[id_] = seq
                match = re.match(r'>([^ ]+)(:? .+)?$', line)
                id_ = match.group(1)
                seq = ''
            else:
                seq += line.upper()
        if id_:
            chroms[id_] = seq

        self._chroms = chroms

    @property
    def chroms(self):
        return set(self._chroms.keys())

    def rename_chroms(self, chrom_name_map):

        chroms = dict()

        for old_name, seq in self._chroms.items():
            if old_name in chrom_name_map:
                chrom_name = chrom_name_map[old_name]
            else:
                chrom_name = old_name
            chroms[chrom_name] = seq
        self._chroms = chroms
        return self

    @property
    def version(self):
        return self._version

    def slice(self, chrom, start, stop):
        return self._chroms[chrom][start:stop]

    def length(self, chrom):
        return len(self._chroms[chrom])

    def to_fasta(self, output_file: Path, width: int = 80):

        with output_file.open('wt') as fh:
            for chrom_name, seq in self._chroms.items():
                fh.write(f'>{chrom_name}\n')
                begin_pos = 0
                end_pos = begin_pos + width

                while begin_pos < len(seq):
                    fragment = seq[begin_pos:end_pos]
                    fh.write(f'{fragment}\n')
                    begin_pos = end_pos
                    end_pos = begin_pos + width

                    # python handles exceeding end_pos well

    # def ploidy(self, chrom, pos):
    #     v = GenomicRegion(chrom.upper, pos, pos)
