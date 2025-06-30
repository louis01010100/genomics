from .genome import Genome
from .variant import Variant
import polars as pl

def assign_coord_id(df, genome_file, n_threads = 1):
    target_columns = ['chrom', 'pos', 'ref', 'alt', 'id']

    for column in target_columns:
        if column not in df.columns:
            print(target_columns, df.columns)
            raise Exception(column)

    genome = Genome(genome_file)

    bag = list()


    for record in df.to_dicts():
        chrom = record['chrom']
        pos = record['pos']
        ref = record['ref']
        alt = record['alt']
        id_ = record['id']
        seq = genome.seq(chrom)

        v = Variant(chrom, pos, ref, alt, id_ = id_).to_vcf(seq).normalize(seq)

        bag.append({
            'chrom': v.chrom,
            'pos': v.pos,
            'ref': v.ref,
            'alt': v.alt,
            'id': v.id,
        })

    data = pl.from_dicts(bag)

    bag = list()
    for keys, df in data.group_by(['chrom', 'pos', 'ref']):
        coord_id = '_'.join([str(x) for x in keys])
        df = df.with_columns(pl.lit(coord_id).alias('coord_id'))
        bag.append(df)

    result = pl.concat(bag)

    return result
