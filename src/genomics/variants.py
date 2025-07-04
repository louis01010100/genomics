from .genome import Genome
import logging
from .variant import Variant, sync
from .vcf import fetch_variants
import polars as pl
from pathos.multiprocessing import ProcessPool

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

# data1_file must contains the VCF 4-tuple
# dat2_file must be an indexed VCF file with the ID column populated
def varmatch(data1_file, data2_file, genome_file, output_file, n_threads, batch_size):

    def batches(records, data2_file, genome_file):

        genome = Genome(genome_file)

        bag = list()
        for record in records:
            bag.append(record)

            if len(bag) >= batch_size:
                yield {
                        'records': bag, 
                        'genome_file': genome_file, 
                        'data2_file': data2_file
                }
                bag = list()

        if len(bag) > 0:
            yield {
                    'records': bag, 
                    'genome_file': genome_file, 
                    'data2_file': data2_file
            }

    def check_relation(v1, v2):

        assert v1.chrom == v2.chrom, f'{v1}\n{v2}'
        assert v1.pos == v2.pos, f'{v1}\n{v2}'
        assert v1.ref == v2.ref, f'{v1}\n{v2}'

        if v1.alts == v2.alts:
            return 'same'

        if set(v1.alts) & set(v2.alts):
            return 'overlap'

        return 'different'

    def process(batch):


        genome = Genome(batch['genome_file'])
        data2_file = batch['data2_file']
        result_collector = list()

        for record in batch['records']:
            chrom_seq = genome.seq(record['chrom'])

            ref_snv = Variant(
                chrom=record['chrom'],
                pos=int(record['pos']),
                ref=record['ref'],
                alt=record['alt'],
            )

            snvs = fetch_variants(
                record['chrom'],
                record['pos'],
                data2_file,
            )

            bag = dict()
            for snv in snvs:
                if snv.alts == ['.']:
                    continue
                try:
                    synced_ref, synced_other = sync(
                        vx = ref_snv,
                        vy = snv,
                        chrom_seq = chrom_seq
                    )
                except RecursionError as e:
                    print(e)
                    print(ref_snv)
                    print(snv)
                    continue

                result = check_relation(synced_ref, synced_other)

                if result not in bag:
                    bag[result] = list()

                bag[result].append({
                    'snv': synced_ref,
                    'matched_snv': synced_other
                })

            snv_same = list()
            snv_overlap = list()

            if 'same' in bag:
                snv_same.extend(bag['same'])
                snv_overlap.extend(bag['same'])

            if 'overlap' in bag:
                snv_overlap.extend(bag['overlap'])

            if snv_same:
                snv_same = [x['matched_snv'] for x in snv_same]
                snv_same = concat_snv(snv_same)

                record['matched_id'] = snv_same['id']
                record['matched_pos'] = snv_same['pos']
                record['matched_ref'] = snv_same['ref']
                record['matched_alt'] = snv_same['alt']
            else:
                record['matched_id'] = None
                record['matched_pos'] = None
                record['matched_ref'] = None
                record['matched_alt'] = None

            if snv_overlap:
                snv_ext = [x['matched_snv'] for x in snv_overlap]
                snv_ext = concat_snv(snv_ext)
                record['extended_matched_id'] = snv_ext['id']
                record['extended_matched_pos'] = snv_ext['pos']
                record['extended_matched_ref'] = snv_ext['ref']
                record['extended_matched_alt'] = snv_ext['alt']
            else:
                record['extended_matched_id'] = None
                record['extended_matched_pos'] = None
                record['extended_matched_ref'] = None
                record['extended_matched_alt'] = None

            result_collector.append(record)
        return result_collector

    bag = []

    data1 = pl.read_csv(data1_file, has_header = True, separator = '\t')

    if n_threads == 1:
        for batch in batches(data1.to_dicts(), data2_file):
            result = process(batch)
            bag.extend(result)

    else:
        with ProcessPool(n_threads) as pool:
            records = data1.to_dicts()
            n_done = 0
            n_total = len(records)
            for results in pool.uimap(process, batches(records, data2_file, genome_file)):
                bag.extend(results)
                n_done += len(results)
                print(f'progress: {(n_done / n_total) * 100:.2f}%',
                      end='\r',
                      flush=True)

    result =  pl.from_dicts(bag, infer_schema_length=None)
    result.write_csv(output_file, include_header = True, separator = '\t')


def concat_snv(snv: list):
    if len(snv) == 1:
        return {
            'id': snv[0].id,
            'pos': snv[0].pos,
            'ref': snv[0].ref,
            'alt': snv[0].alt,
        }

    id_ = list()
    pos = list()
    ref = list()
    alt = list()

    for x in snv:

        id_.append(x.id)
        pos.append(x.pos)
        ref.append(x.ref)
        alt.append(x.alt)

    id_ = ';'.join([str(x) for x in id_])
    pos = ';'.join([str(x) for x in pos])
    ref = ';'.join([str(x) for x in ref])
    alt = ';'.join([str(x) for x in alt])

    return {
        'id': id_,
        'pos': pos,
        'ref': ref,
        'alt': alt,
    }
