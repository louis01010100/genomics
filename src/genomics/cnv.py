import polars as pl
from genomics.gregion import GenomicRegion
from genomics.gregion import create_database
from importlib import resources


COMPLEX_REGIONS_FILE = resources.files('genomics') / 'complex_regions-hg38.tsv'
CHROM_LENGTH_HG38_FILE = resources.files('genomics') / 'chrom_lengths-hg38.tsv'

REQUIRED_COLUMNS={
    'chrom',
    'start',
    'end',
    'state',
}


def validate(
        predictions_file,
        truths_file,
        chrom_length_file = None,
        reciprocal_overlap_cutoff = 0.5, 
        boundary_difference_cutoff = 10000, 
        window_size = 100000,
        step_size = 1000,
    ):

    predictions = pl.read_csv(predictions_file, has_header = True, separator = '\t', infer_schema = False)
    truths = pl.read_csv(truths_file, has_header = True, separator = '\t', infer_schema = False)

    assert REQUIRED_COLUMNS.issubset(predictions.columns), predictions.columns

    predictions = predictions.with_columns(
            pl.col('start').cast(pl.Int64), 
            pl.col('end').cast(pl.Int64),
            pl.arange(0, length(predictions)).alias('idx'),
    )
    truths = truths.with_columns(
            pl.col('start').cast(pl.Int64), 
            pl.col('end').cast(pl.Int64)
            pl.arange(0, length(truths)).alias('idx'),
    )

    complex_regions = create_database(
        pl.read_csv(COMPLEX_REGIONS_FILE, has_header = True, separator = '\t').to_dicts()
    )

    prediction_regions, prediction_fragemnts = annotate_complex_regions(predictions, complex_regions)
    truth_regions, truth_fragments = annotate_complex_regions(truths, complex_regions)

    data_ppv = _validate_cnv(
            predictions, 
            create_database(truths), 
            'prediction', 
            'truth', 
            reciprocal_overlap_cutoff, 
            boundary_difference_cutoff,
    )

    data_sensitivity = _validate_cnv(
            truths, 
            create_database(predictions), 
            'truth', 
            'prediction', 
            reciprocal_overlap_cutoff,
            boundary_difference_cutoff,
    )

    calculate_moving_average(data_ppv, chrom_length_file, window_size, step_size)



    return {
        'data_ppv': data_ppv,
        'data_sensitivity': data_sensitivity,
    }


def moving_average(data_ppv, chrom_length_file, window_size, step_size):
    chrom_lenghs = load_chrom_lengths(chrom_length_file)


    pass



def load_chrom_lengths(input_file):

    bag = list()
    with input_file.open('rt') as fh:
        for line in fh:
            items = line.strip().split('\t')
            chrom = items[0]
            length = items[1]

            if chrom == 'chrom':
                continue

            bag.append({'chrom': chrom, 'length': int(length)})

    return pl.from_dicts(bag)



def _validate_cnv(queries, database, qname, dbname, reciprocal_overlap_cutoff, boundary_difference_cutoff):
    bag = list()
    for query in queries:
        matches = database.find_overlap(query)
        query_region = GenomicRegion(query['chrom'], query['start'], query['end'])


        if len(matches) == 0:
            result = dict()
            result[f'{qname}_idx'] = query['idx']
            result[f'chrom'] = query['chrom']
            result[f'{qname}_start'] = query['start']
            result[f'{qname}_end'] = query['end']
            result[f'{qname}_length'] = query['end'] - query['start']
            result[f'{qname}_cn_state'] = query['cn_state']
            result[f'{dbname}_idx'] = None
            result[f'{dbname}_start'] = None
            result[f'{dbname}_end'] = None
            result[f'{dbname}_length'] = None
            result[f'{dbname}_cn_state'] = None
            result['reciprocal_overlap'] = None
            result['boundary_difference'] = None
            result['reciprocal_overlap_test'] = None
            result['breakpoint_tolerance_test'] = None
            bag.append(result)
            continue

        for match in matches:
            result = dict()
            db_region = GenomicRegion(match['chrom'], match['start'], match['end'])
            reciprocal_overlap = query_region.get_reciprocal_overlap(db_region)
            boundary_difference = query_region.get_max_boundary_difference(db_region)

            result[f'{qname}_idx'] = query['idx']
            result[f'chrom'] = query['chrom']
            result[f'{qname}_start'] = query['start']
            result[f'{qname}_end'] = query['end']
            result[f'{qname}_length'] = query['end'] - query['start']
            result[f'{qname}_cn_state'] = query['cn_state']
            result[f'{dbname}_idx'] = match['idx']
            result[f'{dbname}_start'] = match['start']
            result[f'{dbname}_end'] = match['end']
            result[f'{dbname}_length'] = match['end'] - match['start']
            result[f'{dbname}_cn_state'] = match['cn_state']
            result['reciprocal_overlap'] = reciprocal_overlap
            result['boundary_difference'] = boundary_difference
            result['reciprocal_overlap_test'] = 'PASS' if reciprocal_overlap >= reciprocal_overlap_cutoff else 'FAIL'
            result['breakpoint_tolerance_test'] = 'PASS' if boundary_difference <= boundary_difference_cutoff else 'FAIL'

            bag.append(result)

    report = pl.from_dicts(bag, infer_schema_length = None)

    return report

def annotate_complex_regions(cnvs, complex_region):


    cnvs = pl.read_csv(cnv_file, has_header = True, separator = '\t', infer_schema = False)
    expected_columns = {'chrom', 'start', 'end'}
    assert expected_columns.issubset(set(cnvs.columns), cnv.columns

    cnvs = cnv.with_columns(
            pl.col('start').cast(pl.Int64).alias('start'),
            pl.col('end').cast(pl.Int64).alias('end'),
    )

    bag_region = list()
    bag_fragment = list()

    for record in cnvs.to_dicts():
        matches = complex_regions.find_overlap(record)

        if len(matches) == 0:
            record['complex_region'] = None
            bag_region.append(record)
            bag_fragment.append(record)
        else:


            for match in matches:
                query = record.copy()
                match_region = GenomicRegion(match['chrom'], match['start'], match['end'])
                query_region = GenomicRegion(query['chrom'], query['start'], query['end'])

                query['complex_region'] = match['category']
                query['reciprocal_overlap'] = query_region.get_reciprocal_overlap(match_region)

                bag_region.append(query)

                data = record.copy()
                if record['start'] < match['start']:
                    data['start'] = record['start']
                    data['end'] = match['start']
                    bag_fragment.append(data)
                if record['end'] > match['end']:
                    data['start'] = match['end']
                    data['end'] = record['end']
                    bag_fragment.append(data)
    regions = pl.from_dicts(bag_region)
    fragments = pl.from_dicts(bag_fragment)

    return regions, fragments



