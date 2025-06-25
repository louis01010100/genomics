import polars as pl
from genomics.gregion import GenomicRegion
from genomics.gregion import create_genomic_regions

def validate_cnv(predictions, truths, reciprocal_overlap_cutoff = 0.5, boundary_difference_cutoff = 10000):

    data_ppv = _validate_cnv(
            predictions, 
            create_genomic_regions(truths), 
            'prediction', 
            'truth', 
            reciprocal_overlap_cutoff, 
            boundary_difference_cutoff,
    )

    data_sensitivity = _validate_cnv(
            truths, 
            create_genomic_regions(predictions), 
            'truth', 
            'prediction', 
            reciprocal_overlap_cutoff,
            boundary_difference_cutoff,
    )


    return {
        'data_ppv': data_ppv,
        'data_sensitivity': data_sensitivity,
    }



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
            result[f'{qname}_cn_state'] = query['cn_state']
            result[f'{dbname}_idx'] = None
            result[f'{dbname}_start'] = None
            result[f'{dbname}_end'] = None
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
            result[f'{qname}_cn_state'] = query['cn_state']
            result[f'{dbname}_idx'] = match['idx']
            result[f'{dbname}_cn_state'] = match['cn_state']
            result[f'{dbname}_start'] = match['start']
            result[f'{dbname}_end'] = match['end']
            result['reciprocal_overlap'] = reciprocal_overlap
            result['boundary_difference'] = boundary_difference
            result['reciprocal_overlap_test'] = 'PASS' if reciprocal_overlap >= reciprocal_overlap_cutoff else 'FAIL'
            result['breakpoint_tolerance_test'] = 'PASS' if boundary_difference <= boundary_difference_cutoff else 'FAIL'


            bag.append(result)

    report = pl.from_dicts(bag, infer_schema_length = None)

    return report




