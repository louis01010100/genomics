
from genomics.gregion import GenomicRegion
def validate_cnv(predictions, truth, reciprocal_overlap_cutoff = 0.5, boundary_difference_cutoff = 10000):

    bag = list()
    for p in predictions:
        region_p = GenomicRegion(x['chrom'], x['start'], x['end'])
        founds = truth.find_overlap(p['chrom'], p['start'], p['end'])

        for found in founds:
            region_t = GenomicRegion(truth['chrom'], truth['start'], truth['end'])

            reciprocal_overlap = region_p.get_reciprocal_overlap(region_t)
            boundary_difference = region_p.get_max_boundary_difference(region_t)

            x = p.copy()

            x['truth_id'] = found['id']
            x['truth_start'] = found['start']
            x['truth_end'] = found['end']
            x['truth_state'] = found['state']
            x['reciprocal_overlap'] = reciprocal_overlap
            x['boundary_difference'] = boundary_difference
            x['same_cn_state'] = p.same_state(found)
            x['reciprocal_overlap_test'] = 'PASS' if reciprocal_overlap >= reciprocal_overlap_cutoff else 'FAIL'
            x['breakpoint_tolerance_test'] = 'PASS' if boundary_difference <= boundary_difference_cutoff else 'FAIL'

            bag.append(x)

    report = pl.from_dicts(bag)

    return report




