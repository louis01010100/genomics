
def validate_cnv(predictions, truth, reciprocal_overlap_cutoff = 0.5, breakpoint_tolerance_threshold = 10000):

    bag = list()
    for p in predictions:
        founds = truth.find_overlap(p['chrom'], p['start'], p['end'])

        for found in founds:
            reciprocal_overlap = p.reciprocal_overlap(found)
            breakpoint_delta = p.max_boundary_delta(found)

            x = p.copy()

            x['truth_id'] = found['id']
            x['truth_start'] = found['start']
            x['truth_end'] = found['end']
            x['truth_state'] = found['state']
            x['reciprocal_overlap'] = reciprocal_overlap
            x['breakpoint_delta'] = breakpoint_delta
            x['same_cn_state'] = p.same_state(found)
            x['reciprocal_overlap_test'] = 'PASS' if reciprocal_overlap >= reciprocal_overlap_cutoff else 'FAIL'
            x['breakpoint_tolerance_test'] = 'PASS' if breakpoint_delta <= breakpoint_tolerance_threshold else 'FAIL'

            bag.append(x)

    report = pl.from_dicts(bag)

    return report


def calculate_reciprocal_overlap(x, y):

    gregion(x['chrom'], x['start'], x['end']).
    pass

