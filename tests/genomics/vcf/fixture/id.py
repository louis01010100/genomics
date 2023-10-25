#!/usr/bin/env python
from pathlib import Path


def main():
    input_file = Path('sample.vcf')
    output_file = Path('sample2.vcf')

    with input_file.open('rt') as ifh, output_file.open('wt') as ofh:
        idx = 0
        for line in ifh:
            if line.startswith('#'):
                ofh.write(line)
                continue

            items = line.split('\t', 3)
            id_ = f'AX-{idx:06d}'
            items[2] = id_

            ofh.write('\t'.join(items))
            idx += 1


if __name__ == '__main__':
    main()
