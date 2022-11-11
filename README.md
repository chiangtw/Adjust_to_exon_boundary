## Usage

```
usage: adjust_to_exon_boundary.py [-h] [--dist DIST] [--na_value NA_VALUE] anno_db region_file

positional arguments:
  anno_db
  region_file          4-columns TSV file: (chrm, pos1, pos2, strand)

options:
  -h, --help           show this help message and exit
  --dist DIST          Maximum distances allowed for adjusting positions. (default: 5)
  --na_value NA_VALUE  Placeholder for na values. (default: 'NA')

```

### Example

```
# Generate the annotation db
circmimi_tools gendb annotation.gtf annotation.db
```

```
cat circRNAs.tsv | ./adjust_to_exon_boundary.py annotation.db - > circRNAs.adjust.tsv
```

