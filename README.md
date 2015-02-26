# gRNA_Tool
Tool For gRNA analysis


## CLI

```bash
usage: fastq_cli.py [-h] [--fwd FWD] [--rev REV] [-o OUTPUT] wtlib

CLI tool for gRNA analysis

positional arguments:
  wtlib                 wtlib file

optional arguments:
  -h, --help            show this help message and exit
  --fwd FWD             forward read file
  --rev REV             reverse read file
  -o OUTPUT, --output OUTPUT
                        output csv file name
```

```
python fastq_cli.py wtLib.fa --fwd=forward.fastq.gz --rev=reverse.fastq.gz
```
