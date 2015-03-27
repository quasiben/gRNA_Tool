# gRNA_Tool
Tool For gRNA analysis


## CLI

```bash
usage: fastq_cli.py [-h] [--fwd FWD] [--rev REV] [-d DIR] [--parallel] wtlib

CLI tool for gRNA analysis

positional arguments:
  wtlib              wtlib file

optional arguments:
  -h, --help         show this help message and exit
  --fwd FWD          forward read file
  --rev REV          reverse read file
  -d DIR, --dir DIR  read directory of files. Files must be of match: *R1*.gz,
                     *R2*.gz
  --parallel         multicore support with Dask
```
## Examples

```
python fastq_cli.py wtLib.fa --fwd=forward.fastq.gz --rev=reverse.fastq.gz
python fastq_cli.py library.fa -d 150226_D00108_0322_AC6F9UANXX_Project_mcmanusm-MB11/Sample_Baseline/ --parallel

```
