# ARELIA

ARELIA is a python script that efficiently calculates residue and column reliability scores for given multiple sequence alignments, and masks the alignments based on the reliability scores. Currently, only amino acid alignments are supported.
* A free web server: http://arelia-web.appspot.com/


## Installation

ARELIA is cross-platform software if the following dependency is satisfied.

* `Python 2.7 or higher` (available at https://www.python.org/).
* `Numpy 1.6 or higher` (available at http://www.numpy.org/).

Download `arelia.py` somewhere in your `PATH`, and make `arelia.py` executable. For example, in UNIX-like systems:
```
chmod +x /path/to/your/arelia.py
```
Check it by printing a help page.
```
arelia.py -h
```
Or, just use it with `python` command.
```
python /path/to/your/arelia.py -h
```

## Examples

Get a residue-masked MSA.
```
arelia.py MSA_FILE_IN > MSA_FILE_OUT
```
or
```
arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT
```

Get a column-masked MSA.
```
arelia.py MSA_FILE_IN -msa_col MSA_FILE_OUT
```

Get a residue-masked MSA and residue reliability scores.
```
arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT -scr_res SCORE_FILE_OUT
```

Process MSA files in a directory recursively.
```
arelia.py MSA_DIR_IN -msa_res MSA_DIR_OUT -scr_res SCORE_DIR_OUT
```

Set cutoff (for alignment masking), gap penalty, and window sizes.
```
arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT -cutoff 0.3 -gap -5.0 -W 5 10 15 30
```

Set input and output formats.
```    
arelia.py MSA_DIR_IN -msa_res MSA_DIR_OUT -infmt fasta -outfmt phylip
```





