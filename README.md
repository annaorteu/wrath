# wrath
WRATH: WRapped Analysis of Tagged Haplotypes

<img src="wrath_logo.png" alt="logo" width="50%"/>


---
## Running wrath

The main script `wrath.sh` runs the main version of *WRATH*.

A typical command looks like:

```bash
wrath.sh -g genome_file.fa -c chromosome_name  -w 50000  -s list_of_bam_files.txt -t 15
```

### Requirements

Command line programs:
- faidx
- samtools
- bedtools

Python:
- NumPy
- Seaborn
- [pandas](https://pandas.pydata.org/) 
- [sklearn](https://scikit-learn.org/stable/index.html)

```bash
pip install -U scikit-learn
```

import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering


### Input
