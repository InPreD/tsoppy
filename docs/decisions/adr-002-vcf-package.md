# vcf package

## Context and Problem Statement

There are several vcf reading/writing packages available, we want to select one to work with across the `tsoppy` code, if possible. 

Considering:
 - performance, 
 - ease of use, 
 - active maintenance/development, 
 - proper version tagging/releasing, 
 - numeric issues so that we really do what we aim to...

Three packages listed below were selected based on active development/maintenance and proper version tagging/releasing.

package | code | documentation
---|---|---
CyVCF2 | [github](https://github.com/brentp/cyvcf2) | [documentation](https://brentp.github.io/cyvcf2/)
pysam  | [github](https://github.com/pysam-developers/pysam) | [readthedocs](https://pysam.readthedocs.io/en/latest/)
vcfpy  | [github](https://github.com/bihealth/vcfpy) | [readthedocs](https://vcfpy.readthedocs.io/)

The packages were further tested for performance and usability, the numeric issue came up after the `pysam` package 
returned after reading in the test vcf slightly different values than what the input vcf file contained. The difference is caused
by binary representation of floating points in `pysam` while `CyVCF2` and `vcfpy` use decimal representation meaning that the two
packages read and store exactly the value that is written in the input vcf. There is a concern that using any threshold value
to filter variants may by accident cause values above/below the threshold to pass even though they should not. We have no way to
correct this behaviour so at the moment, the two other packages are favoured.

### Testing performance in BRCA1 region

In this case, allele frequency for alternative allele of identified variants was printed out.

Used code:

```python
import time
from pathlib import Path

import cyvcf2
import pysam
import vcfpy

# test vcf file path
p = Path('../TSO500_V2_testing_output/Logs_Intermediates/DnaDragenCaller/TVD0016-D01-P01-A00/TVD0016-D01-P01-A00.hard-filtered.gvcf.gz')
BRCA1='chr17:41196312-41277381'


###CyVCF2

start_time = time.perf_counter()

vcf = cyvcf2.VCF(p)

for variant in vcf(BRCA1):
    afs = variant.format('AF')
    if afs is not None:
        print(afs[0,0])

end_time = time.perf_counter()
print(f"CyVCF2: {end_time - start_time}")


###pysam

start_time = time.perf_counter()

variant_file = pysam.VariantFile(p, 'r')
variants_in_region = variant_file.fetch(region=BRCA1)

for variant_record in variants_in_region:
    for sample_name in variant_record.samples:
        variant_record_sample = variant_record.samples[sample_name]
        afs = variant_record_sample.get('AF')
        if afs is not None:
            print(afs[0])

variant_file.close()

end_time = time.perf_counter()
print(f"pysam: {end_time - start_time}")


###vcfpy

start_time = time.perf_counter()

reader = vcfpy.Reader.from_path(p)
variants_in_region = reader.fetch(BRCA1)

for variant in variants_in_region:
    call_for_sample = variant.call_for_sample
    
    for sample_name in call_for_sample.keys():
        data = call_for_sample[sample_name].data
        if 'AF' in data.keys():
            print(data['AF'][0])

reader.close()

end_time = time.perf_counter()
print(f"vcfpy: {end_time - start_time}")
```

Output: 

```
0.0033
0.0028
0.01
0.0
0.0
0.0189
0.0162
0.008
0.0016
0.002
0.0215
0.0027
0.0081
0.1274
0.0142
0.0566
0.0046
0.0045
0.0053
0.0043
0.1286
0.0207
0.112
0.0124
0.0018
CyVCF2: 0.01925382000626996
0.0032999999821186066
0.00279999990016222
0.009999999776482582
0.0
0.0
0.01889999955892563
0.016200000420212746
0.00800000037997961
0.0015999999595806003
0.0020000000949949026
0.0215000007301569
0.0027000000700354576
0.008100000210106373
0.1273999959230423
0.0142000000923872
0.05660000070929527
0.004600000102072954
0.0044999998062849045
0.0052999998442828655
0.00430000014603138
0.12860000133514404
0.02070000022649765
0.1120000034570694
0.012400000356137753
0.0017999999690800905
pysam: 0.07340225600637496
0.0033
0.0028
0.01
0.0
0.0
0.0189
0.0162
0.008
0.0016
0.002
0.0215
0.0027
0.0081
0.1274
0.0142
0.0566
0.0046
0.0045
0.0053
0.0043
0.1286
0.0207
0.112
0.0124
0.0018
vcfpy: 0.22303411300526932
```

### Testing performance on all records

In this case, only assignment of allele frequency ('AF') was performed, no printing.

Code:
```python
import time
from pathlib import Path

import cyvcf2
import pysam
import vcfpy

# test vcf file path
p = Path('../TSO500_V2_testing_output/Logs_Intermediates/DnaDragenCaller/TVD0016-D01-P01-A00/TVD0016-D01-P01-A00.hard-filtered.gvcf.gz')


###CyVCF2

start_time = time.perf_counter()

vcf = cyvcf2.VCF(p)

for variant in vcf():
    afs = variant.format('AF')

end_time = time.perf_counter()
print(f"CyVCF2: {end_time - start_time}")


###pysam

start_time = time.perf_counter()

variant_file = pysam.VariantFile(p, 'r')
variants_in_region = variant_file.fetch()

for variant_record in variants_in_region:
    for sample_name in variant_record.samples:
        variant_record_sample = variant_record.samples[sample_name]
        afs = variant_record_sample.get('AF')

variant_file.close()

end_time = time.perf_counter()
print(f"pysam: {end_time - start_time}")


###vcfpy

start_time = time.perf_counter()

reader = vcfpy.Reader.from_path(p)

for variant in reader:
    call_for_sample = variant.call_for_sample
    
    for sample_name in call_for_sample.keys():
        data = call_for_sample[sample_name].data

reader.close()

end_time = time.perf_counter()
print(f"vcfpy: {end_time - start_time}")
```


Output:

```
CyVCF2: 2.140532824967522
pysam: 7.6007183830370195
vcfpy: 15.800526580016594
```

When using `vcfpy`, one has to handle internal structures of the package classes to get out 
what is needed - not a good level of abstraction, making things unnecesarily complicated.

`CyVCF2` is fastest and seem to provide reading, writing as well as functionality to modify the header and the records.
This seems to be all the functionality we may need.

## Decision

We will start using `CyVCF2` package for working with VCF files due to its performance, ease of use, and provided sufficient functionality.

## Consequences

- install `CyVCF2` in the repo
- use `CyVCF2` for working with `vcf` files

### Positive

- performance
- code readability

### Negative

- learning curve

### Neutral



