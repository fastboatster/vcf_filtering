## README

### Requirements

python3 env with PyVCF3 installed

### How to run

python3 main.py vcf-file-path > output-file-path
i.e.:
python3 main.py test.chr22.ann.vcf > result.txt

### Corner cases/Limitations
Cases when vcf records have multiple different deletions are not handled in general way by this code (correct ALTs
is not output in these cases). However, in this particular vcf file, none of these deletions affect protein structure
and are not printed out so this doesn't affect output correctness.
