# annotate_vcf
Repo that contains code to annotate vcf files, with example input and output files and a file describing the fields in the output.

Each variant in the vcf file is annotated with relevant information including variant type, effect, sequencing coverage at site etc.

### To run the code simply type: <br />
*python annotate_vcf.py -i <input_vcf_file> -o <output_annotation_file>*

### Dependencies: <br />
1. python 3.9 <br />
2. Ensembl Variant Effect Predictor (VEP): More information can be found at: https://uswest.ensembl.org/info/docs/tools/vep/index.html
3. Python modules: <br />
numpy, subprocess, vcf, requests, collections, argparse and pytest for unit tests

### To run unit tests: <br />
*pytest test_annotate_vcf.py*

### An example vcf file is provided with example annotation. <br />
example vcf input file: Challenge_data.vcf <br />
example annotation output file produced by running annotate_vcf.py: challenge_annotation.txt

### Description of output <br />
output_description.txt contains information about the columns in the output file.
