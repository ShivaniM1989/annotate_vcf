# annotate_vcf
Repo that contains code to annotate vcf files.

Each variant in the vcf file is annotated with relevant information including variant type, effect, sequencing coverage at site etc.

To run the code simply type: <br />
python annotate_vcf.py -i <input_vcf_file> -o <output_annotation_file>

Dependencies: <br />
1. python 3.9 <br />
2. Ensembl variance effect predictor : More information can be found at: https://uswest.ensembl.org/info/docs/tools/vep/index.html

An example vcf file is provided with example annotation. <br />
example vcf input file: Challenge_data.vcf <br />
example annotation output file produced by running annotate_vcf.py: challenge_annotation.txt
