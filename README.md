# annotate_vcf
Repo that contains code to annotate vcf files
Each variant in the vcf file is annotated with relevant information including variant type, effect, sequencing coverage at site etc.

To run the code simply type:
python annotate_vcf.py -i <input_vcf_file> -o <output_annotation_file>

Dependencies:
python 3.9
Ensembl variance effect predictor : More information can be found at: https://uswest.ensembl.org/info/docs/tools/vep/index.html

An example vcf file is provided with example annotation.
example vcf input file: Challenge_data.vcf
example annotation output file produced by running annotate_vcf.py: challenge_annotation.txt
