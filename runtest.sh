input_file="test/vcf/samplechr22_76g_PASS.vcf.gz"
output_directory="$(dirname $input_file)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run main.nf \
	--input_file $input_file \
	--output_dir $output_directory \
	--gff_ref="test/gff/hsa.gff3" \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
