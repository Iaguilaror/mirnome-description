#!/usr/bin/env nextflow

/*================================================================
The AGUILAR LAB presents...

  The miRNome variants population description pipeline

- A variant finder tool

==================================================================
Version: 0.0.1
Project repository: https://github.com/Iaguilaror/mirnome-description
==================================================================
Authors:

- Bioinformatics Design
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

- Bioinformatics Development
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

- Nextflow Port
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

=============================
Pipeline Processes In Brief:
.
Pre-processing:
  _pre1_buildbeds

Core-processing:
  _002_filtervcf
  _003_vcf2bed
  _004_annotate

Pos-processing
 _pos1_gatherbeds

Anlysis
 _an1_

ENDING
 _register_configs

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
  The miRNome variants population description pipeline
  - A variant finder tool
  v${version}
  ==========================================

	Usage:

  nextflow run main.nf --input_file <path to inputs> [--output_dir path to results ]

	  --input_file    <- To-do;
				To-do;
				To-do;
	  --output_dir     <- directory where results, intermediate and log files will bestored;
				default: same dir where --query_fasta resides
	  -resume	   <- Use cached results if the executed project has been run before;
				default: not activated
				This native NF option checks if anything has changed from a previous pipeline execution.
				Then, it resumes the run from the last successful stage.
				i.e. If for some reason your previous run got interrupted,
				running the -resume option will take it from the last successful pipeline stage
				instead of starting over
				Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show Pipeline version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "miRNome variants"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.input_file = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "Pipeline v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.10.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
  TODO (iaguilar) check the extension of input queries; see getExtension() at https://www.nextflow.io/docs/latest/script.html#check-file-attributes
*/

/* Check if inputs provided
    if they were not provided, they keep the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.input_file ) {
  log.error " Please provide the following params: --input_file \n\n" +
  " For more information, execute: nextflow run main.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.input_file).getParent()

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable (pipeline_name) defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The miRNome variants population description pipeline
- A variant finder tool
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['input file']			= params.input_file
pipelinesummary['input gff']			= params.gff_ref
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/*
	READ INPUTS
*/

/* Load vcf file into channel */
Channel
  .fromPath( "${params.input_file}" )
  .set{ vcf_input }

/* Load gff file into channel */
Channel
  .fromPath( "${params.gff_ref}" )
  .set{ gff_input }

/* _pre1_buildbeds */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/01-buildbeds/*")
	.toList()
	.set{ mkfiles_pre1 }

process _pre1_buildbeds {
	label 'standard'

	publishDir "${results_dir}/_pre1_buildbeds/",mode:"symlink"

	input:
  file gff from gff_input
  file mk_files from mkfiles_pre1

	output:
  file "*_mirBase_22.bed" into results_pre1_buildbeds
  file "complete_mirna.bed" into results_pre1_completebeds

	"""
	bash runmk.sh
	"""

}

/* _002_filtervcf */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/02-filtervcf/*")
	.toList()
	.set{ mkfiles_002 }

process _002_filtervcf {
	label 'standard'

	publishDir "${results_dir}/_002_filtervcf/",mode:"symlink"

	input:
  file vcf from vcf_input
  file bedref from results_pre1_completebeds
  file mk_files from mkfiles_002

	output:
  file "*.filtered.vcf" into results_002_filtervcf

	"""
  export BEDFILE=${bedref}
	bash runmk.sh
	"""

}

/* _003_vcf2bed */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/03-vcf2bed/*")
	.toList()
	.set{ mkfiles_003 }

process _003_vcf2bed {
	label 'standard'

	publishDir "${results_dir}/_003_vcf2bed/",mode:"symlink"

	input:
  file vcf from results_002_filtervcf
  file mk_files from mkfiles_003

	output:
  file "*.bed" into results_003_vcf2bed

	"""
	bash runmk.sh
	"""

}

/* _004_annotate */
/* gather filtered vcfbed and mirnome beds */
results_pre1_buildbeds
  .combine( results_003_vcf2bed )
  .set{ inputs_004 }

/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/04-annotatebybed/*")
	.toList()
	.set{ mkfiles_004 }

process _004_annotate {
	label 'standard'

	publishDir "${results_dir}/_004_annotate/",mode:"symlink"

	input:
  file inputbeds from inputs_004
  file mk_files from mkfiles_004

	output:
  file "*.variants.bed" into results_004_annotate

	"""
  export BASEBED="\$(ls *.filtered.bed)"
	bash runmk.sh
	"""

}

/* _pos1_gatherbeds */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/05-gatherbeds/*")
	.toList()
	.set{ mkfiles_pos1 }

process __pos1_gatherbeds {
	label 'standard'

	publishDir "${results_dir}/_pos1_gatherbeds/",mode:"symlink"

	input:
  file inputbeds from results_004_annotate
  file mk_files from mkfiles_pos1

	output:
  file "allvariants.bed"
  file "variants_of_interest_by_type.csv"
  file "*.svg"

	"""
  export BEDDIR="."
	bash runmk.sh
	"""

}
