# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Septempber-25-2020
# Email: amir.shams84@gmail.com
# Project: Eric-QC
# Aim: python script to build QC execution shell script
# python Eric_QC.py Eric-QC_metadata.csv
# swarm -g 50 -t 30 --time 04:00:00 --logdir ./logs/ --partition quick --sbatch "--mail-type=FAIL" -f Eric-QC_metadata_NGS_test_execution.swarm
# ################################### IMPORT ##################################### #


import os
import sys
import pandas
import glob
import numpy
import collections
# ################################### FUNCTIONS ################################## #


def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False


def is_dataframe_empty(samplesheet_DF):
	"""
	"""
	if samplesheet_DF.empty is True or len(samplesheet_DF.index) == 0:
		#
		return True
	else:
		return False


def parse_samplesheet_to_dict(samplesheet_file_Path, sample_column=None):
	"""
	parse input samplesheet with all parameters into dict for future access
	"""
	if is_file_exist(samplesheet_file_Path) is True:
		#
		pass
	else:
		#if is_file_exist(samplesheet_file_Path) is True:
		print("FATAL ERROR:samplesheet file path is not accessible or the premission to read is not granted!!!\nSamplesheet file Path: " + samplesheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	#
	samplesheet_Dict = {}
	#
	if sample_column is None:
		#set the default sample_id as index
		sample_column = "title"
	else:
		#set the samplesheet index
		pass
	#
	samplesheet_DF = pandas.read_csv(
		samplesheet_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter=",",
		index_col=None
	)
	#
	if is_dataframe_empty(samplesheet_DF) is True:
		#dataframe is empty
		print("FATAL ERROR: Samplesheet file is empty!!!\nSamplesheet file Path: " + samplesheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	else:
		#
		pass

	samplesheet_DF.replace(numpy.nan, '', regex=True)
		
	#
	transposed_samplesheet_DF = samplesheet_DF.transpose()
	samplesheet_Dict = transposed_samplesheet_DF.to_dict()
	#
	return samplesheet_Dict


def build_study_metadata_dict(samplesheet_Dict, study_metadata_Dict):
	"""
	"""
	for index, sample_Dict in samplesheet_Dict.items():
		##
		study_metadata_Dict[sample_Dict["title"]] = {}
		study_metadata_Dict[sample_Dict["title"]]["input"] = sample_Dict["input"]
		study_metadata_Dict[sample_Dict["title"]]["output"] = sample_Dict["output"]
		if str(sample_Dict["single_end"]) == "1":
			#
			layout = "single"
		else:
			#
			layout = "paired"
		#
		study_metadata_Dict[sample_Dict["title"]]["layout"] = layout
		
	else:
		##for index, sample_Dict in samplesheet_Dict.items():
		pass

	return True


def parse_datadir(DATA_DIR, LAYOUT):
	"""
	"""
	fastq_path_List = []
	fastq_path_List = sorted(glob.glob(DATA_DIR + "/" + "*.fastq*", recursive=True))
	fastq_path_List.extend(sorted(glob.glob(DATA_DIR + "/" + "*.fq*", recursive=True)))
	#
	if len(fastq_path_List) == 0:
		print("FATAL ERROR:no fastq file detected, please check your inputdir")

		raise Exception("ABORTING!!!")
		sys.exit(2)
	else:
		pass
	#
	fastq_file_path_Dict = {}
	for each_fastq_path in fastq_path_List:
		##
		each_fastq_basename = os.path.basename(each_fastq_path)
		if LAYOUT != "single":
			#
			reverse_flag = False
			for each_delimiter in REVERSE_DELIMITER_LIST:
				##
				if each_delimiter in each_fastq_basename:
					#
					reverse_flag = True
				else:
					#
					pass
			else:
				##for each_delimiter in REVERSE_DELIMITER_LIST:
				pass
			if reverse_flag is True:
				#
				continue
			else:
				#if reverse_flag is True:
				pass

			for each_delimiter in FORWARD_DELIMITER_LIST:
				##
				if each_delimiter in each_fastq_basename:
					#
					each_fastq_basename = each_fastq_basename.replace(each_delimiter, "")
					break
				else:
					#
					pass
			else:
				##for each_delimiter in FORWARD_DELIMITER_LIST:
				pass
		else:
			#
			pass
		sample_name = "NO_SAMPLE_NAME"

		sample_name = each_fastq_basename.split(".fastq", 1)[0].split(".fq", 1)[0]
			
		#
		if sample_name not in fastq_file_path_Dict:
			#
			fastq_file_path_Dict[sample_name] = [each_fastq_path]
		else:
			fastq_file_path_Dict[sample_name].append(each_fastq_path)

	else:
		##for each_fastq_path in fastq_path_List:
		pass

	return fastq_file_path_Dict


def process_fastq_path(forward_fastq_path):
	"""
	"""
	forward_fastq_name = os.path.basename(forward_fastq_path).split(".fastq", 1)[0].split(".fq", 1)[0]
	reverse_fastq_name = forward_fastq_name
	unstranded_fastq_name = forward_fastq_name.split(FORWARD_DELIMITER_LIST[0], 1)[0].split(FORWARD_DELIMITER_LIST[1], 1)[0].split(REVERSE_DELIMITER_LIST[0], 1)[0].split(REVERSE_DELIMITER_LIST[1], 1)[0]
	for fwd, rev in zip(FORWARD_DELIMITER_LIST, REVERSE_DELIMITER_LIST):
		##
		reverse_fastq_name = reverse_fastq_name.replace(fwd, rev)
	else:
		##
		pass
	reverse_fastq_path = forward_fastq_path.replace(forward_fastq_name, reverse_fastq_name)
	return (forward_fastq_path, reverse_fastq_path, unstranded_fastq_name)
# ################################### WORKFLOW ################################## #


def build_kraken_script(study_metadata_Dict, study):
	"""
	"""
	kraken_script = "###############################################\n"
	kraken_script += "#KRAKEN-KRONA\n"
	KRAKEN2_DB = "/fdb/kraken/20180907_standard_kraken2"
	#
	LAYOUT = study_metadata_Dict[study]["layout"]
	fastq_file_path_Dict = parse_datadir(study_metadata_Dict[study]["input"], study_metadata_Dict[study]["layout"])
		
	for each_sample in fastq_file_path_Dict:
		##
		forward_fastq_path, reverse_fastq_path, unstranded_fastq_name = process_fastq_path(fastq_file_path_Dict[each_sample][0])
		if LAYOUT == "paired":
			#
			kraken_layout_command_string = forward_fastq_path + " " + reverse_fastq_path + " --paired"
		else:
			#
			kraken_layout_command_string = forward_fastq_path
		#
		kraken_report_path = study_metadata_Dict[study]["output"] + "/" + study + "/kraken"
		kraken_report = kraken_report_path + "/" + each_sample + ".kraken.txt"
		kraken_log = kraken_report_path + "/" + each_sample + ".kraken.log"
		kraken_sample_script = "mkdir -p {kraken_report_path}; module load kraken/2.0.8-beta; kraken2 --threads $SLURM_CPUS_PER_TASK --db {kraken_db} --report {kraken_report} --gzip-compressed ".format(kraken_report_path=kraken_report_path, kraken_db=KRAKEN2_DB, kraken_report=kraken_report) + kraken_layout_command_string
		kraken_script += kraken_sample_script + " 1> {kraken_log} 2>&1 \n\n".format(kraken_log=kraken_log)
	else:
		##for each_sample in fastq_file_path_Dict:
		pass
	return kraken_script


def build_fastqc_script(study_metadata_Dict, study):
	"""
	"""
	fastqc_script = "###############################################\n"
	fastqc_script += "#FASTQC\n"

	LAYOUT = study_metadata_Dict[study]["layout"]
	fastq_file_path_Dict = parse_datadir(study_metadata_Dict[study]["input"], study_metadata_Dict[study]["layout"])
		
	for each_sample in fastq_file_path_Dict:
		##
		forward_fastq_path, reverse_fastq_path, unstranded_fastq_name = process_fastq_path(fastq_file_path_Dict[each_sample][0])

		fastqc_report_path = study_metadata_Dict[study]["output"] + "/" + study + "/fastqc"
		fastqc_log = study_metadata_Dict[study]["output"] + "/" + study + "/fastqc/" + each_sample + ".fastqc.log"
		if LAYOUT == "paired":
			#
			fastqc_sample_script = "mkdir -p {fastqc_report_path}; module load fastqc; fastqc -o {fastqc_report_path} -f fastq --threads $SLURM_CPUS_PER_TASK ".format(fastqc_report_path=fastqc_report_path) + forward_fastq_path + " 1> {fastqc_log} 2>&1 \n".format(fastqc_log=fastqc_log)
			fastqc_sample_script += "\nmkdir -p {fastqc_report_path}; module load fastqc; fastqc -o {fastqc_report_path} -f fastq --threads $SLURM_CPUS_PER_TASK ".format(fastqc_report_path=fastqc_report_path) + reverse_fastq_path + " 1>> {fastqc_log} 2>&1 \n".format(fastqc_log=fastqc_log)
		else:
			#
			fastqc_sample_script = "mkdir -p {fastqc_report_path}; module load fastqc; fastqc -o {fastqc_report_path} -f fastq --threads $SLURM_CPUS_PER_TASK ".format(fastqc_report_path=fastqc_report_path) + forward_fastq_path + " 1> {fastqc_log} 2>&1 \n".format(fastqc_log=fastqc_log)
		#
		fastqc_script += fastqc_sample_script + "\n\n"
	else:
		##for each_sample in fastq_file_path_Dict:
		pass
	
	return fastqc_script


def build_fastq_screen_script(study_metadata_Dict, study):
	fastq_screen_script = "###############################################\n"
	fastq_screen_script += "#FASTQ_SCREEN\n"
	
	LAYOUT = study_metadata_Dict[study]["layout"]
	fastq_file_path_Dict = parse_datadir(study_metadata_Dict[study]["input"], study_metadata_Dict[study]["layout"])
		
	for each_sample in fastq_file_path_Dict:
		##
		forward_fastq_path, reverse_fastq_path, unstranded_fastq_name = process_fastq_path(fastq_file_path_Dict[each_sample][0])

		fastq_screen_report_path = study_metadata_Dict[study]["output"] + "/" + study + "/fastq_screen"
		fastq_screen_log = study_metadata_Dict[study]["output"] + "/" + study + "/fastq_screen/" + each_sample + ".fastq_screen.log"
		if LAYOUT == "paired":
			#
			fastq_screen_sample_script = "mkdir -p {fastq_screen_report_path}; module load fastq_screen; fastq_screen --force --aligner bowtie2 --subset 100000 --threads $SLURM_CPUS_PER_TASK --bowtie2 --local --outdir {fastq_screen_report_path} ".format(fastq_screen_report_path=fastq_screen_report_path) + forward_fastq_path + " 1> {fastq_screen_log} 2>&1 \n".format(fastq_screen_log=fastq_screen_log)
			fastq_screen_sample_script += "\nmkdir -p {fastq_screen_report_path}; module load fastq_screen; fastq_screen --force --aligner bowtie2 --subset 100000 --threads $SLURM_CPUS_PER_TASK --bowtie2 --local --outdir {fastq_screen_report_path} ".format(fastq_screen_report_path=fastq_screen_report_path) + reverse_fastq_path + " 1> {fastq_screen_log} 2>&1 \n".format(fastq_screen_log=fastq_screen_log)
		else:
			#
			fastq_screen_sample_script = "mkdir -p {fastq_screen_report_path}; module load fastq_screen; fastq_screen --force --aligner bowtie2 --subset 100000 --threads $SLURM_CPUS_PER_TASK --bowtie2 --local --outdir {fastq_screen_report_path} ".format(fastq_screen_report_path=fastq_screen_report_path) + forward_fastq_path + " 1> {fastq_screen_log} 2>&1 \n".format(fastq_screen_log=fastq_screen_log)
		#
		
		fastq_screen_script += fastq_screen_sample_script + "\n\n"
	else:
		##for each_sample in fastq_file_path_Dict:
		pass
	
	return fastq_screen_script


def build_fastp_script(study_metadata_Dict, study):
	"""
	fastp \\
	--thread {threads} \\
	-i $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq.gz \\
	{params.paired_trim_fastp_trim_fastp_concat_cat_qc_fastp} \\
	--json $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.json \\
	--html $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.html
	"""

	fastp_script = "###############################################\n"
	fastp_script += "#FASTP\n"
	
	LAYOUT = study_metadata_Dict[study]["layout"]
	fastq_file_path_Dict = parse_datadir(study_metadata_Dict[study]["input"], study_metadata_Dict[study]["layout"])
		
	for each_sample in fastq_file_path_Dict:
		##
		forward_fastq_path, reverse_fastq_path, unstranded_fastq_name = process_fastq_path(fastq_file_path_Dict[each_sample][0])
		if LAYOUT == "paired":
			#
			fastp_layout_command_string = "-i " + forward_fastq_path + " -I " + reverse_fastq_path
		else:
			#
			fastp_layout_command_string = "-i " + forward_fastq_path
		#
		fastp_report_path = study_metadata_Dict[study]["output"] + "/" + study + "/fastp"
		fastp_report_json = fastp_report_path + "/" + each_sample + ".fastp.json"
		fastp_report_html = fastp_report_path + "/" + each_sample + ".fastp.html"
		fastp_log = fastp_report_path + "/" + each_sample + ".fastp.log"
		fastp_sample_script = "mkdir -p {fastp_report_path}; module load fastp; fastp --thread $SLURM_CPUS_PER_TASK --json {fastp_report_json} --html {fastp_report_html} ".format(fastp_report_path=fastp_report_path, fastp_report_json=fastp_report_json, fastp_report_html=fastp_report_html) + fastp_layout_command_string + " 1> {fastp_log} 2>&1 \n".format(fastp_log=fastp_log)
		fastp_script += fastp_sample_script + "\n\n"
	else:
		##for each_sample in fastq_file_path_Dict:
		pass

	return fastp_script


def build_sortmeRNA_script(study_metadata_Dict, study):
	"""
	sortmerna \\
	--threads {threads} \\
	--workdir $TEMP_PATH \\
	--ref {SORTMERNA_DB}/silva-euk-28s-id98.fasta \\
	--ref {SORTMERNA_DB}/rfam-5.8s-database-id98.fasta \\
	--ref {SORTMERNA_DB}/silva-arc-16s-id95.fasta \\
	--ref {SORTMERNA_DB}/silva-euk-18s-id95.fasta \\
	--ref {SORTMERNA_DB}/rfam-5s-database-id98.fasta \\
	--ref {SORTMERNA_DB}/silva-bac-23s-id98.fasta \\
	--ref {SORTMERNA_DB}/silva-bac-16s-id90.fasta \\
	--ref {SORTMERNA_DB}/silva-arc-23s-id98.fasta \\
	--reads {input.unmap_fastq} \\
	{params.sortmerna_layout_sortmerna_command} \\
	--aligned $TEMP_PATH/${{GENERAL_TAG}} --fastx \\
	"""

	sortmerna_script = "###############################################\n"
	sortmerna_script += "#SORTMERNA\n"

	SORTMERNA_DB = "/usr/local/apps/sortmeRNA/3.0.3/sortmerna-3.0.3/rRNA_databases"

	LAYOUT = study_metadata_Dict[study]["layout"]
	fastq_file_path_Dict = parse_datadir(study_metadata_Dict[study]["input"], study_metadata_Dict[study]["layout"])
		
	for each_sample in fastq_file_path_Dict:
		##
		forward_fastq_path, reverse_fastq_path, unstranded_fastq_name = process_fastq_path(fastq_file_path_Dict[each_sample][0])
		if LAYOUT == "paired":
			#
			sortmerna_layout_command_string = " --reads " + forward_fastq_path + " --reads " + reverse_fastq_path
		else:
			#
			sortmerna_layout_command_string = " --reads " + forward_fastq_path
		#
		sortmerna_report_path = study_metadata_Dict[study]["output"] + "/" + study + "/sortmerna"
		sortmerna_report = sortmerna_report_path + "/" + each_sample
		sortmerna_log = sortmerna_report_path + "/" + each_sample + ".sortmerna.log"
		sortmerna_sample_script = "mkdir -p {sortmerna_report_path}; module load sortmerna; sortmerna --threads $SLURM_CPUS_PER_TASK --log --aligned {sortmerna_report} ".format(sortmerna_report_path=sortmerna_report_path, sortmerna_report=sortmerna_report) + sortmerna_layout_command_string
		sortmerna_sample_script += " --ref {sortmerna_DB}/silva-euk-28s-id98.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/rfam-5.8s-database-id98.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/silva-arc-16s-id95.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/silva-euk-18s-id95.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/rfam-5s-database-id98.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/silva-bac-23s-id98.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/silva-bac-16s-id90.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_sample_script += "--ref {sortmerna_DB}/silva-arc-23s-id98.fasta ".format(sortmerna_DB=SORTMERNA_DB)
		sortmerna_script += sortmerna_sample_script + " 1> {sortmerna_log} 2>&1 \n\n".format(sortmerna_log=sortmerna_log)
	else:
		##for each_sample in fastq_file_path_Dict:
		pass
	return sortmerna_script


def build_multiqc_script(study_metadata_Dict, study):
	"""
	"""
	multiqc_script = "#!/bin/bash\n"
	multiqc_script += ""
	multiqc_script += ""
	multiqc_script += "###############################################\n"
	multiqc_script += "#MULTIQC\n"
	multiqc_script += ""
	multiqc_report_path = study_metadata_Dict[study]["output"] + "/" + study
	multiqc_script += "mkdir -p {multiqc_report_path}\n".format(multiqc_report_path=multiqc_report_path)
	multiqc_script += "module load multiqc\n"
	multiqc_script += "cd {multiqc_report_path}\n".format(multiqc_report_path=multiqc_report_path)
	multiqc_script += "multiqc --force --exclude general_stats --filename {study}_multiqc_report . \n".format(study=study)
		
	return multiqc_script
# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
# LAYOUT
FORWARD_DELIMITER_LIST = ["_1", "_R1"]
REVERSE_DELIMITER_LIST = ["_2", "_R2"]
# ------------------------------------
# ################################### EXEC ####################################### #

samplesheet_Dict = parse_samplesheet_to_dict(sys.argv[1])
study_metadata_Dict = {}
build_study_metadata_dict(samplesheet_Dict, study_metadata_Dict)
########
#


for each_study in study_metadata_Dict:
	##
	print("processing swarm script for " + each_study)
	swarm_script_String = ""
	execution_script_string = "#!/bin/bash\n\n"

	swarm_script_String += build_kraken_script(study_metadata_Dict, each_study)
	swarm_script_String += build_fastqc_script(study_metadata_Dict, each_study)
	swarm_script_String += build_fastq_screen_script(study_metadata_Dict, each_study)
	swarm_script_String += build_fastp_script(study_metadata_Dict, each_study)
	#execution_script_String += build_sortmeRNA_script(study_metadata_Dict, each_study)
	
	swarm_file_name = sys.argv[1].replace(".csv", "_" + each_study + "_script.swarm")
	with open(swarm_file_name, 'w') as execution_swarm:
		execution_swarm.write(swarm_script_String)

	multiqc_string = build_multiqc_script(study_metadata_Dict, each_study)

	multiqc_file_name = sys.argv[1].replace(".csv", "_" + each_study + "_multiqc.sh")
	with open(multiqc_file_name, 'w') as multiqc_sh:
		multiqc_sh.write(multiqc_string)
	
	execution_script_string += "swarm_jobID=$(swarm -g 50 -t 30 --time 04:00:00 --logdir ./logs/ --partition quick --sbatch '--mail-type=FAIL' -f " + swarm_file_name + ")\n"
	
	execution_script_string += "\n\nsbatch --mem=10G --cpus-per-task=10 --partition=quick --time=04:00:00 --dependency afterany:$swarm_jobID " + multiqc_file_name + "\n"

	with open(sys.argv[1].replace(".csv", "_" + each_study + "_execution.sh"), 'w') as execution_script:
		execution_script.write(execution_script_string)
else:
	##
	pass
print("Done!!")
# ################################### FINITO #################################### #