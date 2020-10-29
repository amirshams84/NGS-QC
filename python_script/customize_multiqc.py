# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: October-25-2020
# Email: amir.shams84@gmail.com
# Project: NGS-QC
# Aim: python script to customize generated multiqc html file, to remove inaccurate metrics
# python cutomize_multiqc.py $PATH/NGS_test_multiqc_report.html 
# ################################### IMPORT ##################################### #



import os
import sys
import re
from bs4 import BeautifulSoup
# ################################### FUNCTIONS ################################## #


def remove_fastqc_sequence_counts(multiqc_html_file_Path):
	"""
	"""
	html_File = open(multiqc_html_file_Path, "r")
	html_soup_Object = BeautifulSoup(html_File, 'html.parser')
	
	for each_div in html_soup_Object.find_all("div", class_="mqc-section-plot"):
		##
		try:
			target_level_1 = each_div.find("div", class_="mqc_hcplot_plotgroup").find("div", class_="btn-group hc_switch_group")
		except AttributeError:
			continue
		
		try:
			target_level_2 = target_level_1.find("button", class_=["btn", "btn-default" "btn-sm"]).attrs["data-target"]
		except AttributeError:
			continue
		#print(target_level_2)
		if "fastqc_sequence_counts" in target_level_2:
			#
			#print(each_div)
		
			new_tag = html_soup_Object.new_tag('h5')
			new_tag.string = "-- This plot removed due to the inaccuracy of the calculation method --"
			each_div.replace_with(new_tag)

			#print("####################")
		else:
			pass
	else:
		##
		pass
	html_File.close()
	with open(multiqc_html_file_Path, "w", encoding='utf-8') as file:
		file.write(str(html_soup_Object.html))

	return True

# ################################### GLOBAL VARIABLES ########################### #
# ################################### EXECUTION ################################## #


remove_fastqc_sequence_counts(sys.argv[1])

# ################################### FINITO ##################################### #