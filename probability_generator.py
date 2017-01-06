#/usr/bin/python
#Author: Ben Wilson <ben.wilson.87@gmail.com>
#PI: Dmitri Petrov
#Description: A Python script for sampling lineages to generate empirical probabilities

import sys, string, math, time, glob, ast, os, re, numpy
start = time.time()
os.chdir('../Output/') #cd to Output file directory (optional)
datafiles = glob.glob('outputFinal*') #Grab all txt files, make sure you point to proper directory
output_stats = open('fixation_statistics.csv','w') #Open summary output file

#Write file headers
output_stats.write("p_soft"+","+"p_hard"+","+"p_ext"+","+"heterozygosity"+","+"b_w"+","+"d_w"+","+"b_m"+","+"d_m"+","+"wU"+","+"Kmax"+","+"n_simulations"+","+"n_sweeps"+"\n")

for datafile in datafiles: #For each parameter range file containing lineage distributions
	#Split filename values, assign to header variables
	values = re.split("_", datafile)
	Kmax, b_w, d_w, b_m, d_m, wU = values[2], values[4], values[6], values[8], values[10], values[12]

	#Counters for each parameter range
	number_simulations = 0
	number_sweeps = 0
	no_sweeps = 0
	heterozygosities = []
	hard_sweeps = 0
	soft_sweeps = 0

	#Open each file, iterate through each lineage list to characterize sweeps
	filehandle = open(datafile)
	for lineage_dist in filehandle:
		number_simulations+=1
		lineage_dist= ast.literal_eval(lineage_dist)
		if(lineage_dist[0] == 0 and len(lineage_dist) == 1):
			no_sweeps+=1
		else:
			mutant_lineage_dist= lineage_dist[1:]
			if(len(mutant_lineage_dist) == 1):
				hard_sweeps+=1
			else:
				soft_sweeps+=1

			#Calculate heterozygosity
			homozygosity = 0.0
			sample_size = 2
			sum_mutant_lineage_dist= 0.0
			for mutant_frequency in mutant_lineage_dist:
				sum_mutant_lineage_dist+=mutant_frequency
			if sum_mutant_lineage_dist> 0 :
				for i in range(len(mutant_lineage_dist)):
					mutant_lineage_dist[i] /= sum_mutant_lineage_dist

			for mutant_frequency in mutant_lineage_dist:
				homozygosity += (mutant_frequency)**sample_size
			heterozygosities.append(1.0-homozygosity)


	heterozygosity = numpy.mean(heterozygosities) if len(heterozygosities) > 0 else 'NA'
	p_soft = soft_sweeps/float(number_simulations)
	p_hard = hard_sweeps/float(number_simulations)
	p_ext = no_sweeps/float(number_simulations) #Divide by number of lines in file
	number_sweeps = number_simulations - no_sweeps
	output_stats.write(str(p_soft)+","+str(p_hard)+","+str(p_ext)+","+str(heterozygosity)+","+str(b_w)+","+str(d_w)+","+str(b_m)+","+str(d_m)+","+str(wU)+","+str(Kmax)+","+str(number_simulations)+","+str(number_sweeps)+"\n")
output_stats.close()
elapsed = time.time() - start
print "Complete: "+str(elapsed)
