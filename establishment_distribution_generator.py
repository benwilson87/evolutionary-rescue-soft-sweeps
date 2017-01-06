#/usr/bin/python
#Author: Ben Wilson <ben.wilson.87@gmail.com>
#PI: Dmitri Petrov
#Description: A Python sampling script for sampling lineages 
import sys, string, math, time, glob, ast, os, re, numpy
start = time.time()
os.chdir('../Output/') #cd to Output file directory
datafiles = glob.glob('outputEstablishment*') #Grab all txt files
output_stats = open('establishment_distribution.csv','w') #Open summary output file

#Write file headers

output_stats.write("t_1"+","+"t_2"+","+"b_w"+","+"d_w"+","+"b_m"+","+"d_m"+","+"wU"+","+"Kmax"+"\n")
		
for datafile in datafiles: #For each parameter range file containing establishment distribution
	#Split filename values, assign to header variables
	values = re.split("_", datafile)	
	Kmax, b_w, d_w, b_m, d_m, wU = values[2], values[4], values[6], values[8], values[10], values[12]	

	#Open each file, iterate through each SFS to characterize sweeps
	filehandle = open(datafile)	
	for establishments in filehandle:
		establishment_list = ast.literal_eval(establishments)
		if(len(establishment_list) == 0):
			t_1 = 'NA'
			t_2 = 'NA'
			output_stats.write(str(t_1)+","+str(t_2)+","+str(b_w)+","+str(d_w)+","+str(b_m)+","+str(d_m)+","+str(wU)+","+str(Kmax)+"\n")
		elif(len(establishment_list) == 1):
			t_1 = establishment_list[0]
			t_2 = 'NA'
			output_stats.write(str(t_1)+","+str(t_2)+","+str(b_w)+","+str(d_w)+","+str(b_m)+","+str(d_m)+","+str(wU)+","+str(Kmax)+"\n")
			
		else:				
			t_1 = establishment_list[0]
			t_2 = establishment_list[1] 
			output_stats.write(str(t_1)+","+str(t_2)+","+str(b_w)+","+str(d_w)+","+str(b_m)+","+str(d_m)+","+str(wU)+","+str(Kmax)+"\n")
output_stats.close()
elapsed = time.time() - start
print "Complete: "+str(elapsed)
