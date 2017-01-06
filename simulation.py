#!/usr/bin/python
#Author: Ben Wilson <benwilson87@gmail.com>
#PI: Dmitri Petrov
#Description: Python functions for birth-death simulations of evolutionary rescue
import sys, string, numpy, math, time, random, itertools
#from multiprocessing import Pool (optional, for parallel processes)

#-----CLASS DEFINITION------
class Population: #Class definition for population object

	#-----------METHODS---------
	def __init__(self, initial_population_size_int, maximum_carrying_capacity_int, wildtype_birth_rate_float, wildtype_death_rate_float, mutant_birth_rate_float, mutant_death_rate_float, mutation_rate_float, tracking_binary, simulation_id_int): #Initialization

		#Static attributes
		self.w_0 = initial_population_size_int #Assign initial population size to wildtype population
		self.K_max = maximum_carrying_capacity_int #Assign maximum carrying capacity
		self.b_mutant = mutant_birth_rate_float #Assign genotype-intrinsic birth rate of mutant
		self.d_mutant = mutant_death_rate_float #Assign genotype-intrinsic death rate of mutant
		self.b_wildtype = wildtype_birth_rate_float #Assign genotype-intrinsic birth rate of wildtype
		self.d_wildtype = wildtype_death_rate_float #Assign genotype-intrinsic death rate of wildtype
		self.mu = mutation_rate_float #Assign mutation rate toward resistant mutant
		self.tracking = tracking_binary #Assign binary value to allow for tracking of lineage trajectories
		self.simulation_id = simulation_id_int #Assign unique id to each simulation

		#Dynamic attributes (I initially built the population object to contain lists for lineages and event rates in order to check the validity of the birth-death algorithm;
		#however, this could be done much more efficiently as a dictionary without the drawbacks of dynamic memory allocation. Feel free to change this for the sake of speed.
		self.population_size = self.w_0 #Initialize population size
		self.wildtype_size = self.w_0 #Initialize variable for wildtype population size tracking
		self.mutant_size = 0 #Initialize variable for mutant population size tracking
		self.lineages = [self.w_0] #Create list to separate individual lineages, initialize with wildtype
		self.event_time = 0.0 #Create float to store time of each event (mutation, birth, or death)
		self.event_rates = [self.mu*self.w_0, self.b_wildtype*self.w_0, self.d_wildtype*self.w_0] #Create and initialize list to store Poisson-distributed rates for each possible event: mutation is always first position, births/deaths are always paired starting with wildtype lineage

		#Tracking attributes for data analysis
		self.number_of_mutation_events = 0 #Create variable to store number of mutations that occur
		self.mutant_lineage_births = [] #Create empty list to store mutant birth times
		self.crossover_time = 0 #Create variable to hold time when wildtype and mutant subpopulations are same size
		self.halfway_fixed_time = 0 #Create variable to hold time when mutant subpopulation is halfway to carrying capacity
		self.fixation_time = 0 #Create variable to hold time when mutant subpopulation fixes at original carrying capacity
		self.lineages_half_sweep = [] #Empty list to store halfway state of population during sweep
		self.wildtype_filehandle = "outputWildtype_Kmax_"+str(self.K_max)+"_bw_"+str(self.b_wildtype)+"_dw_"+str(self.d_wildtype)+"_bm_"+str(self.b_mutant)+"_dm_"+str(self.d_mutant)+"_wU_"+str(self.w_0*self.mu) #Filehandle to track wildtype lineage size through time
		self.mutant_filehandle = "outputMutant_Kmax_"+str(self.K_max)+"_bw_"+str(self.b_wildtype)+"_dw_"+str(self.d_wildtype)+"_bm_"+str(self.b_mutant)+"_dm_"+str(self.d_mutant)+"_wU_"+str(self.w_0*self.mu) #Filehandle to track mutant lineage size through time


	def set_mutation_rate(self): #Set mutation rate (time-varying)
		mutation_rate = self.mu*self.lineages[0] #Assign product of beneficial mutation rate and wildtype population size
		return mutation_rate #Return mutation rate

	def set_birth_rate_wildtype(self): #Set birth rate for wildtype (non-resistant)
		birth_rate_wildtype = self.b_wildtype*self.lineages[0] #Assign product of birth rate and lineage size
		return birth_rate_wildtype #Return wildtype birth rate

	def set_death_rate_wildtype(self): #Set wildtype death rate
		death_rate_wildtype = self.d_wildtype*self.lineages[0] #Assign product of death rate and lineage size
		return death_rate_wildtype #Return wildtype death rate

	def set_birth_rate_mutant(self, index_int, population_size_int): #Set birth rate for mutant (resistant)
		birth_rate_mutant = self.b_mutant*self.lineages[index_int]*(1.0 - population_size_int/float(self.K_max)) #Assign product of birth rate and lineage size
		return birth_rate_mutant #Return mutant birth rate

	def set_death_rate_mutant(self, index_int): #Set death rate for mutant
		death_rate_mutant = self.d_mutant*self.lineages[index_int] #Assign product of death rate and lineage size, scaled by minimum population density
		return death_rate_mutant #Return mutant death rate

	def set_time_of_next_event(self): #Use competing Poisson processes to sample time of next event
		rate_sum = sum(self.event_rates) #Sum all the rates of all events that could occur
		dt = random.expovariate(rate_sum) #Draw the waiting time until an event occurs from an exponential distribution
		self.event_time += dt #Add the waiting time to the current (previous event) time

	def event_sampling(self): #Weighted sampling of Poisson-distributed events based on relative rates
		rnd = random.random()*sum(self.event_rates) #Take a random floating point number in the length of the sum of the rates
		event_rate_index = 0 #Initialize the index of the rate corresponding to the event (mutation default)
		for i, rate in enumerate(self.event_rates): #For each rate and index pair in the list of event rates
			rnd -= rate #Subtract rate from random floating point, in order of increasing index order
			if rnd < 0: #If result is negative, floating point corresponds to that rate index,
				event_rate_index = i #assign index of corresponding event
				break
		return event_rate_index #Return index of event with the corresponding event

	def mutation_event(self):
		self.lineages[0]-=1 #kill a wildtype individual,
		self.lineages.append(1) #append a mutant birth to the lineage list
		self.wildtype_size-=1 #update wildtype population size
		self.mutant_size+=1 #update mutant population size
		self.number_of_mutation_events+=1 #update number of mutations variable
		if(self.tracking and self.simulation_id%100 == 1):
			output_wildtype = open(self.wildtype_filehandle,"a")
			output_wildtype.write(str(self.event_time)+","+str(self.wildtype_size)+","+str(self.simulation_id)+"\n")
			output_wildtype.close()
			output_mutant = open(self.mutant_filehandle,"a")
			output_mutant.write(str(self.event_time)+","+str(self.mutant_size)+","+str(self.simulation_id)+"\n")
			output_mutant.close()
		self.mutant_lineage_births.append(self.event_time) #Add time that mutant lineage emerges
		self.event_rates[0] = self.set_mutation_rate() #Update mutation rate
		self.event_rates[1] = self.set_birth_rate_wildtype() #Update birth rate of wildtype lineage for new lineage size
		self.event_rates[2] = self.set_death_rate_wildtype() #Update death rate of wildtype lineage for new lineage size
		self.event_rates.extend((self.set_birth_rate_mutant(-1,self.population_size), self.set_death_rate_mutant(-1))) #append two new entries to event rate list (mutant birth and death rates)

	def wildtype_birth_death(self, event_int):
		self.lineages[0]+=event_int #update corresponding wildtype lineage size
		self.population_size+=event_int #update population size
		self.wildtype_size+=event_int #update wildtype population size
		if(self.tracking and self.simulation_id%100 == 1):
			output_wildtype = open(self.wildtype_filehandle,"a")
			output_wildtype.write(str(self.event_time)+","+str(self.wildtype_size)+","+str(self.simulation_id)+"\n")
			output_wildtype.close()
		self.update_density_dependent_rates() #update all density-dependent mutant birth rates
		self.event_rates[0] = self.set_mutation_rate() #update mutation rate
		self.event_rates[1] = self.set_birth_rate_wildtype() #update wildtype birth rate
		self.event_rates[2] = self.set_death_rate_wildtype() #update wildtype death rate


	def mutant_birth_death(self, lineage_index, event_int):
		self.lineages[lineage_index]+=event_int #update corresponding mutant lineage size
		self.population_size+=event_int #update population size
		self.mutant_size+=event_int #update mutant population size
		if(self.tracking and self.simulation_id%100 == 1):
			output_mutant = open(self.mutant_filehandle,"a")
			output_mutant.write(str(self.event_time)+","+str(self.mutant_size)+","+str(self.simulation_id)+"\n")
			output_mutant.close()
		self.update_density_dependent_rates() #update all density-dependent mutant birth rates
		self.event_rates[(2*lineage_index + 2)] = self.set_death_rate_mutant(lineage_index) #update mutant death rate

	def update_density_dependent_rates(self):
		rate_index = 3 #Start with index of first mutant birth rate
		while(rate_index < len(self.event_rates)): #While mutant birth rates exist,
			self.event_rates[rate_index] = self.set_birth_rate_mutant((rate_index-1)/2,self.population_size) #update with new population size
			rate_index+=2

	def population_update(self, event_index_int): #Update the individual lineages after an event has occurred
		lineage_index  = (event_index_int-1)/2 #Assign index for lineage to be updated
		if(event_index_int == 0): #If a mutation event has occurred,
			self.mutation_event() #call mutation event method

		elif(event_index_int == 1 or event_index_int ==2): #Else if wildtype birth/death event has occurred,
			event = 1 if (event_index_int == 1) else -1 #determine if birth or death,
			self.wildtype_birth_death(event) #call wildtype birth/death method

		else: #Else if mutant birth/death event has occurred,
			event = 1 if (event_index_int%2 == 1) else -1 #determine if birth or death,
			self.mutant_birth_death(lineage_index, event) #call mutant birth/death method

		self.remove_extinct_lineages() #Remove any mutant lineages that have gone extinct
		self.population_tracking_update() #Update various measures that track population dynamics

	def population_tracking_update(self):
		if(self.population_size - self.lineages[0] >= self.lineages[0] and self.crossover_time == 0): #If mutant subpopulation is as large or larger than wildtype subpopulation,
			self.crossover_time = float(self.event_time) #update crossover time
		if(self.population_size - self.lineages[0] >= (0.5*self.K_max)*(1.0-self.d_mutant/self.b_mutant) and self.halfway_fixed_time == 0): #If mutant subpopulation is at least half of its equilibrium population size,
			self.halfway_fixed_time, self.lineages_half_sweep =  float(self.event_time), list(self.lineages) #update halfway time and store population state at halfway point
		if(self.population_size - self.lineages[0] >= (0.99*self.K_max)*(1.0-self.d_mutant/self.b_mutant) and self.fixation_time == 0): #If mutant subpopulation is fixed to equilibrium size,
			self.fixation_time = float(self.event_time) #update fixation time

	def remove_extinct_lineages(self): #Remove any lineages that go extinct during process
		lineage_index = len(self.lineages)-1 #Beginning with final index
		while(lineage_index): #Garbage collection for dead lineages that do not need to be tracked ad infinitum
			if(self.lineages[lineage_index] == 0): #If mutant lineage has gone extinct,
				del self.lineages[lineage_index] #delete corresponding lineage
				del self.mutant_lineage_births[lineage_index - 1] #delete corresponding birth time for extinct mutant
				del self.event_rates[(2*lineage_index + 2)] #delete corresponding death rate
				del self.event_rates[(2*lineage_index + 1)] #delete corresponding birth rate
			lineage_index-=1

	def birth_death_process(self): #Dynamic process that simulates birth-death model
		self.set_time_of_next_event() #Sample time of first event, append time
		event = self.event_sampling() #Assign event
		self.population_update(event) #Update population dynamics with corresponding event

#----------Main method----------
def main(parameter_list):

	#Parameters to be set manually through STDIN
	K_max = int(parameter_list[0]) #Assign maximum carrying capacity
	b_wildtype = float(parameter_list[1]) #Assign wildtype birth rate
	d_wildtype = float(parameter_list[2]) #Assign wildtype death rate
	b_mutant = float(parameter_list[3]) #Assign mutant birth rate
	d_mutant = float(parameter_list[4]) #Assign mutant death rate
	mu = float(parameter_list[5]) #Assign mutation rate toward resistant mutant
	w_0 = int(parameter_list[6]) #Asssign initial population size to wildtype population
	number_simulations = int(parameter_list[7]) #Assign integer for number of simulations for each set of parameters

	#Simulate same population for certain number of times
	tracking = 0 #Assign binary value for lineage tracking
	while(number_simulations):
		numpy.random.seed() #Seed random number generator (optional: include a process ID or the like to make sure parallel processes are independent)
		population = Population(w_0, K_max, b_wildtype, d_wildtype, b_mutant, d_mutant, mu, tracking, number_simulations) #Create a population, evolve until extinction or adaptation occurs

		while(population.population_size > 0): #While population is not extinct,
			population.birth_death_process() #run birth death process
			if(population.fixation_time != 0): #If mutant population sweeps,
				break #stop simulation

		#Output data to filehandles
		number_of_mutations_filehandle = "outputNumberMutations_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		established_mutant_lineages_filehandle = "outputEstablishment_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		crossover_filehandle = "outputCrossover_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		half_sweep_filehandle = "outputHalfway_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		fixation_filehandle = "outputFixation_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		halfway_lineage_filehandle = "outputHalfwayLineages_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)
		final_lineage_filehandle = "outputFinalLineages_Kmax_"+str(population.K_max)+"_bw_"+str(population.b_wildtype)+"_dw_"+str(population.d_wildtype)+"_bm_"+str(population.b_mutant)+"_dm_"+str(population.d_mutant)+"_wU_"+str(population.w_0*population.mu)


		output_mutations = open(number_of_mutations_filehandle,"a")
		output_mutations.write(str(population.number_of_mutation_events)+"\n")
		output_mutations.close()

		output_establishment = open(established_mutant_lineages_filehandle,"a")
		output_establishment.write(str(population.mutant_lineage_births)+"\n")
		output_establishment.close()

		output_crossover = open(crossover_filehandle, "a")
		output_crossover.write(str(population.crossover_time)+"\n")
		output_crossover.close()

		output_halfway = open(half_sweep_filehandle, "a")
		output_halfway.write(str(population.halfway_fixed_time)+"\n")
		output_halfway.close()

		output_fixation = open(fixation_filehandle, "a")
		output_fixation.write(str(population.fixation_time)+"\n")
		output_fixation.close()

		output_halfway_lineages = open(halfway_lineage_filehandle, "a")
		output_halfway_lineages.write(str(population.lineages_half_sweep)+"\n")
		output_halfway_lineages.close()

		output_final_lineages = open(final_lineage_filehandle, "a")
		output_final_lineages.write(str(population.lineages)+"\n")
		output_final_lineages.close()

		number_simulations-=1
	return 0



#----------Run main method---------
#Create lists of parameters to explore
K_max = 110000
b_wildtype = 0.0
d_wildtype = 1.0
b_mutant = 1.1
d_mutant = 1.0
mu = 0.001
w_0 = 10**4
num_simulations = 100
parameters = [K_max,b_wildtype,d_wildtype,b_mutant,d_mutant,mu,w_0,num_simulations]
start = time.time()

if __name__ == '__main__':
	main(parameters)
elapsed = time.time() - start
print str(elapsed)
