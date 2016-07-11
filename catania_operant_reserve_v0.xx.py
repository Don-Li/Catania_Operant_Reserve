# -*- coding: utf-8 -*-
"""
Current patch notes:
v.04    11/7/16
    Changed functions:
        herrnsteins_hyperbola()
            Now takes the base operant level as a specified global variable.
"""
import numpy as np
import os
from scipy.optimize import curve_fit

def reserve_replenishment(delay_vector,increment_max,function_type):
    if function_type == "exponential":
        reserve_replenishment = np.sum(increment_max*np.exp((-delay_vector*0.001811)+1))
    if function_type == "reciprocal":
        reserve_replenishment = np.sum(increment_max/delay_vector)
    if function_type == "hyperbolic":
        reserve_replenishment = np.sum(2*increment_max/(1+delay_vector))
    if function_type == "linear":
        reserve_replenishment = np.sum(-delay_vector + increment_max)
    return reserve_replenishment

def variable_interval_schedule(VI_interval, time):
    next_rft_available = int((np.random.exponential(VI_interval))) + time
    return next_rft_available

def behaviour_emission(reserve_now, reserve_max):
    behaviour = False
    if reserve_now/reserve_max > np.random.random(1):
        behaviour = True
    return behaviour

def recording_matrices(conditions_matrix, levels_per_factor, number_of_sessions):
    increment_max_conditions = levels_per_factor[1]
    rft_rate_levels = levels_per_factor[2]
    decrement_value_conditions = levels_per_factor[0]
    total_entries = number_of_sessions * time_to_run
    record_matrix_0 = np.empty(shape = [total_entries, len(conditions_matrix)], dtype = float) #This is for experimental details, e.g. increment_max
    record_matrix_1 = np.zeros(shape = [total_entries, 5], dtype = int) #This is for Time, Cum_rfts, Cum_reserve1, Rft, Reserve1_output
    record_matrix_2 = np.empty(shape = [total_entries, 1] , dtype = float) #This is for the current reserve_value
    record_matrix_0[:,0] = np.repeat(conditions_matrix[0], increment_max_conditions*rft_rate_levels*time_to_run)
    record_matrix_0[:,1] = np.tile(np.repeat(conditions_matrix[1], rft_rate_levels*time_to_run), decrement_value_conditions)
    record_matrix_0[:,2] = np.tile(np.repeat(conditions_matrix[2], time_to_run), increment_max_conditions*decrement_value_conditions)
    record_matrix_1[:,0] = np.tile(range(0,time_to_run), number_of_sessions)
    return (record_matrix_0,record_matrix_1, record_matrix_2)

def schedule_block(reserve_max, increment_max, decrement_value, operant_level, initial_reserve_level, time_to_run, random_interval, function_type, record_matrix_1, record_matrix_2, session):
    time = 0
    reinforcement_arranged = False
    last_reinforcer_time_1 = 0
    last_reinforcer_time_2 = 0
    next_reinforcer_time = None
    reserve_now = initial_reserve_level

    while time < time_to_run:
        record_matrix_2[time+(session*time_to_run)] = reserve_now
        emitted_behaviour = behaviour_emission(reserve_now, reserve_max)
        if reinforcement_arranged == False:
            next_reinforcer_time = variable_interval_schedule(random_interval, time)
            reinforcement_arranged = True
        if emitted_behaviour == True:
            record_matrix_1[time+(session*time_to_run),4] = 1 #Reserve1_output column
            reserve_now -= decrement_value
            if reserve_now <= 0:
                reserve_now = 0
            if time >= next_reinforcer_time:
                record_matrix_1[time+(session*time_to_run),3] = 1 #Rft column
                last_reinforcer_time_2 = last_reinforcer_time_1
                last_reinforcer_time_1 = time+1
                delay_vector = np.where(record_matrix_1[last_reinforcer_time_2+(session*time_to_run):last_reinforcer_time_1+(session*time_to_run),4][::-1])[0]+1
                reserve_now += reserve_replenishment(delay_vector, increment_max, function_type)
                if reserve_now >= reserve_max:
                    reserve_now = reserve_max
                reinforcement_arranged = False
            if time < next_reinforcer_time:
                record_matrix_1[time+(session*time_to_run),3] = 0 #Rft column
        if emitted_behaviour == False:
            record_matrix_1[time+(session*time_to_run),4] = 0
            record_matrix_1[time+(session*time_to_run),3] = 0
        record_matrix_1[time+(session*time_to_run),1:3] = np.array([np.sum(record_matrix_1[(session*time_to_run):((session)*time_to_run)+time+1,3]),np.sum(record_matrix_1[(session*time_to_run):((session)*time_to_run)+time+1,4])])
        time += 1
    return None

def experiment_design_values(variable_interval_conditions, decrement_value_conditions, increment_max_conditions):
    conditions = [variable_interval_conditions, decrement_value_conditions, increment_max_conditions]
    condition_numbers = np.empty(3, dtype = int) #3 is the number of factors
    for i in range(0,3):
        condition_numbers[i] = len(conditions[i])
    number_of_sessions = np.product(condition_numbers)
    conditions_matrix = [decrement_value_conditions, increment_max_conditions, variable_interval_conditions]
    levels_per_factor = np.empty(3, dtype = int)
    for i in range(0,3):
        levels_per_factor[i] = len(conditions_matrix[i])
    return (conditions_matrix, levels_per_factor, number_of_sessions)

def run(variable_interval_conditions, decrement_value_conditions, increment_max_conditions, record_matrix, number_of_sessions):
    session = 0
    for d_condition in decrement_value_conditions:
        decrement_value = d_condition
        for i_condition in increment_max_conditions:
            increment_max = i_condition
            for v_condition in variable_interval_conditions:
                variable_interval = v_condition
                schedule_block(reserve_max, increment_max, decrement_value,operant_level, initial_reserve_level, time_to_run, variable_interval, function_type, record_matrix[1], record_matrix[2], session)
                session += 1
    return None

def calculate_statistics(conditions_matrix, levels_per_factor, number_of_sessions, record_matrix_1, bin_size):
    increment_max_levels = levels_per_factor[1]
    rft_rate_levels = levels_per_factor[2]
    decrement_value_levels = levels_per_factor[0]
    statistics_matrix[:,0] = np.repeat(conditions_matrix[0], increment_max_levels*rft_rate_levels)
    statistics_matrix[:,1] = np.tile(np.repeat(conditions_matrix[1], rft_rate_levels), decrement_value_levels)
    statistics_matrix[:,2] = np.tile(conditions_matrix[2], increment_max_levels*decrement_value_levels)
    within_condition_indices_end = np.array(range(1,int(time_to_run/bin_size*number_of_sessions)+1))*(bin_size)
    within_condition_indices_start = within_condition_indices_end-bin_size
    session_counts = np.empty((int(time_to_run/bin_size*number_of_sessions),2))
    for i in range(0, number_of_sessions):
        for j in range(0,len(within_condition_indices_end)):
            session_counts[j] = np.sum(record_matrix_1[within_condition_indices_start[j]:within_condition_indices_end[j]][:,3:5],axis = 0)
        statistics_matrix[i,3:5] = np.average(session_counts[i*int(time_to_run/bin_size):i*int(time_to_run/bin_size)+int(time_to_run/bin_size)], axis = 0)
    return None

def herrnsteins_hyperbola(reinforcement_rate, k ,re):
    predicted_response_rate =  reinforcement_rate * k / (reinforcement_rate + re) + base_operant_level
    return predicted_response_rate

def fit_hyperbola(statistics_matrix, fitted_hyperbola_matrix, levels_per_factor,conditions_matrix):
    decrement_level_levels = levels_per_factor[0]
    increment_max_levels = levels_per_factor[1]
    initial_guesses = np.array([0,0])
    fitted_hyperbola_matrix[:,0] = np.repeat(conditions_matrix[0], increment_max_levels)
    fitted_hyperbola_matrix[:,1] = np.tile(conditions_matrix[1], decrement_level_levels)
    k = 0
    for i in conditions_matrix[0]: #decrement levels
        for j in conditions_matrix[1]: #increment_max levels
            response_rate = statistics_matrix[(statistics_matrix[:,0]==i) & (statistics_matrix[:,1]==j),4]
            reinforcement_rate = statistics_matrix[(statistics_matrix[:,0]==i) & (statistics_matrix[:,1]==j),3]
            fitted_hyperbola_matrix[k, 2:] =  curve_fit(f = herrnsteins_hyperbola, xdata = reinforcement_rate, ydata = response_rate, p0 = initial_guesses, bounds = (0,np.inf), method = "dogbox")[0]
            k += 1
    return None

def output_to_csv(record_matrix, working_directory, variable_interval_conditions, decrement_value_conditions, increment_max_conditions):
    os.chdir(working_directory)
    VI_string = 'VI '+''.join(str(variable_interval_conditions)).strip('[]') + "_"
    decrement_value_condition_string = "D "+''.join(str(decrement_value_conditions)).strip('[]') + "_"
    increment_max_condition_string = "I"+''.join(str(increment_max_conditions)).strip('[]') + "_"
    file_name = "COR RAW " + VI_string + decrement_value_condition_string + increment_max_condition_string + ".csv"
    combined_record_matrix = np.concatenate((record_matrices), axis = 1)
    file_headers = ["decrement_value", "increment_max", "VI", "Time", "Cum_rfts", "Cum_target", "Rft", "Reserve1_output","Reserve1_now"]
    file_header_string = ','.join(file_headers)
    np.savetxt(file_name, combined_record_matrix, delimiter = ",", header = file_header_string)
    
    file_name_statistics_matrix = "COR STATISTICS " + VI_string + decrement_value_condition_string + increment_max_condition_string + ".csv"
    file_headers_statistics_matrix = ["decrement_value number", "increment_max", "Arranged rft rate", "Observed rft rate", "Response rate"]
    file_header_string_statistics_matrix = ','.join(file_headers_statistics_matrix)
    np.savetxt(file_name_statistics_matrix, statistics_matrix, delimiter = ",", header = file_header_string_statistics_matrix)
    
    file_name_hyperbola_matrix = "COR HYPERBOLA " + VI_string + decrement_value_condition_string + increment_max_condition_string + ".csv"
    file_headers_hyperbola_matrix = ["decrement_value number", "increment_max", "k", "re"]
    file_header_string_hyperbola_matrix = ','.join(file_headers_hyperbola_matrix)
    np.savetxt(file_name_hyperbola_matrix, fitted_hyperbola_matrix, delimiter = ",", header = file_header_string_hyperbola_matrix)
    return None

#Global variables
reserve_max = 100000
initial_reserve_level = 75000
time_to_run = 25000
function_type = "exponential"
bin_size = 500
base_operant_level = 0

#parameters
increment_max_conditions = np.array([.03])*reserve_max
decrement_value_conditions = np.array([.01])*reserve_max
operant_level = np.array([0])
variable_interval_conditions = np.array([1,2,3,5,8,10,18,25,68,112,200])
experiment_design_value = experiment_design_values(variable_interval_conditions, decrement_value_conditions, increment_max_conditions)
conditions_matrix = experiment_design_value[0]
levels_per_factor = experiment_design_value[1]
number_of_sessions = experiment_design_value[2]

#initialise structures
record_matrices = recording_matrices(conditions_matrix, levels_per_factor, number_of_sessions)
statistics_matrix = np.empty([number_of_sessions,len(conditions_matrix)+2])
fitted_hyperbola_matrix = np.zeros(shape = (np.product(levels_per_factor[0:2]),4), dtype = float)

working_directory = "D:\Dropbox\Dropbox\Don\Catania Operant Reserve wd"
run(variable_interval_conditions, decrement_value_conditions, increment_max_conditions, record_matrices, number_of_sessions)
calculate_statistics(conditions_matrix, levels_per_factor, number_of_sessions, record_matrices[1], bin_size)
fit_hyperbola(statistics_matrix, fitted_hyperbola_matrix, levels_per_factor, conditions_matrix)
output_to_csv(record_matrices, working_directory, variable_interval_conditions, decrement_value_conditions, increment_max_conditions)

