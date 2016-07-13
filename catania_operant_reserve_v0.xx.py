# -*- coding: utf-8 -*-
"""
Current patch notes:
v.05    13/7/16
    Added functions:
        cartesian_product()
            Used for generating combinations of factors. Not yet implemented in any functions.
            This function should go inside: recording_matrices(), calculate_statistics() and fit_hyperbola()
        rate_multiplied_gamma_cdf()
            Calculates values for a rate multiplied gamma cumulative density function
        fit_gamma_cdf()
            Fits response rates to observed rft rates with the rate_multiplied_gamma_cdf()
    Changed functions:
        reserve_replenishment()
            Augmented the functions with scaling values.
        recording_matrices()
            Now uses the cartesian_product() function to produce the experiment details matrix.
            The matrix previously used for the probability of responding (i,e, reserve now) merged into behavioural outcomes matrix.
        calculate_statistics()
            filling in the experimental details now done by cartesian_product().
            Now also gives standard deviations for the rates.
        fit_hyperbola()
            fitting in the experimental details now done by cartesian_product()
            calculates r squared and puts it into the fitted hyperbola matrix
            also calculates a vector of residuals from fitting the hyperbola and puts them into statistics_matrix
        output_to_csv()
            writing of file names now done in a loop
            the files are updated to reflect the additional information that is added from other functions.
            Now also writes a file for the fit_gamma_cdf()
            We can now choose if we want to have the raw file.
        schedule_block()
            Now updates the record matrices with the new additional information
            Now takes the replenishment scaling value as an argument
        run()
            Now iterates through the design matrix rather than having nested loops.
        reserve_replenishment()
            Now sums over positive weighting values when linear function is used.
    Removed functions:
        experiment_design_values()
"""
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.special import gdtr

def cartesian_product(conditions_array, number_of_sessions, levels_per_factor, entries_per_session = False):
    cartesian_substrate = conditions_array[:]
    if entries_per_session != False:
        cartesian_substrate.append(list(range(0, entries_per_session)))
    number_of_rows = np.prod([len(condition) for condition in cartesian_substrate])
    cartesian_product = np.zeros(shape = (number_of_rows, len(cartesian_substrate)))
    rep_vector = [1]+[len(x) for x in cartesian_substrate]+[1]
    for j in range(0, len(cartesian_substrate)):
        cartesian_product[:,j] = np.tile(np.repeat(cartesian_substrate[j],np.prod(rep_vector[j+2:])), np.prod(rep_vector[:(j+1)]))
    return cartesian_product

def reserve_replenishment(delay_vector,increment_max,function_type,scaling_value):
    if function_type == "exponential":
        reserve_replenishment = np.sum(increment_max*np.exp(-(delay_vector-+1)/scaling_value))
    if function_type == "reciprocal":
        reserve_replenishment = np.sum(increment_max/(delay_vector**(1/scaling_value)))
    if function_type == "hyperbolic":
        reserve_replenishment = np.sum(2*increment_max/(1+delay_vector**(1/scaling_value)))
    if function_type == "linear":
        weight_vector = (-delay_vector - scaling_value*(delay_vector-1) + (increment_max+1))
        reserve_replenishment = np.sum(weight_vector[weight_vector>0])
    return reserve_replenishment

def variable_interval_schedule(VI_interval, time):
    next_rft_available = int((np.random.exponential(VI_interval))) + time
    return next_rft_available

def behaviour_emission(reserve_now, reserve_max):
    behaviour = False
    if reserve_now/reserve_max > np.random.random(1):
        behaviour = True
    return behaviour

def recording_matrices(conditions_array, levels_per_factor, number_of_sessions, entries_per_session):
    total_entries = number_of_sessions * time_to_run
    record_matrix_experiment_details = cartesian_product(conditions_array, number_of_sessions, levels_per_factor, entries_per_session) #This is for experimental details, e.g. increment_max
    record_matrix_behavioural_outcomes = np.zeros(shape = [total_entries, 6], dtype = int) #This is for Cum_rfts, Cum_reserve1, Rft, Reserve1_output, reserve_value1,reserve_1_increment
    return (record_matrix_experiment_details,record_matrix_behavioural_outcomes)

def schedule_block(reserve_max, increment_max, decrement_value, initial_reserve_level, time_to_run, random_interval, function_type, record_matrix_behavioural_outcomes, session, scaling_value):
    time = 0
    rft_arranged = False
    last_rft_time_1 = 0
    last_rft_time_2 = 0
    next_rft_time = None
    reserve_now = initial_reserve_level

    while time < time_to_run:
        record_matrix_behavioural_outcomes[time+(session*time_to_run),4] = reserve_now
        emitted_behaviour = behaviour_emission(reserve_now, reserve_max)
        if rft_arranged == False:
            next_rft_time = variable_interval_schedule(random_interval, time)
            rft_arranged = True
        if emitted_behaviour == True:
            record_matrix_behavioural_outcomes[time+(session*time_to_run),3] = 1 #Reserve1_output column
            reserve_now -= decrement_value
            if reserve_now <= 0:
                reserve_now = 0
            if time >= next_rft_time:
                record_matrix_behavioural_outcomes[time+(session*time_to_run),2] = 1 #Rft column
                last_rft_time_2 = last_rft_time_1
                last_rft_time_1 = time+1
                delay_vector = np.where(record_matrix_behavioural_outcomes[last_rft_time_2+(session*time_to_run):last_rft_time_1+(session*time_to_run),3][::-1])[0]+1
                reserve_replenishment_value = reserve_replenishment(delay_vector, increment_max, function_type, scaling_value)
                if (reserve_max-reserve_now) < reserve_replenishment_value:
                    reserve_replenishment_value = reserve_max-reserve_now
                reserve_now += reserve_replenishment_value
                record_matrix_behavioural_outcomes[time+(session*time_to_run),5] = reserve_replenishment_value
                rft_arranged = False
            if time < next_rft_time:
                record_matrix_behavioural_outcomes[time+(session*time_to_run),2] = 0 #Rft column
        if emitted_behaviour == False:
            record_matrix_behavioural_outcomes[time+(session*time_to_run),2] = 0 #rft column
            record_matrix_behavioural_outcomes[time+(session*time_to_run),3] = 0 #reserve1_output column
        record_matrix_behavioural_outcomes[time+(session*time_to_run),0:2] = np.array([np.sum(record_matrix_behavioural_outcomes[(session*time_to_run):((session)*time_to_run)+time+1,2]),np.sum(record_matrix_behavioural_outcomes[(session*time_to_run):((session)*time_to_run)+time+1,3])]) #cum_rft and cum_reserve1 columns
        time += 1
    return None

def run(design_matrix, record_matrix, number_of_sessions):
    session = 0
    for i in range(0,len(design_matrix)):
        decrement_value = design_matrix[i,0]
        increment_max = design_matrix[i,1]
        replenishment_scaling = design_matrix[i,2]
        variable_interval = design_matrix[i,3]
        schedule_block(reserve_max, increment_max, decrement_value, initial_reserve_level, time_to_run, variable_interval, function_type, record_matrix[1], session, replenishment_scaling)
        session += 1
    return None

def calculate_statistics(conditions_array, levels_per_factor, number_of_sessions, record_matrix_behavioural_outcomes, bin_size):
    #columns: decrement, increment, replenishment_scaling, VI, rft rate, response rate, rft rate stdev, response rate stdev, residuals
    within_condition_indices_end = np.array(range(1,int(time_to_run/bin_size*number_of_sessions)+1))*(bin_size)
    within_condition_indices_start = within_condition_indices_end-bin_size
    session_counts = np.empty((int(time_to_run/bin_size*number_of_sessions),2))
    for i in range(0, number_of_sessions):
        for j in range(0,len(within_condition_indices_end)):
            session_counts[j] = np.sum(record_matrix_behavioural_outcomes[within_condition_indices_start[j]:within_condition_indices_end[j]][:,2:4],axis = 0)
        statistics_matrix[i,4:6] = np.average(session_counts[i*int(time_to_run/bin_size):i*int(time_to_run/bin_size)+int(time_to_run/bin_size)], axis = 0)
        statistics_matrix[i,6:8] = np.std(session_counts[i*int(time_to_run/bin_size):i*int(time_to_run/bin_size)+int(time_to_run/bin_size)], axis = 0)
    return None

def herrnsteins_hyperbola(reinforcement_rate, k ,re):
    predicted_response_rate =  reinforcement_rate * k / (reinforcement_rate + re) + base_operant_level
    return predicted_response_rate

def fit_hyperbola(statistics_matrix, fitted_hyperbola_matrix, levels_per_factor):
    initial_guesses = np.array([1,1])
    sessions_per_curve = levels_per_factor[-1] #Use each set of rft rate conditions for one curve.
    for i in range(0, len(fitted_hyperbola_matrix)):
        reinforcement_rates = statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),4]
        response_rates = statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),5]
        fitted_hyperbola_matrix[i, 3:5] =  curve_fit(f = herrnsteins_hyperbola, xdata = reinforcement_rates, ydata = response_rates, p0 = initial_guesses, bounds = (0,np.inf), method = "dogbox")[0]
        statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),8] = response_rates - herrnsteins_hyperbola(reinforcement_rates,fitted_hyperbola_matrix[i,3], fitted_hyperbola_matrix[i,4]) #Put the residuals in the statistics matrix
        sum_of_squares_total = np.sum(np.square(response_rates-np.average(response_rates)))
        sum_of_squares_residual = np.sum(np.square(statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),8]))
        r_squared = 1-sum_of_squares_residual/sum_of_squares_total
        fitted_hyperbola_matrix[i, 5] = r_squared
    return None

def rate_multiplied_gamma_cdf(x,shape, scale, rate_multiplier):
    alpha = scale
    beta = shape
    F_x = gdtr(alpha,beta,x)*rate_multiplier
    return F_x

def fit_gamma_cdf(statistics_matrix, fitted_gamma_cdf_matrix, levels_per_factor):
    initial_guesses = np.array([1,1,1])
    sessions_per_curve = levels_per_factor[-1] #Use each set of rft rate conditions for one curve.
    for i in range(0, len(fitted_hyperbola_matrix)):
        reinforcement_rates = statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),4]
        response_rates = statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),5]
        initial_guesses[2] = max(response_rates)
        fitted_gamma_cdf_matrix[i, 3:6] = curve_fit(f=rate_multiplied_gamma_cdf, xdata = reinforcement_rates, ydata = response_rates, p0 = initial_guesses, method = "lm")[0]
        statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),9] = response_rates - rate_multiplied_gamma_cdf(reinforcement_rates,fitted_gamma_cdf_matrix[i,3], fitted_gamma_cdf_matrix[i,4],fitted_gamma_cdf_matrix[i,5]) #Put the residuals in the statistics matrix
        sum_of_squares_total = np.sum(np.square(response_rates-np.average(response_rates)))
        sum_of_squares_residual = np.sum(np.square(statistics_matrix[sessions_per_curve*i : sessions_per_curve*(i+1),9]))
        r_squared = 1-sum_of_squares_residual/sum_of_squares_total
        fitted_gamma_cdf_matrix[i, 6] = r_squared
    return None

def output_to_csv(record_matrix, working_directory, conditions_array, write_raw): #conditions_array order: decrement, increment, replenishment scaling, VI
    os.chdir(working_directory)
    condition_initials = ["D ","I ","S ","V "] #decrement, increment, replenishment scaling, VI
    string_of_conditions = ""
    conditions_listing = conditions_array[:]
    conditions_listing[2] = []
    for i in range(0,len(conditions_array[2])):
        conditions_listing[2].append(str(conditions_array[2][i])[:6])
    for i in range(0,len(condition_initials)):
        string_of_conditions += condition_initials[i]+''.join(str(conditions_listing[i]).strip('[]')) + "_"
    file_name = "COR RAW " + string_of_conditions + ".csv"
    combined_record_matrix = np.concatenate((record_matrices), axis = 1)
    file_headers = ["decrement_value", "increment_max","replenishment_scaling", "VI", "time", "cum_rfts", "cum_target", "rft", "reserve1_output","reserve1_now", "reserve1_increment"]
    file_header_string = ','.join(file_headers)
    if write_raw == True:
        np.savetxt(file_name, combined_record_matrix, delimiter = ",", header = file_header_string)

    file_name_statistics_matrix = "COR STATISTICS " + string_of_conditions + ".csv"
    file_headers_statistics_matrix = ["decrement_value number", "increment_max", "replenishment_scaling", "VI", "Observed rft rate", "Response rate", "rft rate stdev", "response rate stdev","hyperbola residuals", "gamma residuals"]
    file_header_string_statistics_matrix = ','.join(file_headers_statistics_matrix)
    np.savetxt(file_name_statistics_matrix, statistics_matrix, delimiter = ",", header = file_header_string_statistics_matrix)

    file_name_hyperbola_matrix = "COR HYPERBOLA FITS " + string_of_conditions + ".csv"
    file_headers_hyperbola_matrix = ["decrement_value number", "increment_max", "replenishment_scaling", "k", "re", "hyperbola_r_squared"]
    file_header_string_hyperbola_matrix = ','.join(file_headers_hyperbola_matrix)
    np.savetxt(file_name_hyperbola_matrix, fitted_hyperbola_matrix, delimiter = ",", header = file_header_string_hyperbola_matrix)

    file_name_gamma_cdf_matrix = "GAMMA FITS " + string_of_conditions + ".csv"
    file_headers_gamma_cdf_matrix = ["decrement_value number", "increment_max", "replenishment_scaling", "shape", "scale", "ceiling_rate","gamma_r_squared"]
    file_header_string_gamma_cdr_matrix = ','.join(file_headers_gamma_cdf_matrix)
    np.savetxt(file_name_gamma_cdf_matrix, fitted_gamma_cdf_matrix, delimiter = ",", header = file_header_string_gamma_cdr_matrix)
    return None

#Global variables
reserve_max = 100000
initial_reserve_level = 75000
time_to_run = 25000
function_type = "hyperbolic"
bin_size = 500
base_operant_level = 0

#conditions
increment_max_conditions = np.array([.03])*reserve_max
decrement_value_conditions = np.array([.01])*reserve_max
variable_interval_conditions = np.array([1,2,3,5,8,10,18,25,68,112,200])
replenishment_scaling_conditions = np.array([6.345,5.713,5.081,4.451,3.823,3.200,2.587,2.009,1.440,0.914,0.443])

#condition information
conditions_array = [decrement_value_conditions, increment_max_conditions, replenishment_scaling_conditions, variable_interval_conditions]
levels_per_factor = [len(conditions_array[conditions]) for conditions in range(0,len(conditions_array))]
number_of_sessions = np.product(levels_per_factor)
design_matrix = cartesian_product(conditions_array, number_of_sessions, levels_per_factor, entries_per_session = False)

#file writing
write_raw = False

#initialise structures
record_matrices = recording_matrices(conditions_array, levels_per_factor, number_of_sessions, time_to_run)
statistics_matrix = np.append(design_matrix,np.zeros(shape=(len(design_matrix),6)),axis=1) #Add 5 columns for rft rate(col 3) and response rate(col 4), stdev rft rate, and stdev response rate, residuals from hyperbola, residuals from gamma.
base_fitted_values_matrix = cartesian_product(conditions_array[:-1], number_of_sessions, levels_per_factor, entries_per_session = False) #design matrix sans schedule conditions
fitted_hyperbola_matrix = np.append(base_fitted_values_matrix,np.zeros(shape=(len(base_fitted_values_matrix),3)),axis=1) #Add 3 rows, hyperbola parameters: [k, re, and r^2]
fitted_gamma_cdf_matrix = np.append(base_fitted_values_matrix,np.zeros(shape=(len(base_fitted_values_matrix),4)),axis=1) #Add 4 columns, for shape, scale, ceiling_rate, and r^2

working_directory = "D:\Dropbox\Dropbox\Don\Catania Operant Reserve wd"
run(design_matrix, record_matrices, number_of_sessions)
calculate_statistics(conditions_array, levels_per_factor, number_of_sessions, record_matrices[1], bin_size)
fit_hyperbola(statistics_matrix, fitted_hyperbola_matrix, levels_per_factor)
fit_gamma_cdf(statistics_matrix, fitted_gamma_cdf_matrix, levels_per_factor)
output_to_csv(record_matrices, working_directory, conditions_array,write_raw)