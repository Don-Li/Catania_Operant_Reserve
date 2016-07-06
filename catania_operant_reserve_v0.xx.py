# -*- coding: utf-8 -*-
"""
Current patch notes:
v.01    6/7/16
    Added functions:
        reserve_replenishment()
            For a given delay-to-reinforcement gradient and vector of delays and maximum increment, return the amount that the reserve is increased by for a reinforcer.
        recording_matrices
            Sets up a matrix to record repsonses and reinforcers
        schedule_block
            Run the algorithm on the schedule
    Added variables:
        increment_max
            Gives the maximum contribution of a behaivour to the reserve_replenishment function.
        decrement_value
            How much the reserve is decremented when a behaviour occurs.
        reserve_max
            maximum reserve value
        operat_level
            specifies baseline probability of responding
        initial_reserve_level
            reserve level at beginning of experiment
        time_to_run
            number of iterations of the algorithm
"""
import numpy as np

reserve_max = 100000
increment_max = reserve_max* .03
decrement_value = reserve_max* .01
operant_level = 0
initial_reserve_level = 75000
time_to_run = 25000
random_interval = 8
function_type = "reciprocal"

def reserve_replenishment(delay_vector,increment_max,function_type):
    if function_type == "exponential":
        reserve_replenishment = np.sum(increment_max*np.exp(-delay_vector+1))
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

def recording_matrices(time_to_run):
    record_matrix_0 = np.empty(shape = (time_to_run))
    record_matrix_1 = np.zeros(shape = (time_to_run,2))
    record_matrix_0[:] = np.array(range(0,time_to_run))
    return (record_matrix_0, record_matrix_1)

def schedule_block(reserve_max, increment_max, decrement_value, operant_level, initial_reserve_level, time_to_run, random_interval, function_type, record_matrix_1):
    time = 0
    reinforcement_arranged = False
    last_reinforcer_time_1 = 0
    last_reinforcer_time_2 = 0
    next_reinforcer_time = None
    reserve_now = initial_reserve_level

    while time < time_to_run:
        emitted_behaviour = behaviour_emission(reserve_now, reserve_max)
        if reinforcement_arranged == False:
            next_reinforcer_time = variable_interval_schedule(random_interval, time)
            reinforcement_arranged = True
        if emitted_behaviour == True:
            record_matrix_1[time,1] = 1
            reserve_now -= decrement_value
            if time >= next_reinforcer_time:
                record_matrix_1[time,0] = 1
                last_reinforcer_time_2 = last_reinforcer_time_1
                last_reinforcer_time_1 = time
                delay_vector = np.where(record_matrix_1[last_reinforcer_time_2:last_reinforcer_time_1+1,1][::-1])[0]+1
                reserve_now += reserve_replenishment(delay_vector, increment_max, function_type)
                reinforcement_arranged = False
            if time < next_reinforcer_time:
                record_matrix_1[time,0] = 0
        if emitted_behaviour == False:
            record_matrix_1[time,1] = 0
        time += 1
    return ()


record_matrices = recording_matrices(time_to_run)
schedule_block(reserve_max, increment_max, decrement_value, operant_level, initial_reserve_level, time_to_run, random_interval, function_type, record_matrices[1])
