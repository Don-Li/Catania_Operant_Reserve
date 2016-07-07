Catania's Operant Reserve.

Catania's Operant Reserve is an algorithm for free-operant behaviour on reinforcement schedules. It appears in Catania (2005).

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

v.02    7/7/16
    Added functions:
        experiment_design_values()
            Gives information about the design, e.g. number of decrement_value conditions.
        calculate_statistics()
            Calculates rates for responses and reinforcers
        herrnsteins_hyperbola()
            Calculates value for herrnsteins hyperbola
        fit_hyperbola()
            Fits herrnsteins hyperbola to results from parametric experiment
        output_to_csv()
            Outputs results including raw records, statistics, and fitted values to csv
    Changed functions:
        recording_matrices()
            Now supports parametric experiment
        schedule_block()
            Now supports parametric experiment
        run()
            Now supports parametric experiment

v.03    7/7/16
    Changed functions:
        output_to_csv()
            hyerbola csv file now correctly outputs hyperbola parameters


References:



