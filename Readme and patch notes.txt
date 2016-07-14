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

v.04    11/7/16
    Changed functions:
        herrnsteins_hyperbola()
            Now takes the base operant level as a specified global variable.

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

v.06    14/7/16
    Changed functions:
        output_to_csv()
            Now makes nicer function names.

References:



