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






References:



