# Code for "A Macroeconomic Model with FinanciallyConstrained Intermediaries and Producers"
## Vadim Elenev, Tim Landvoigt, Stijn Van Nieuwerburgh<sup>1</sup>


Fast Steps to Reproduce Benchmark Model
=======================================

-   In Matlab, execute the script `main_create_env`.

    -   If you have the Matlab Symbolic Math Toolbox and use version
        2018b or earlier, leave the flag `useJacobian` in line 43 on.
        Otherwise set to false to use numerical differentiation.

    -   Note: generating the analytic Jacobian for the benchmark model
        takes approximately 5 minutes with version 2018b, and can take
        longer for the other experiments.

    -   `main_create_env` will create a file `env_bench_ini0` that
        contains the experiment definition for the benchmark economy on
        the initial coarse grid.

-   Execute script `main_run_exper` to run the benchmark on the coarse
    grid for up to 100 iterations.

    -   You can set the number of parallel workers to be started in
        line 12.

    -   Set to zero if you want to run it with a single process.

    -   On a computer with sixteen cores (and 16 parallel workers) this
        should take about 30 minutes.

    -   `main_run_exper` creates a results file named
        `res_[current_date_time]` that contains the converged policy
        functions.

    -   Rename this file to `res_20200910_bench_i100.mat`.

-   Execute `main_create_env` again using results file from previous
    step as input.

    -   Configure `main_create_env` as follows (leave other variables as
        is):

        -   `xpername = 'econ';`

        -   `guess_mode = 'guess';`

        -   `guess_path = 'res_20200910_bench_i100';`

    -   `main_create_env` will create a file `env_bench` that contains
        the experiment definition for the benchmark economy on the fine
        grid, using the resulting policy functions from the coarse grid
        iterations in `res_20200910_bench_i100.mat` as initial guess.

-   Execute script `main_run_exper` to run the benchmark on the fine
    grid for up to 30 iterations.

    -   Configure `main_run_exper` as follows (leave other variables as
        is):

        -   `exper_path = 'env_bench.mat';`

        -   `maxit = 30;`

        -   `price_zns = true;`

    -   Set to zero if you want to run it with a single process.

    -   On a computer with sixteen cores (and 16 parallel workers) this
        should take about 45 minutes.

    -   `main_run_exper` creates a results file named
        "`res_[current_date_time]`" that contains the converged policy
        functions.

    -   Rename this file to `res_20200910_bench_s130.mat`.

-   Simulate the model using `sim_stationary` and `sim_trans_cluster`.

    -   `sim_stationary` simulates the model policies contained in
        `res_20200910_bench_s130.mat` for 10,000 periods and writes out
        the resulting time-series and several statistics. The main
        output is a file named `sim_res_20200910_bench_s130.mat`.

    -   `sim_trans_cluster` reads both `res_20200910_bench_s130.mat` and
        `sim_res_20200910_bench_s130.mat`, and simulates generalized
        IRFs.

    -   To plot IRFs, run `plot_trans`.
	
**For More Details See readme_ECTA.pdf and readme_replication.txt**
	

<sup>1</sup>: Elenev: Johns Hopkins University, Carey Business School; email:
    <velenev@jhu.edu>. Landvoigt: University of Pennsylvania Wharton
    School, NBER, and CEPR; email: <timland@wharton.upenn.edu>. Van
    Nieuwerburgh: Columbia University Graduate School of Business, NBER,
    and CEPR; email: <svnieuwe@gsb.columbia.edu>.
