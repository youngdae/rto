# RTO

We implement a real-time optimization strategy based on warm-start for a moving horizon of multiperiod ACOPFs.

# How to run
```shell
$ julia --project=. src/mpc.jl "data/case{NAME}" "data/case{NAME}/halfhour_30" T H LS RS W OPT cutline cutgen perturb QP sol pfsolve result
```
where
* data/case{NAME}: str, the name of the case file
* data/case{NAME}/halfhour_30: str, the name of the scenario file
* T: int, the length of a time horizon, e.g., T=10 for 10 time periods in a horizon.
* H: int, the number of times the horizon rolls forward, e.g., H=2 means we solve a multiperiod ACOPF once and move forward one time period and solve it again.
* LS: float, load scale in (0.0,1.0], e.g., 1.0 means we do not scale the load profile.
* RS: float, ramp scale in (0.0,1.0], e.g., RS=0.01 means that we allow generator ramp rate to be a 1 percent.
* W: str in {"cold", "shift_copy", "shift_phase1"}. This is for warm-starting of a rolling horizon.
* OPT: int, the option file number IPOPT will be using. OPT=2 will use the default file ipopt.opt.
* cutline: int, 1 if we want to cut a line, 0 otherwise.
* cutgen: int, 1 if we want to turn off a generator, 0 otherwise.
* perturb: int, 1 if we want to perturb the entire load in a time horizon, 0 otherwise.
* QP: int, 1 if we want to approximately solve the multiperiod problem using QP approximation, 0 otherwise.
* SOL: int, 1 if we want to load from a file, 0 otherwise.
* pfsolve: int, 1 if we want to solve a power flow problem for cold-start, 0 otherwise.
* result: str, the name of the result filename

For example,
```shell
$ julia --project=. src/mpc.jl "data/case9" "data/case9/halfhour_30" 10 1 1.0 0.01 "cold" 2 0 1 0 0 0 0 "case9_result"
```
will solve case9 problem with T=10 and one generator turned off.

One could also run the problem in REPL by defining ARGS appropriately.

Alternatively, you could use gensol.sh script in the src directory. After setting appropriate values in it, you could run
```shell
$ src/gensol.sh
```
