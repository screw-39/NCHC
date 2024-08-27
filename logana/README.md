# logan

This repository is used to analysis and visualize a jobs log on slurm, which is a free and open-source job scheduler, generates by using `sacct`. 

![image](./images/example_heatmap.png)

# Usage

Just need three step to do the analysis.

1. Use following command in the cluster to generate the log file.
```!bash
sacct -pa --noconvert --units=M --delimiter="%|%" --format=AllocCPUS,AllocNodes,CPUTimeRAW,End,ExitCode,Flags,Group,JobID,JobName,NCPUS,NNodes,NodeList,NTasks,Partition,ReqCPUS,ReqMem,ReqNodes,Start,State,Submit,SystemCPU,TimelimitRaw,TotalCPU,UserCPU > example.log
```

2. Put the log file in `./log/`. There is already have an example log for the example. 

3. Edit the log path in `./src/work_flow.py` and execute it.

Wait for the analysis and you could find the result in `./images`.

# How to add self analysis

This repository divide in four parts in `./src`.

`extract.py` : Translates a log file into dataframe.

`transform.py` : Has some basic functions and different analysis definite in functions.

`visual.py` : Functions for generate different images with different analysis.

`work_flow.py` : Set the work flow of a analysis.

If you want to add an analysis by yourself, you can add you anaysis function in `transform.py`. And then use the tool in `visual.py`, or add one by yourself, to generate the image. Finally, add a complete analysis work flow in `work_flow.py`.
