import multiprocessing
from extract import *
from transform import *
from visualize import *
import time

'''
-------------transfrom-------------
usage(log, unit=300):   log(dataframe) -> {'cpu_use_rate'(dataframe), 'cpu_occupy'(dataframe), 'cpu_occupy_backfill'(dataframe)}
wait_time(log):         log(dataframe) -> log[NCPUS, wait_time(second)]
work_time(log):         log(dataframe) -> log[NCPUS, work_time(second)]
cancel_time(log):       log(dataframe) -> log[NCPUS, cancel_time(second)]
submit_partition(log):  log(dataframe) -> {map[Partition * submit_time(second in weekday)], x_sub, y_sub}
ncpu_job_count(log):    log(dataframe) -> log[#cpu, job_count(cumulative)]

-------------visualize-------------
plot_usage_heatmap(df, df2, title)
plot_time_scatter(data, title)
plot_submit_heatmap(log, title)
plot_cumulative(log, title)

-------------work flow-------------
data -> usage -> plot_usage_heatmap
data -> wait_time -> plot_time_scatter
data -> work_time -> plot_time_scatter
data -> cancel_time -> plot_time_scatter
data -> submit_partition -> plot_submit_heatmap
data -> ncpu_job_count -> plot_cumulative

'''
def work_flow_usage(log):
    usage_data = usage(log)
    plot_usage_heatmap(usage_data['cpu_occupy'].T, usage_data['cpu_occupy_backfill'].T, 'example_heatmap')

def work_flow_wait_time(log):
    wait_time_data = wait_time(log)
    plot_time_scatter(wait_time_data, 'example_wait_time')

def work_flow_work_time(log):
    work_time_data = work_time(log)
    plot_time_scatter(work_time_data, 'example_work_time')

def work_flow_cancel_time(log):
    cancel_time_data = cancel_time(log)
    plot_time_scatter(cancel_time_data, 'example_canceltime')

def work_flow_submit_partition(log):
    submit_partition_data = submit_partition(log)
    plot_submit_heatmap(submit_partition_data, 'example_submit')

def work_flow_ncpu_job_count(log):
    ncpu_job_count_data = ncpu_job_count(log)
    plot_cumulative(ncpu_job_count_data, 'example_ncpu_count')

def run_functions_in_parallel(log):
    functions = [
        work_flow_usage,
        work_flow_wait_time,
        work_flow_work_time,
        work_flow_cancel_time,
        work_flow_submit_partition,
        work_flow_ncpu_job_count
        ]
    process = []
    for function in functions:
        p = multiprocessing.Process(target=function, args=(log,))
        process.append(p)
        p.start()
    for p in process:
        p.join()

def run(data, multiprocessing=0):
    log = extract(data)
    if multiprocessing:
        run_functions_in_parallel(log)

    else:
        work_flow_usage(log)
        work_flow_wait_time(log)
        work_flow_work_time(log)
        work_flow_cancel_time(log)
        work_flow_submit_partition(log)
        work_flow_ncpu_job_count(log)

if __name__ == '__main__':
    starttime = time.time()
    run('../log/example.log')
    print('USED TIME: {:0.3f} seconds'.format(time.time() - starttime))
    
    
    
    
    
    

     
    
    
    
    
    
