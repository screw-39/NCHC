from .src import *
from datetime import datetime, timedelta
from airflow import DAG
from airflow.operators.python_operator import PythonOperator
# BranchPythonOperator
from airflow.operators.dummy_operator import DummyOperator

'''
-------------transfrom-------------
usage(log, unit=300):   log(dataframe) -> {'cpu_use_rate'(dataframe), 'cpu_occupy'(dataframe), 'cpu_occupy_backfill'(dataframe)}
wait_time(log):         log(dataframe) -> log[NCPUS, wait_time(second)]
work_time(log):         log(dataframe) -> log[NCPUS, work_time(second)]
cancel_time(log):       log(dataframe) -> log[NCPUS, cancel_time(second)]
submit_partition(log):  log(dataframe) -> {map[Partition * submit_time(second in weekday)], x_sub, y_sub}
ncpu_job_count(log):    log(dataframe) -> log[#cpu, job_count(cumulative)]

-------------visualize-------------
plot_usage_heatmap(*log1, *log2, title)
plot_time_scatter(*log, title)
plot_submit_heatmap(*log, title)
plot_cumulative(*log, title)

-------------work flow-------------
data -> usage -> plot_usage_heatmap
data -> wait_time -> plot_time_scatter
data -> work_time -> plot_time_scatter
data -> cancel_time -> plot_time_scatter
data -> submit_partition -> plot_submit_heatmap
data -> ncpu_job_count -> plot_cumulative
'''

def extract_log(data):
    log = extract(data)
    return log

def work_flow_usage(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    usage_data = usage(log)
    plot_usage_heatmap(usage_data['cpu_occupy'].T, usage_data['cpu_occupy_backfill'].T, 'example_heatmap')

def work_flow_wait_time(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    wait_time_data = wait_time(log)
    plot_time_scatter(wait_time_data, 'example_wait_time')

def work_flow_work_time(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    work_time_data = work_time(log)
    plot_time_scatter(work_time_data, 'example_work_time')

def work_flow_cancel_time(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    cancel_time_data = cancel_time(log)
    plot_time_scatter(cancel_time_data, 'example_canceltime')

def work_flow_submit_partition(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    submit_partition_data = submit_partition(log)
    plot_submit_heatmap(submit_partition_data, 'example_submit')

def work_flow_ncpu_job_count(**context):
    log = context['task_instance'].xcom_pull(task_ids='extract_log')
    ncpu_job_count_data = ncpu_job_count(log)
    plot_cumulative(ncpu_job_count_data, 'example_ncpu_count')

default_args = {
    'owner': 'Chi-Feng, CHENG',
    'start_date': datetime(2024, 10, 1),
    #@monthly
    'schedule_interval': '@once',
    'retries': 2,
    'retry_delay': timedelta(minutes=1)
}

with DAG(dag_id='airflow_present', default_args=default_args) as dag:

    workflow_start = DummyOperator(
        task_id='workflow_start'
    )

    extract_log = PythonOperator(
        task_id='extract_log',
        python_callable=extract,
        op_args=['../log/example.log'],
        provide_context=True
    )

    work_flow_usage = PythonOperator(
        task_id='work_flow_usage',
        python_callable=work_flow_usage,
        provide_context=True
    )

    work_flow_wait_time = PythonOperator(
        task_id='work_flow_wait_time',
        python_callable=work_flow_wait_time,
        provide_context=True
    )

    work_flow_work_time = PythonOperator(
        task_id='work_flow_work_time',
        python_callable=work_flow_work_time,
        provide_context=True
    )

    work_flow_cancel_time = PythonOperator(
        task_id='work_flow_cancel_time',
        python_callable=work_flow_cancel_time,
        provide_context=True
    )

    work_flow_submit_partition = PythonOperator(
        task_id='work_flow_submit_partition',
        python_callable=work_flow_submit_partition,
        provide_context=True
    )

    work_flow_ncpu_job_count = PythonOperator(
        task_id='work_flow_ncpu_job_count',
        python_callable=work_flow_ncpu_job_count,
        provide_context=True
    )

    workflow_start >> extract_log >> work_flow_usage >> work_flow_wait_time >> work_flow_work_time >> work_flow_cancel_time >> work_flow_submit_partition >> work_flow_ncpu_job_count
    
    
    
    
    
    

     
    
    
    
    
    
