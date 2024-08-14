from datetime import datetime
import pandas as pd

def time_translator(time):
    # translate time(yyyy-mm-ddThh:mm:ss) into time(sec)
    return int(datetime.strptime(time, "%Y-%m-%dT%H:%M:%S").timestamp())

def is_valid_datetime_format(date_string):
    try:
        datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S")
        return True
    except ValueError:
        return False

def job_processor(i, row, trans_firsttime):
    """
    Process the job imformation into:
        jobstart (int):                *Begin unit* of the job
        jobend (int):                  *End unit* of the job
        nodelist (list[str]):          All nodes the job works
        allocated_cpu (list[float]):   Number of allocated CPU
        cpu_useratio (list[float]):    Utilization of CPU
        wait_time (int):               Waiting time of the job
    
    Parameters:
        row (DataFrame):               One job imformation from dataframe
        trans_firsttime (int):         Begin time of timeline(sec)
        unittime (int):                Length of one unit time(sec)
    """

    subtime = time_translator(row.Submit) - trans_firsttime
    #job_id = row.JobID
    work_time = work_time = time_translator(row.End) - time_translator(row.Start)
    user = f'user{i % 5 + 1}'
    if (i % 5 + 1) < 4:
        account = 'account1'
    else:
        account = 'account2'
    partition = row.Partition
    qos = 'normal'
    ntask = row.AllocCPUS
    n_per_node = int(int(row.AllocCPUS)/int(row.AllocNodes))
    if int(row.AllocCPUS)%int(row.AllocNodes) != 0:
        n_per_node += 1

    return f'-dt {subtime} -e submit_batch_job | -J jobid_{i} -sim-walltime {work_time}.00 --uid={user} -t 96:00:00 -n {ntask} --ntasks-per-node={n_per_node} -A {account} -p {partition} -q {qos} pseudo.job'

def main():
    """
    Load dataframe and time setting -> Read row and calculate the value -> Output the job event
    You can set the dataframe(log file) path, time setting and workflow in main function.

    *If you want to analyze more information from the dataframe, it is easy to achieve by modifying the job_processor() and workflow in the main().*

    Parameters:
        log (DataFrame):  Target dataFrame
        firsttime (date): Begin of the data (yyyy-mm-ddThh:mm:ss)
        fn (path name): Path of the events file
        i : Count of the job
    """
    # File path and time setting
    log = pd.read_parquet('C:/Users/2403037/Documents/sidework/slurm/20240801_F1_log.parquet')
    log = log[log.Submit > '2024-07-01T00:00:00']
    firsttime = '2024-07-01T00:00:00'
    fn = 'C:/Users/2403037/Documents/sidework/slurm/jobs0701_31.events'

    #print(firsttime)
    trans_firsttime = time_translator(firsttime)
    i = 1

    with open(fn, 'w') as f:
        for index, row in log.iterrows():
            if row.Group != '': #one job count once
                if is_valid_datetime_format(row.Start) and is_valid_datetime_format(row.End):
                    if row.NodeList != 'None assigned':
                        job = job_processor(i, row, trans_firsttime)
                        f.write(f'{job}\n')
                        i += 1

            # if i == 101:
            #     break
    print(i-1)





if __name__ == '__main__':
    main()