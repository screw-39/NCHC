from datetime import datetime
import re
import numpy as np
import pandas as pd

'''
This part will catch the formated dataframe from extract.py and transform the data into muiltiple types.
Muiltiple types of data transform form this part will give to next part(visual.py) to visualize the log.
---------------------------------------------------------------------------------------------------------
usage(log, unit=300):   log(dataframe) -> {'cpu_use_rate'(dataframe), 'cpu_occupy'(dataframe), 'cpu_occupy_backfill'(dataframe)}
wait_time(log):         log(dataframe) -> log[NCPUS, wait_time(second)]
work_time(log):         log(dataframe) -> log[NCPUS, work_time(second)]
cancel_time(log):       log(dataframe) -> log[NCPUS, cancel_time(second)]
submit_partition(log):  log(dataframe) -> {map[Partition * submit_time(second in weekday)], x_sub, y_sub}
ncpu_job_count(log):    log(dataframe) -> log[#cpu, job_count(cumulative)]
'''

def time_translator(time):
    # translate time(yyyy-mm-ddThh:mm:ss) into time(sec)
    return int(datetime.strptime(time, "%Y-%m-%dT%H:%M:%S").timestamp())

def weekday_time_translator(time):
    # translate time(yyyy-mm-ddThh:mm:ss) into weekday time(sec)
    # 解析時間標籤
    dt = datetime.strptime(time, "%Y-%m-%dT%H:%M:%S")
    
    # 取得星期數，星期一為0，日為6
    weekday = dt.weekday()

    # 計算星期一到當下的總秒數
    seconds_since_monday = weekday * 24 * 60 * 60 + dt.hour * 60 * 60 + dt.minute * 60 + dt.second

    return seconds_since_monday

def time_to_seconds(time):
    # translate time(hh:mm:ss) into time(sec)
    t = datetime.strptime(time, "%H:%M:%S")
    return t.hour * 3600 + t.minute * 60 + t.second

def trans_x_time(df, unit, unittime, trans_firsttime):
    # cpu_use_rate.index * unittime + trans_firsttime
    # translate df.index(time block) into time(%Y-%m-%dT%H:%M:%S)
    time = [datetime.fromtimestamp(i * unittime + trans_firsttime).strftime('%Y-%m-%dT%H:%M:%S') for i in range(unit)]
    df.index = time

def init_the_time(firsttime, lasttime, unit):
    # set the time data, return trans_firsttime, trans_lasttime, unittime
    trans_lasttime = time_translator(lasttime)
    trans_firsttime = time_translator(firsttime)
    unittime = (trans_lasttime - trans_firsttime) // unit
    return trans_firsttime, trans_lasttime, unittime

def extend_node_list(nodelist):
    # icpnp[101-103] -> [icpnp101, icpnp102, icpnp103]
    # 使用正規表達式來解析字串
    pattern = r'([a-zA-Z]+)(\[\d+-\d+(?:,\d+-\d+)*\])'
    matches = re.findall(pattern, nodelist)

    # 將符合的結果整理成 list
    result = []
    for match in matches:
        ranges = match[1][1:-1].split(',')
        for r in ranges:
            start, end = map(int, r.split('-'))
            for i in range(start, end + 1):
                result.append(f"{match[0]}{i}")
    return result

def is_valid_datetime_format(date_string):
    try:
        datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S")
        return True
    except ValueError:
        return False
    
def make_data(unit):
    # Creat the zero matrix
    data = np.zeros((unit, 598))
    # Translate into dataFrame
    data = pd.DataFrame(data)
    # Set the node name
    columns = []
    for i in range(101, 157):
        columns.append(f'icpnq{i}')
    for i in range(201, 257):
        columns.append(f'icpnq{i}')
    for i in range(301, 357):
        columns.append(f'icpnq{i}')
    for i in range(401, 457):
        columns.append(f'icpnq{i}')
    for i in range(501, 557):
        columns.append(f'icpnq{i}')
    for i in range(601, 657):
        columns.append(f'icpnq{i}')
    for i in range(701, 757):
        columns.append(f'icpnq{i}')
    for i in range(101, 157):
        columns.append(f'icpnp{i}')
    for i in range(201, 257):
        columns.append(f'icpnp{i}')
    for i in range(301, 349):
        columns.append(f'icpnp{i}')
    for i in range(1, 7):
        columns.append(f'gpn0{i}')
    for i in range(1, 41):
        columns.append(f'ncpn{i}')
    # Renew column name of dataFrame
    data.columns = columns
    return data

def make_submit_map(log):
    '''
    creat a blank dataframe for submit_partition map
    size: [partition * time]
    '''
    Partition = ['development', 'ct112', 'ct448', 'ct1k', 'ct2k', 'ct4k', 'ct8k', 'arm-dev', 'arm144', 'arm576', 'arm1440', 'dc-ENT112024-01', 'visual-dev', 'visual', 'vscode', 'jupyter']
    Submit = []
    
    for time in log.Submit:
        if time not in Submit:
            Submit.append(time)

    # Creat the zero matrix
    heatmap = np.zeros((len(Submit), len(Partition)))
    heatmap = pd.DataFrame(heatmap)
    heatmap.columns = Partition
    heatmap.index = Submit

    return heatmap

def add_value_to_data(data, timestart, timeend, nodelist, value):
    """
    Add value on data[timestart:timeend+1, nodelist]
    
    Parameters:
        data (DataFrame):         DataFrame to catch the value(from make_data())
        timestart (int):          Begin time of the job
        timeend (int):            End time of the job
        nodelist (list):          Node list of the job
        value (float or list):    Value you want to add
    """
    if isinstance(value, list):
        for i in range(len(value)):
            data.loc[timestart+i, nodelist] += value[i]
    else:
        data.loc[timestart:timeend+1, nodelist] += value
    return data

def job_processor(row, trans_firsttime, unittime):
    """
    Process the job imformation into:
        {'jobstart', 'jobend', 'nodelist', 'allocated_cpu', 'backfill_allocated_cpu', 'cpu_useratio'}
        jobstart (int):                *Begin unit* of the job
        jobend (int):                  *End unit* of the job
        nodelist (list[str]):          All nodes the job works
        allocated_cpu (list[float]):   Number of allocated CPU
        cpu_useratio (list[float]):    Utilization of CPU    <-----------------------(wait to fix): fix function of useratio and time_to_seconds('%D-' in the time)
        allocated_cpu(backfill):       Number of allocated CPU(backfill)
    
    Parameters:
        row (DataFrame):               One job imformation from dataframe
        trans_firsttime (int):         Begin time of timeline(sec)
        unittime (int):                Length of one unit time(sec)
    """
    
    try:
        jobstart = (time_translator(row.Start) - trans_firsttime) // unittime
        if ((time_translator(row.Start) - trans_firsttime) % unittime) != 0:
            jobstart += 1
    except:
        jobstart = None
    try:
        jobend = (time_translator(row.End) - trans_firsttime) // unittime
        if ((time_translator(row.End) - trans_firsttime) % unittime) != 0:
            jobend += 1
    except:
        jobend = None

    nodelist = row.NodeList
    if '[' in nodelist: nodelist = extend_node_list(nodelist)
    
    # allocated_cpu part
    try: 
        #row.AllocCPUS / int(row.AllocNodes)
        allocated_value = int(row.AllocCPUS) / int(row.AllocNodes)
        allocated_cpu = [allocated_value for _ in range(jobend - jobstart + 1)]
        
        if (jobend - jobstart) > 0:
            allocated_cpu[0] = allocated_value*(1 - (((time_translator(row.Start) - trans_firsttime) % unittime) / unittime))
            allocated_cpu[-1] = allocated_value*(((time_translator(row.End) - trans_firsttime) % unittime) / unittime)
        else:
            allocated_cpu[0] = allocated_value*((time_translator(row.End) - time_translator(row.Start)) / unittime)
        
    except:
        allocated_cpu = [0]

    # allocated_cpu_backfill part
    if 'Backfill' in row.Flags:
        try: 
            #row.AllocCPUS / int(row.AllocNodes)
            allocated_value_Backfill = int(row.AllocCPUS) / int(row.AllocNodes)
            allocated_cpu_Backfill = [allocated_value_Backfill for _ in range(jobend - jobstart + 1)]
            
            if (jobend - jobstart) > 0:
                allocated_cpu_Backfill[0] = allocated_value_Backfill*(1 - (((time_translator(row.Start) - trans_firsttime) % unittime) / unittime))
                allocated_cpu_Backfill[-1] = allocated_value_Backfill*(((time_translator(row.End) - trans_firsttime) % unittime) / unittime)
            else:
                allocated_cpu_Backfill[0] = allocated_value_Backfill*((time_translator(row.End) - time_translator(row.Start)) / unittime)
            
        except:
            allocated_cpu_Backfill = [0]

    else:
        allocated_cpu_Backfill = [0]
    
    # cpu_useratio part
    # try:
    #     useratio_value = (time_to_seconds(row.TotalCPU) / (time_translator(row.End) - time_translator(row.Start)) / int(row.AllocCPUS))
    #     cpu_useratio = [useratio_value for _ in range(jobend - jobstart + 1)]
        
    #     if (jobend - jobstart) > 0:
    #         cpu_useratio[0] = useratio_value*(1 - (((time_translator(row.Start) - trans_firsttime) % unittime) / unittime))
    #         cpu_useratio[-1] = useratio_value*(((time_translator(row.End) - trans_firsttime) % unittime) / unittime)
    #     else:
    #         cpu_useratio[0] = useratio_value*((time_translator(row.End) - time_translator(row.Start)) / unittime)

    # except:
    #     cpu_useratio = [0]

    result = {'jobstart':jobstart, 'jobend':jobend, 'nodelist':nodelist, 'allocated_cpu':allocated_cpu, 'backfill_allocated_cpu':allocated_cpu_Backfill}
    return result

def usage(log, unit=300):
    '''
    log(dataframe) -> {'cpu_use_rate'(dataframe), 'cpu_occupy'(dataframe), 'cpu_occupy_backfill'(dataframe)}
    unit(int):         How many unit do you want to split
    '''

    #clean the log for processing
    log = log.query('Group != ""').query('Start != "None"').query('End != "None"').query('Start != "Unknown"').query('End != "Unknown"').query('NodeList != "None assigned"')
    log = log.loc[:, ['Start', 'End', 'NodeList', 'AllocCPUS', 'AllocNodes', 'TotalCPU', 'Flags']]

    #NODEMAX = 598
    firsttime = log.End.sort_values(ascending=1).iloc[0]
    
    lasttime = log.End.sort_values(ascending=0).iloc[0]
    trans_firsttime, trans_lasttime, unittime = init_the_time(firsttime, lasttime, unit)
    
    # Prepare dataframe for catching the value from job_processor
    # more picture -> more make_data()
    cpu_occupy = make_data(unit)
    cpu_occupy_backfill = make_data(unit)
    cpu_use_rate = make_data(unit)

    
    for index, row in log.iterrows():
        
        result = job_processor(row, trans_firsttime, unittime)
        # result = [jobstart(int), jobend(int), nodelist(list), allocated_cpu(list), backfill_allocated_cpu(list), cpu_useratio(list)]
        
        try:
            # add_value_to_data(data, timestart, timeend, nodelist, value)
            cpu_use_rate = add_value_to_data(cpu_use_rate, result['jobstart'], result['jobend'], result['nodelist'], result['cpu_useratio'])
        except:
            pass

        try:
            # add_value_to_data(data, timestart, timeend, nodelist, value)
            cpu_occupy = add_value_to_data(cpu_occupy, result['jobstart'], result['jobend'], result['nodelist'], result['allocated_cpu'])
        except:
            pass

        try:
            # add_value_to_data(data, timestart, timeend, nodelist, value)
            cpu_occupy_backfill = add_value_to_data(cpu_occupy_backfill, result['jobstart'], result['jobend'], result['nodelist'], result['backfill_allocated_cpu'])
        except:
            pass

    #transform the x-axis of dfs ( unit -> yyyy-mm-ddThh:mm:ss )
    #trans_x_time(cpu_use_rate.index * unittime + trans_firsttime)
    trans_x_time(cpu_use_rate, unit, unittime, trans_firsttime)
    trans_x_time(cpu_occupy, unit, unittime, trans_firsttime)
    trans_x_time(cpu_occupy_backfill, unit, unittime, trans_firsttime)

    return {'cpu_use_rate':cpu_use_rate, 'cpu_occupy':cpu_occupy, 'cpu_occupy_backfill':cpu_occupy_backfill}


def wait_time(log):
    '''
    log(dataframe) -> log[NCPUS, wait_time(second)]
    '''
    log = log.query('Group != ""').query('Start != "None"').query('Submit != "None"').query('Start != "Unknown"').query('Submit != "Unknown"').query('NodeList != "None assigned"')
    log = log.loc[:, ['NCPUS', 'Submit', 'Start']]
    log.Start = log.Start.apply(time_translator)
    log.Submit = log.Submit.apply(time_translator)
    log.NCPUS = log.NCPUS.apply(int)
    log['wait_time'] = log['Start'] - log['Submit']
    log = log.loc[:, ['NCPUS', 'wait_time']].sort_values(by=['NCPUS'])

    average = log.groupby('NCPUS').describe().iloc[:,1]
    average = pd.DataFrame({'NCPUS':average.index, 'average_time':average})

    return {'log':log, 'average':average}

def work_time(log):
    '''
    log(dataframe) -> log[NCPUS, work_time(second)]
    '''
    log = log.query('Group != ""').query('Start != "None"').query('End != "None"').query('Start != "Unknown"').query('End != "Unknown"').query('NodeList != "None assigned"')
    log = log.loc[:, ['NCPUS', 'Start', 'End']]
    log.Start = log.Start.apply(time_translator)
    log.End = log.End.apply(time_translator)
    log.NCPUS = log.NCPUS.apply(int)
    log['work_time'] = log['End'] - log['Start']
    log = log.loc[:, ['NCPUS', 'work_time']].sort_values(by=['NCPUS'])

    average = log.groupby('NCPUS').describe().iloc[:,1]
    average = pd.DataFrame({'NCPUS':average.index, 'average_time':average})

    return {'log':log, 'average':average}

def cancel_time(log):
    '''
    log(dataframe) -> log[NCPUS, cancel_time(second)]
    '''
    log = log.query('Group != ""').query('Start == "None"').query('Submit != "None"').query('End != "None"').query('Submit != "Unknown"').query('End != "Unknown"')
    log = log.loc[:, ['NCPUS', 'Submit', 'End']]
    log.Submit = log.Submit.apply(time_translator)
    log.End = log.End.apply(time_translator)
    log.NCPUS = log.NCPUS.apply(int)
    log['cancel_time'] = log['End'] - log['Submit']
    log = log.loc[:, ['NCPUS', 'cancel_time']].sort_values(by=['NCPUS'])

    average = log.groupby('NCPUS').describe().iloc[:,1]
    average = pd.DataFrame({'NCPUS':average.index, 'average_time':average})

    return {'log':log, 'average':average}

def submit_partition(log):
    '''
    log(dataframe) -> {map[Partition * submit_time(second in weekday)], x_sub, y_sub}
    '''
    log = log.query('Group != ""').query('Submit != "None"').query('Submit != "Unknown"')
    log = log.loc[:, ['Submit', 'Partition']]
    log.Submit = log.Submit.apply(weekday_time_translator)
    log = log.sort_values('Submit')
    map = make_submit_map(log)
    for ind, row in log.iterrows():
        try:
            map.loc[row.Submit, row.Partition] += 1
        except:
            continue
    y_sub = pd.DataFrame.from_dict({'Partition':[i for i in map.columns], 'y':[int(i) for i in map.sum().values]})
    x_sub = pd.DataFrame.from_dict({'Submit_time':[i for i in map.index], 'y':[r.sum() for i,r in map.iterrows()]})

    return {'map':map, 'x_sub':x_sub, 'y_sub':y_sub}

def ncpu_job_count(log):
    '''
    log(dataframe) -> log[#cpu, job_count(cumulative)]
    '''
    log = log.query('Group != ""').query('Submit != "None"').query('Submit != "Unknown"').query('NCPUS != ""')
    log = log.loc[:, ['NCPUS']]
    log = log.sort_values('NCPUS')
    log['count'] = [1 for _ in range(len(log.NCPUS))]
    
    cpu = [int(i) for i in log.groupby('NCPUS').describe().index]
    cpu.sort()
    cumulative = []
    count = 0
    for c in cpu:
        count += len(log[log.NCPUS == str(c)]['count'])
        cumulative.append(count)
    
    log = pd.DataFrame({'#CPU':cpu, 'count':cumulative})

    return log