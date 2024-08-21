from datetime import datetime
import re
import numpy as np
import pandas as pd

def time_translator(time):
    # translate time(yyyy-mm-ddThh:mm:ss) into time(sec)
    return int(datetime.strptime(time, "%Y-%m-%dT%H:%M:%S").timestamp())

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
        cpu_useratio (list[float]):    Utilization of CPU (wait to fix)
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
        #row.AllocCPUS / int(row.AllocNodes
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
            #row.AllocCPUS / int(row.AllocNodes
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
        result.append(allocated_cpu_Backfill)
    
    # cpu_useratio part
    try:
        useratio_value = (time_to_seconds(row.TotalCPU) / (time_translator(row.End) - time_translator(row.Start)) / int(row.AllocCPUS))
        useratio_cpu = [useratio_value for _ in range(jobend - jobstart + 1)]
        
        if (jobend - jobstart) > 0:
            useratio_cpu[0] = useratio_value*(1 - (((time_translator(row.Start) - trans_firsttime) % unittime) / unittime))
            useratio_cpu[-1] = useratio_value*(((time_translator(row.End) - trans_firsttime) % unittime) / unittime)
        else:
            useratio_cpu[0] = useratio_value*((time_translator(row.End) - time_translator(row.Start)) / unittime)

    except:
        cpu_useratio = [0]

    result = {'jobstart':jobstart, 'jobend':jobend, 'nodelist':nodelist, 'allocated_cpu':allocated_cpu, 'backfill_allocated_cpu':allocated_cpu_Backfill, 'cpu_useratio':cpu_useratio}
    return result

def usage_heatmap(log, unit=300):
    '''
    log(dataframe) -> usage_heatmap(dataframe)
    unit(int):         How many unit do you want to split
    '''

    #clean the log for processing
    log = log.query('Start != "None"').query('End != "None"').query('Start != "Unknown"').query('End != "Unknown"').query('NodeList != "None assigned"')
    
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
        if row.Group != '': #one job count once
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



