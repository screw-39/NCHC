import pandas as pd
import glob

#file_path = glob.glob(r'./hpc*.txt')[0]
file_path = 'C:/Users/2403037/Documents/sidework/logana/log/20240801_hpc_jobs_log.txt'
log_save_name = 'C:/Users/2403037/Documents/sidework/logana/data/20240801_F1_log.parquet'

with open(file_path, 'r', encoding="utf-8") as f:
    lines = f.readlines()

column_name = lines[0].split('%|%')[:-1]
df_log = {}
for  i in range(len(column_name)):
    df_log[column_name[i]] = []

for line in lines[1:]:
    for i in range(len(column_name)):
        df_log[column_name[i]].append(line.split('%|%')[:-1][i])

df_log = pd.DataFrame(df_log)
df_log.to_parquet(log_save_name)