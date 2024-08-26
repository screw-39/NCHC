import pandas as pd

def extract(log_path, save_parquet=0, parquet_path=None):
    """
    log(txt) -> df_log(dataframe)

    log_path(str):       path of log file
    save_parquet(0/1):   save the parquet or not(defult not)
    parquet_name(str):   where you want to save the parquet
    """

    log = log_path

    with open(log, 'r', encoding="utf-8") as f:
        lines = f.readlines()

    column_name = lines[0].split('%|%')[:-1]
    df_log = {}
    for  i in range(len(column_name)):
        df_log[column_name[i]] = []

    for line in lines[1:]:
        for i in range(len(column_name)):
            df_log[column_name[i]].append(line.split('%|%')[:-1][i])

    df_log = pd.DataFrame(df_log)
    
    if save_parquet == True:
        df_log.to_parquet(parquet_path)

    return df_log

if __name__ == '__main__':
    pass