import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

'''
usage(log, unit=300):   log(dataframe) -> {'cpu_use_rate'(dataframe), 'cpu_occupy'(dataframe), 'cpu_occupy_backfill'(dataframe)}
wait_time(log):         log(dataframe) -> log[NCPUS, wait_time(second)]
work_time(log):         log(dataframe) -> log[NCPUS, work_time(second)]
cancel_time(log):       log(dataframe) -> log[NCPUS, cancel_time(second)]
submit_partition(log):  log(dataframe) -> {map[Partition * submit_time(second in weekday)], x_sub, y_sub}
ncpu_job_count(log):    log(dataframe) -> log[#cpu, job_count(cumulative)]
'''


def plot_usage_heatmap(df, title, df2=None):
    # Draw heatmap
    fig = go.Figure(data=go.Heatmap(
        z=df.values,
        x=df.columns,
        y=df.index,
        colorscale=[[0, 'rgb(255,255,255)'], [0.0001, 'rgb(200,200,255)'], [1, 'rgb(0,0,255)']],
        showscale=False,
        colorbar=dict(thickness=20, ticklen=4),
        zmin=0,  # 最小值為0
        zmax=df.values.max()
    ))

    if df2 != None:
        fig.add_trace(
            go.Heatmap(
            z=df2.values,
            x=df2.columns,
            y=df2.index,
            colorscale=[[0, 'rgb(255,255,255)'], [0.0001, 'rgb(255,200,200)'], [1, 'rgb(255,0,0)']],
            showscale=False,
            opacity=0.5,
            colorbar=dict(thickness=20, ticklen=4),
            zmin=0,  # 最小值為0
            zmax=df2.values.max()
        )
        )

    fig.update_layout(
        title={'text':title, 'font':{'size': 70}},
        xaxis_nticks=36,
        plot_bgcolor='White',  # 將背景設置為白色
        #width=500,  # 圖的寬度
        height=1800,
        yaxis={'tickfont':{'size':60}},  # 調整y軸標籤字體大小
        xaxis={'tickfont':{'size':60}},
        #xaxis={'tickfont':{'size':60}, 'range':['2024-07-01T00:00:00','2024-07-17T00:00:00']},  # 調整x軸標籤字體大小
        showlegend=False
    )

    fig.update_traces(hoverongaps=False)  # 不顯示空值的tooltip
    fig.update_traces(zmid=0, colorbar=dict(
        tickvals=[0, df.values.max()],
        ticktext=['0', str(df.values.max())]
    ))

    fig.add_shape(
        type='line', line=dict(dash='solid'),
        name = 'ct112,448 vs ct4k,8k',
        #%Y-%m-%dT%H:%M:%S
        x0 = df.columns[0],
        x1 = df.columns[-1],
        y0 = "icpnp101",
        y1 = "icpnp101"
    )

    fig.add_shape(
        type='line', line=dict(dash='solid'),
        name = 'ct4k,8k vs ct1k,2k',
        #%Y-%m-%dT%H:%M:%S
        x0 = df.columns[0],
        x1 = df.columns[-1],
        y0 = "icpnq256",
        y1 = "icpnq256"
    )

    fig.show()
    pio.write_image(fig, f'{title}.png', width=24*200, height=16*200, scale=2)

def plot_time_scatter(data, y_title, title):
    # data[#CPU, time]
    fig = go.Figure(go.Scatter(
        x=data.iloc[:, 0],
        y=data.iloc[:, 1],
        name=y_title,
        mode='markers',
        marker={'opacity':0.4}
        ))
    
    average_time = data.groupby('#CPU').describe().iloc[:,1]

    fig.update_traces(marker_color='rgb(0, 0, 0)', marker_line_color='rgb(0, 0, 0)',
                        marker_line_width=1.5, opacity=0.3)

    fig.add_trace(go.Scatter(x=average_time[average_time > 1].index, y=average_time[average_time > 1], name='average'))

    fig.update_layout(
        plot_bgcolor='White',
        yaxis={
            'type':'log'
        },
        xaxis={
            'title':'#CPU',
            'type':'log',
            'range':[-0.1,4.8]
        },
        title={'text':title, 'font':{'size':30}},
        xaxis_title={'text':'#CPU', 'font':{'size': 30}},
        yaxis_title={'text':'Seconds', 'font':{'size':30}},
        height=600,
        width=1700
    )

    fig.add_shape(
            type='line', line=dict(dash='dash'),
            name = 'min',
            label = {'text':'min'},
            x0 = 0,
            x1 = 60000,
            y0 = 60,
            y1 = 60
        )

    fig.add_shape(
            type='line', line=dict(dash='dash'),
            name = 'hour',
            label = {'text':'hour'},
            x0 = 0,
            x1 = 60000,
            y0 = 3600,
            y1 = 3600
        )

    fig.add_shape(
            type='line', line=dict(dash='dash'),
            name = 'day',
            label = {'text':'day'},
            x0 = 0,
            x1 = 60000,
            y0 = 3600*24,
            y1 = 3600*24
        )

    fig.add_shape(
            type='line', line=dict(dash='dash'),
            name = 'week',
            label = {'text':'week'},
            x0 = 0,
            x1 = 60000,
            y0 = 3600*24*7,
            y1 = 3600*24*7
        )

    fig.show()
    #pio.write_image(fig, f'{title}.png', scale=2)

def plot_submit_heatmap(log):
    fig = go.Figure(data=go.Heatmap(
        z=test.values,
        x=test.columns,
        y=test.index,
        colorscale=[[0, 'rgb(255,255,255)'], [0.0001, 'rgb(125,125,255)'], [1, 'rgb(0,0,255)']],
        showscale=False,
        colorbar=dict(thickness=20, ticklen=4),
        zmin=0,  # 最小值為0
        zmax=test.values.max(),
        xaxis = 'x',
        yaxis = 'y'
    ))

    fig.add_trace(go.Bar(
            x = df2['y'],
            y = df2.Submit_time,
            xaxis = 'x2',
            yaxis = 'y',
            marker = {
                'line':{
                'color':'rgb(0, 0, 0)',
                'width':1.5
            }
            },
            orientation='h'
        ))

    fig.add_trace(go.Bar(
            x = df.Partition,
            y = df.y,
            xaxis = 'x',
            yaxis = 'y2',
            marker = {
                'color':'rgb(0, 0, 125)'
                }
        ))

    week_ticks = [0, 86400, 172800, 259200, 345600, 432000, 518400]
    week_labels = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']

    fig.update_layout(
        #title={'text':'test', 'font':{'size': 70}},
        #xaxis_nticks=36,
        plot_bgcolor='White',  # 將背景設置為白色
        width=1200,  # 圖的寬度
        height=1800,
        #yaxis={'tickfont':{'size':60}},  # 調整y軸標籤字體大小
        #xaxis={'tickfont':{'size':60}},
        #xaxis={'tickfont':{'size':60}, 'range':['2024-07-01T00:00:00','2024-07-17T00:00:00']},  # 調整x軸標籤字體大小
        showlegend=False,
        xaxis = dict(
            zeroline = False,
            domain = [0,0.83],
            showgrid = False
        ),
        yaxis = dict(
            zeroline = False,
            domain = [0,0.83],
            showgrid = False
        ),
        xaxis2 = dict(
            zeroline = False,
            domain = [0.85,1],
            showgrid = False
        ),
        yaxis2 = dict(
            type = "log",
            zeroline = False,
            domain = [0.85,1],
            showgrid = False,
            # tickvals=week_ticks,  # 自定義 tick 位置
            # ticktext=week_labels 
        ),
        shapes=[
            # 添加一個矩形作為外框
            dict(
                type="rect",
                x0=0, y0=0, x1=0.83, y1=0.83,
                xref="paper", yref="paper",
                line=dict(color="black", width=1)
            ),
            
            dict(
                type="rect",
                x0=0, y0=0.85, x1=0.83, y1=1,
                xref="paper", yref="paper",
                line=dict(color="black", width=1)
            ),

            dict(
                type="rect",
                x0=0.85, y0=0, x1=1, y1=0.83,
                xref="paper", yref="paper",
                line=dict(color="black", width=1)
            ),

            dict(
                type="rect",
                x0=0.85, y0=0, x1=1, y1=0.83,
                xref="paper", yref="paper",
                line=dict(color="black", width=1)
            ),
        ]
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=24*60*60,
        x1=35,  # 結束位置
        y1=24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=2*24*60*60,
        x1=35,  # 結束位置
        y1=2*24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=3*24*60*60,
        x1=35,  # 結束位置
        y1=3*24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=4*24*60*60,
        x1=35,  # 結束位置
        y1=4*24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=5*24*60*60,
        x1=35,  # 結束位置
        y1=5*24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

    fig.add_shape(
        type="line",
        x0=0,  # 開始位置
        y0=6*24*60*60,
        x1=35,  # 結束位置
        y1=6*24*60*60,
        line=dict(
            color="RoyalBlue",
            width=2,
            dash="dashdot",  # 虛線樣式
        ),
        xref="x2",  # 對應主 x 軸
        yref="y1",  # 對應主 y 軸
    )

def plot_cumulative():
    fig = go.Figure()
    fig_area = px.area(x=df.iloc[:,0], y=df.iloc[:,2])
    for trace in fig_area.data:
        fig.add_trace(trace)
    fig.add_trace(go.Bar(x=df.iloc[:,0], y=df.iloc[:,2]))
    #fig.add_trace(go.Area(x=df.iloc[:,0], y=df.iloc[:,2]))
    #fig = px.area(df, x='#CPU', y='cumulative')
    #fig = px.bar(df, x='#CPU', y='cumulative')

    fig.update_layout(
        title={'text':'#CPU v.s. cumulative count', 'font':{'size':30}},
        xaxis_title={'text':"#CPU", 'font':{'size':30}},
        yaxis_title={'text':'cumulative count', 'font':{'size':30}},
        showlegend=False,
        width=1500,
        yaxis=dict(tickfont=dict(size=20)),  # 調整y軸標籤字體大小
        xaxis=dict(tickfont=dict(size=20))
        )
    fig.update_traces(marker_color='rgb(0, 0, 0)', marker_line_color='rgb(0, 0, 0)',
                        marker_line_width=1.5, opacity=0.6)

    fig.add_shape(
            type='line', line=dict(dash='dash'),
            name = '95%',
            x0 = 0,
            x1 = df.iloc[-1,0],
            y0 = i*0.95,
            y1 = i*0.95
        )
        
    fig.show()
    pio.write_image(fig, f'CPU_vs_cumulative_count.png', scale=2)