import sqlite3
import pandas as pd
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

# -----------------------------------------
# 1. 讀取資料庫
# -----------------------------------------
DB_PATH = "./DB/SYMMETRY.db"      # ← 修改成你的 SQLite 檔案路徑
TABLE_TS = "TEST_VOLTAGE"    # TIME-VOLTAGE 資料表
TABLE_PARAM = "TEST_PARAMETER"    # 參數表
ELECTRODE_PARAM = "ELECTRODE_PARAMETER"

def load_test_ids():
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(f"SELECT DISTINCT TEST_ID FROM {TABLE_TS}", conn)
    conn.close()
    return sorted(df["TEST_ID"].tolist())

def load_timeseries(test_id):
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(
        f"SELECT TIME, VOLTAGE FROM {TABLE_TS} WHERE TEST_ID = {test_id} ORDER BY TIME ASC",
        conn
    )
    conn.close()
    return df

def load_parameters(test_id):
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(
        f"SELECT * FROM {TABLE_PARAM} WHERE TEST_ID = {test_id}",
        conn
    )
    conn.close()
    return df

def load_electrodes(test_id):
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(
        f"SELECT * FROM {ELECTRODE_PARAM} WHERE TEST_ID = {test_id}",
        conn
    )
    conn.close()
    return df

# -----------------------------------------
# 2. Dash App
# -----------------------------------------
external_stylesheets = [dbc.themes.CERULEAN]
app = Dash(__name__, external_stylesheets=external_stylesheets)
test_ids = load_test_ids()

app.layout = dbc.Container([
    dbc.Row([
        html.Div('Test Voltage-Time Dashboard', className="text-primary text-center fs-3")
    ]),

    dbc.Row([
        html.Label("Select Test (via id):"),
        html.Div(id="slider-label", style={"fontSize": "20px", "textAlign": "center"}),
        
        dcc.Slider(
            id="test-slider",
            min=0,
            max=len(test_ids) - 1,
            step=1,
            value=0,  # 初始位置
            marks=None,
        ),
        html.Br(),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Graph(figure={}, id="voltage-plot"),
        ], width=5),
        dbc.Col([
            dcc.Graph(figure={}, id="location-plot"),
        ], width=7),
    ]),

], fluid=True)

# -----------------------------------------
# 3. Callbacks
# -----------------------------------------
@app.callback(
    Output("slider-label", "children"),
    Input("test-slider", "value")
)
def update_slider_label(idx):
    id_now = test_ids[idx]
    id_min = test_ids[0]
    id_max = test_ids[-1]

    return f"test id = {id_now} (range: {id_min} → {id_max})"

@app.callback(
    Output("voltage-plot", "figure"),
    Input("test-slider", "value"),
)
def update_graph(selected_index):

    test_id = test_ids[selected_index]
    df = load_timeseries(test_id)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df["TIME"],
        y=df["VOLTAGE"],
        mode="lines",
        name=f"TEST {test_id}"
    ))

    title = f"Voltage vs Time (id = {test_id})"

    fig.update_layout(
        xaxis_title="Time(ms)",
        yaxis_title="Voltage(mV)",
        title=title,
        hovermode="x unified",
        width=600,
        height=500
    )

    return fig

@app.callback(
    Output("location-plot", "figure"),
    Input("test-slider", "value"),
)
def update_location_plot(selected_index):

    test_id = test_ids[selected_index]
    df = load_electrodes(test_id)

    x = df["X"]
    y = df["Y"]
    z = df["Z"]

    fig = go.Figure()

    for i in range(4):

        color = "red" if i%2 == 0 else "blue"  # 高頻=紅色, 低頻=藍色
        label = f"Electrode {i} ({2040-40*i} Hz)"

        fig.add_trace(go.Scatter3d(
            x=[x[i]], y=[y[i]], z=[z[i]],
            mode="markers+text",
            marker=dict(size=6, color=color),
            text=[label],
            textposition="top center",
            name=label
        ))

    fig.add_trace(go.Scatter3d(
        x=[0 for _ in range(31)], y=[0 for _ in range(31)], z=[-15+i for i in range(31)],
        mode="markers",
        marker=dict(size=10),
        name=f"Neuron(soma)"
    ))
    fig.add_trace(go.Scatter3d(
        x=[0 for _ in range(31)], y=[0 for _ in range(31)], z=[15+i for i in range(31)],
        mode="markers",
        marker=dict(size=3),
        name=f"Neuron(apic)"
    ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(title="X", range=[50,-50]),
            yaxis=dict(title="Y", range=[50,-50]),
            zaxis=dict(title="Z", range=[50,-50]),
            aspectmode='cube'
        ),
        scene_camera=dict(
            eye=dict(x=0, y=2.5, z=0)  # ← 改成沿著 y 軸看的視角
        ),
        title="Electrode Locations",
        # width=800,
        # height=800,
        margin=dict(l=0, r=0, b=0, t=60)
    )

    return fig

# -----------------------------------------
# 4. Run server
# -----------------------------------------
if __name__ == "__main__":
    app.run(debug=True)