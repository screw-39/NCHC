import sqlite3
import pandas as pd
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

# -----------------------------------------
# 1. 讀取資料庫
# -----------------------------------------
DB_PATH = "SYMMETRY_TEST.db"      # ← 修改成你的 SQLite 檔案路徑
TABLE_TS = "TEST_VOLTAGE"    # TIME-VOLTAGE 資料表
TABLE_PARAM = "TEST_PARAMETER"    # 參數表
ELECTRODE_PARAM = "ELECTRODE_PARAMETER"

def load_test_ids():
    conn = sqlite3.connect(DB_PATH)
    # df = pd.read_sql_query(f"SELECT DISTINCT TEST_ID FROM {TABLE_TS} WHERE TEST_ID >= 733", conn)
    df = pd.read_sql_query(f"SELECT DISTINCT TEST_ID FROM {TABLE_TS} WHERE TEST_ID >= 12 AND TEST_ID <= 732", conn)
    # test1: id = 12, test2: id = 733
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
        html.Label("Select Test (via θ):"),
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
        html.Label("Or input θ (°):"),
        dcc.Input(
            id="theta-input",
            type="number",
            value=0,
            step=0.5
        ),
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
    if test_ids[idx] > 732:
        theta_now = (test_ids[idx] - 732) * 0.5
        theta_min = (test_ids[0] - 732) * 0.5
        theta_max = (test_ids[-1] - 732) * 0.5
    else:
        theta_now = (test_ids[idx] - 12) * 0.5
        theta_min = (test_ids[0] - 12) * 0.5
        theta_max = (test_ids[-1] - 12) * 0.5

    return f"θ = {theta_now}°   (range: {theta_min}°  →  {theta_max}°)"

# 1. Slider → Input
@app.callback(
    Output("theta-input", "value"),
    Input("test-slider", "value")
)
def sync_input_from_slider(slider_index):
    if test_ids[slider_index] > 732:
        theta = (test_ids[slider_index] - 732) * 0.5
    else:
        theta = (test_ids[slider_index] - 12) * 0.5
    return theta


# 2. Input → Slider
@app.callback(
    Output("test-slider", "value"),
    Input("theta-input", "value")
)
def sync_slider_from_input(theta):
    if theta is None:
        return 0

    # 由 θ 反推 test_id
    # target_test_id = int(theta / 0.5 + 732)
    target_test_id = int(theta / 0.5 + 12)

    # 找 test_ids 裡最接近的
    nearest_index = min(range(len(test_ids)), key=lambda i: abs(test_ids[i] - target_test_id))

    return nearest_index

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

    if test_id > 732:
        title = f"Voltage vs Time (θ = {(test_id - 732) * 0.5}°)"
    else:
        title = f"Voltage vs Time (θ = {(test_id - 12) * 0.5}°)"

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
    df = load_electrodes(test_id).iloc[:2]

    x = df["X"].iloc[:2]
    y = df["Y"].iloc[:2]
    z = df["Z"].iloc[:2]

    fig = go.Figure()

    for i in range(2):

        color = "red" if i == 0 else "blue"  # 高頻=紅色, 低頻=藍色
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
        x=[0 for _ in range(6)], y=[0 for _ in range(6)], z=[-3+i for i in range(6)],
        mode="markers",
        marker=dict(size=10),
        name=f"Neuron(soma)"
    ))
    fig.add_trace(go.Scatter3d(
        x=[0 for _ in range(10)], y=[0 for _ in range(10)], z=[3+i for i in range(10)],
        mode="markers",
        marker=dict(size=3),
        name=f"Neuron(apic)"
    ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(title="X", range=[15,-15]),
            yaxis=dict(title="Y", range=[15,-15]),
            zaxis=dict(title="Z", range=[15,-15]),
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