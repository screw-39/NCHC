import numpy as np
import neuron
import LFPy
import os
import plotly.graph_objects as go
import pandas as pd

neuron.h.load_file("stdrun.hoc")

def instantiate_cell(cellParameters):
    cell = LFPy.Cell(**cellParameters, delete_sections=True)
    cell.set_pos(x=0, y=0, z=0)

    # insert hh mechanism in everywhere, reduced density elsewhere
    for sec in cell.allseclist:
        sec.insert('hh')
        if not 'soma' in sec.name():
            # reduce density of Na- and K-channels to 5% in dendrites
            sec.gnabar_hh = 0.006
            sec.gkbar_hh = 0.0018
            
    return cell

#連續正弦波(從t=200開始)
def generate_sin_wave_pulses(
    width: float,
    t_start: float,
    t_stop: float,
    dt: float,
    stim_elec_params: dict,
    num_electrodes: int
):
    # 創建時間軸
    t_ext = np.arange(0, t_stop, dt)

    # 初始化所有電極的電流矩陣 (每個電極預設為 0)
    current = np.zeros((num_electrodes, len(t_ext)))

    # 針對每個指定的刺激電極
    for el_id, params in stim_elec_params.items():
        amp = params["amp"]
        freq = params["freq"]
        phase = params["phase"]

        pulse_start = t_start   # 計算每個脈衝的開始時間
        pulse_end = pulse_start + width  # 計算每個脈衝的結束時間
            
        # 找出 `t_ext` 中符合這個脈衝時間範圍的索引
        pulse_indices = np.where((t_ext >= pulse_start) & (t_ext < pulse_end))[0]

        # 產生正弦波，不同電極擁有不同的相位、頻率與振幅
        theta = 2 * np.pi * freq * (t_ext[pulse_indices]) + phase
        # theta = np.pi * freq * (t_ext[pulse_indices]) + phase
        sin_wave = amp * ((np.sin(theta)))

        current[el_id, pulse_indices] = sin_wave

    return current, t_ext

def main(j, k, l):
    # ---------- Simulation parameters ----------
    cellParameters = {
        'morphology' : './model/ball_and_stick.hoc',
        'tstart' : 0, # ignore startup transients
        'tstop' : 20,
        'dt' : 2**-6,
        'v_init' : -60, 
        'passive' : False,
    }

    # class RecExtElectrode parameters:
    # 1. 定義這三個角度 (由您的模擬環境提供)
    # theta: 俯仰角 (Pitch) - 決定 P0 的高低
    # ro     : 方位角 (Yaw)   - 決定 P0 的水平方向
    # roll   : 側傾角 (Roll)  - 決定正方形是否要「歪著頭」旋轉 <--- 新增這個
    # 如果沒有 roll 變數，可以設為 ro 或 theta 來製造連動的翻滾效果，或設常數 np.pi/4
    theta = j * np.pi / 18
    ro = k * np.pi / 18
    roll = l * np.pi / 18
    R = 10

    # 2. 定義旋轉矩陣
    # 這裡是 Rz(ro) * Ry(theta) * Rx(roll) 的組合
    # 這種寫法保證了不管怎麼轉，三個軸永遠垂直

    # 預先計算三角函數
    c_r, s_r = np.cos(ro), np.sin(ro)
    c_t, s_t = np.cos(theta), np.sin(theta)
    c_l, s_l = np.cos(roll), np.sin(roll)

    # 建構旋轉矩陣 (3x3 Matrix)
    # 這個矩陣代表：先把物體側傾(roll)，再抬頭(theta)，再水平轉(ro)
    R_matrix = np.array([
    [c_r*c_t,  c_r*s_t*s_l - s_r*c_l,  c_r*s_t*c_l + s_r*s_l],
    [s_r*c_t,  s_r*s_t*s_l + c_r*c_l,  s_r*s_t*c_l - c_r*s_l],
    [-s_t,     c_t*s_l,                c_t*c_l]
    ])

    # 3. 定義初始形狀 (原本躺好的正八面體)
    # 這裡定義三根軸的「原始長度向量」
    # Axis 1 (原本指向上): (10, 0, 0)
    # Axis 2 (原本指向前): (0, 10, 0)
    # Axis 3 (原本指向右): (0, 0, 10)
    p_init = np.array([
    [R, 0, 0],   # P0 原始點
    [-R, 0, 0],  # P1 原始點
    [0, R, 0],   # P2 原始點
    [0, -R, 0],  # P3 原始點
    [0, 0, R],   # P4 原始點
    [0, 0, -R],  # P5 原始點
    ]).T # 轉置以便矩陣相乘

    # 4. 進行旋轉 (矩陣乘法)
    # 這行代碼會同時算出所有 6 個點的新座標
    p_rotated = np.dot(R_matrix, p_init)

    # 5. 取出結果
    x = p_rotated[0, :]
    y = p_rotated[1, :]
    z = p_rotated[2, :]

    # 如果 theta_1 和 ro 是陣列，上面的 ax3_z = 0 需要改成 np.zeros_like(theta_1)
    electrodeParameters = dict(
        x=x,
        y=y,
        z=z,
        N=np.array([[0., 0., 1.] for _ in range(6)]),
        r=20.,  # 5um radius
        n=50,  # nb of discrete point used to compute the potential
        sigma=1,  # conductivity S/m
        method="linesource"
        )

    # create cell:
    cell = instantiate_cell(cellParameters)

    # Set stimulation parameters for one electrode
    width1 = 1    # 脈衝寬度pluse width (ms)
    t_start = 2   # (ms)
    t_stop = cell.tstop
    dt = cell.dt

    amp1 = 0.5519*1e5 # 振幅 (nA)
    amp2 = 0.276*1e5
    frequency = 1000
    delta = 20
    stim_elec_params = {
        0:  {"amp": amp2, "freq": frequency + delta, "phase": np.pi }, #+x
        1:  {"amp": amp2, "freq": frequency, "phase": np.pi },         #-x
        2:  {"amp": amp2, "freq": frequency + delta, "phase": np.pi }, 
        3:  {"amp": amp2, "freq": frequency, "phase": np.pi }, 
    }

    # ---- 對每個 cell 套用外加刺激（每次皆使用「新的」probe，避免快取形狀衝突）----
    # 呼叫函數
    electrode = LFPy.RecExtElectrode(cell=cell, **electrodeParameters)
    current, t_ext = generate_sin_wave_pulses(
                width=width1,t_start=t_start, t_stop=t_stop, dt=dt,
                stim_elec_params=stim_elec_params, num_electrodes=6
            )     
    currents = np.array(current)
    electrode.probe.set_currents(currents)
    v_ext = cell.enable_extracellular_stimulation(electrode, t_ext, n=5)

    # run simulation:
    SPIKES = cell.simulate(
        probes=[electrode],
        rec_vmem=True
    )

    t = []
    voltage = []
    for i, v in enumerate(cell.somav):
        t.append(i)
        voltage.append(v)
    df = pd.DataFrame({'TIME':t, 'VOLTAGE':voltage})
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["TIME"],
        y=df["VOLTAGE"],
        mode="lines",
    ))
    fig.show()


    fig = go.Figure()
    for i in range(6):
        color = "red" if (i == 0) or (i == 2) else "blue"  # 高頻=紅色, 低頻=藍色
        label = f"Electrode {i}"
        fig.add_trace(go.Scatter3d(
            x=[x[i]], y=[y[i]], z=[z[i]],
            mode="markers+text",
            marker=dict(size=6, color=color),
            text=[label],
            textposition="top center",
            name=label
        ))
    fig.add_trace(go.Scatter3d(
    x=[0, 0], y=[0, 0], z=[15, -15],
    line=dict(
        color='yellow',
        width=15
        )
    ))
    fig.add_trace(go.Scatter3d(
    x=[0, 0], y=[0, 0], z=[15, 45],
    line=dict(
        color='orange',
        width=3
        )
    ))
    fig.update_layout(
        scene=dict(
            xaxis=dict(title="X", range=[50,-50]),
            yaxis=dict(title="Y", range=[50,-50]),
            zaxis=dict(title="Z", range=[50,-50]),
            aspectmode='cube'
            ),
        )
    fig.show()


if __name__ == "__main__":
    main(2, 0, -2)