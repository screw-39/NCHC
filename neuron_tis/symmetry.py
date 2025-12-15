import numpy as np
import neuron
import LFPy
import sqlite3
import os

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

def create_db():
    #connect to database
    conn = sqlite3.connect('./DB/SYMMETRY.db')

    #creat new table
    conn.execute('''CREATE TABLE TEST_PARAMETER(        
            TEST_ID         INT      NOT NULL,
            RO              REAL     NOT NULL,
            THETA           REAL     NOT NULL,
            AMPLITUDE       REAL     NOT NULL,
            FREQUENCY       REAL     NOT NULL,
            DELTA           REAL     NOT NULL,
            PRIMARY KEY (TEST_ID)
            );''')
    #submit change
    conn.commit()

    conn.execute('''CREATE TABLE TEST_VOLTAGE(      
            TEST_ID INT     NOT NULL,
            TIME    REAL    NOT NULL,
            VOLTAGE REAL    NOT NULL,
            PRIMARY KEY (TEST_ID, TIME),
            FOREIGN KEY (TEST_ID) REFERENCES TEST_PARAMETER(TEST_ID)
            );''')
    conn.commit()

    conn.execute('''CREATE TABLE ELECTRODE_PARAMETER          
            (TEST_ID        INT      NOT NULL,
            ELECTRODE_ID    INT      NOT NULL,
            X               REAL     NOT NULL,
            Y               REAL     NOT NULL,
            Z               REAL     NOT NULL,
            PRIMARY KEY (TEST_ID, ELECTRODE_ID),
            FOREIGN KEY (TEST_ID) REFERENCES TEST_PARAMETER(TEST_ID)     
            );''')
    conn.commit()
    conn.close()

def main(j, k):
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
    theta = j * np.pi / 18
    ro = k * np.pi / 18
    print(f'theta: {theta}')
    electrodeParameters = dict(
        x=np.array([10*np.sin(ro)*np.cos(theta), 10*np.sin(ro+np.pi)*np.cos(theta+np.pi), 0, 0]),
        y=np.array([10*np.sin(theta)*np.sin(ro), -10*np.sin(theta+np.pi)*np.sin(ro+np.pi), 10, -10]),
        z=np.array([10*np.sin(theta), -10*np.sin(theta), 0, 0]),
        N=np.array([[0., 0., 1.] for _ in range(4)]),
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

    amp1 = 0.1905*1e5 # 振幅 (nA)
    amp2 = 0
    frequency = 1000
    delta = 20
    stim_elec_params = {
        0:  {"amp": amp1, "freq": frequency + delta, "phase": np.pi }, #+x
        1:  {"amp": amp1, "freq": frequency, "phase": np.pi },         #-x
        2:  {"amp": amp1, "freq": frequency + delta, "phase": np.pi }, 
        3:  {"amp": amp1, "freq": frequency, "phase": np.pi }, 
    }

    # ---- 對每個 cell 套用外加刺激（每次皆使用「新的」probe，避免快取形狀衝突）----
    # 呼叫函數
    electrode = LFPy.RecExtElectrode(cell=cell, **electrodeParameters)
    current, t_ext = generate_sin_wave_pulses(
                width=width1,t_start=t_start, t_stop=t_stop, dt=dt,
                stim_elec_params=stim_elec_params, num_electrodes=4
            )     
    currents = np.array(current)
    electrode.probe.set_currents(currents)
    v_ext = cell.enable_extracellular_stimulation(electrode, t_ext, n=5)

    # run simulation:
    SPIKES = cell.simulate(
        probes=[electrode],
        rec_vmem=True
    )
    if not os.path.isfile('./DB/SYMMETRY.db'):
        create_db()

    conn = sqlite3.connect('./DB/SYMMETRY.db')
    c = conn.cursor()

    try:
        TEST_ID = c.execute("SELECT TEST_ID FROM TEST_PARAMETER;").fetchall()[-1][0]
        TEST_ID = int(TEST_ID) + 1
    except:
        TEST_ID = 0

    print(f"TEST_ID: {TEST_ID}")
    print("edit TEST_PARAMETER")
    c.execute(f'''INSERT INTO TEST_PARAMETER (TEST_ID,RO,THETA,AMPLITUDE,FREQUENCY,DELTA)
      VALUES ({TEST_ID}, {ro}, {theta}, {amp1}, {frequency}, {delta} );''')
    
    print("edit ELECTRODE_PARAMETER")
    for i in range(4):
        c.execute(f'''INSERT INTO ELECTRODE_PARAMETER (TEST_ID,ELECTRODE_ID,X,Y,Z)
        VALUES ({TEST_ID}, {i}, {electrodeParameters['x'][i]}, {electrodeParameters['y'][i]}, {electrodeParameters['z'][i]} );''')

    print("edit TEST_VOLTAGE")
    #PLot voltage of soma to see if neuron has spike.
    t = 0
    for v in cell.somav:
        c.execute(f'''INSERT INTO TEST_VOLTAGE (TEST_ID,TIME,VOLTAGE)
        VALUES ({TEST_ID}, {t}, {v} );''')
        t += dt

    conn.commit()
    conn.close()

if __name__ == "__main__":
    for j in range(37):
        for k in range(37):
            main(j, k)