import numpy as np
import plotly.graph_objects as go

def main():
    theta_1 = np.pi * (90/180)
    theta_2 = 0
    ro = np.pi * (50/180)
    x=np.array([10*np.sin(theta_1)*np.cos(ro), -10*np.sin(theta_1)*np.cos(ro), 10*np.sin(theta_1)*np.cos(ro + (np.pi/2)), -10*np.sin(theta_1)*np.cos(ro + (np.pi/2)), 10*np.sin(theta_1 + (np.pi/2))*np.cos(ro + (np.pi/2)), -10*np.sin(theta_1 + (np.pi/2))*np.cos(ro + (np.pi/2))])
    y=np.array([10*np.sin(theta_1)*np.sin(ro), -10*np.sin(theta_1)*np.sin(ro), 10*np.sin(theta_1)*np.sin(ro + (np.pi/2)), -10*np.sin(theta_1)*np.sin(ro + (np.pi/2)), 10*np.sin(theta_1 + (np.pi/2))*np.sin(ro + (np.pi/2)), -10*np.sin(theta_1 + (np.pi/2))*np.sin(ro + (np.pi/2))])
    z=np.array([10*np.cos(theta_1), -10*np.cos(theta_1), 10*np.cos(theta_1), -10*np.cos(theta_1), 10*np.cos(theta_1 + (np.pi/2)), -10*np.cos(theta_1 + (np.pi/2))])
    
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode="markers+text",
        marker=dict(size=6),
        textposition="top center",
    ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(title="X", range=[15,-15]),
            yaxis=dict(title="Y", range=[15,-15]),
            zaxis=dict(title="Z", range=[15,-15]),
            aspectmode='cube'
        ),
    )

    fig.show()


if __name__ == "__main__":
    main()
