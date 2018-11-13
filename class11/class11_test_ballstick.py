from neuron import h, gui
import ballandstick as bs
import util as u
cells=[]
N=3
for i in range(N):
    cell = bs.BallAndStick()
    cells.append(cell)

src = cells[0]
tgt = cells[1]
syn = h.ExpSyn(tgt.dend(0.5))
nc = h.NetCon(src.soma(0.5)._ref_v, syn, sec=src.soma)
nc.weight[0] = 0.05 # uS
nc.delay = 5 #ms
soma_v_vec, dend_v_vec, t_vec = u.set_recording_vectors(cells[0])
