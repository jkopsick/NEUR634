from neuron import h,gui
soma = h.Section(name='soma')
soma.insert('pas')
>>> asyn = h.AlphaSynapse(soma(0.5))
>>> asyn.e
0.0
>>> asyn.gmax
0.0
>>> asyn.gmax = 0.1
# Set onset to 1 so we can see the change from baseline
>>> asyn.onset
0.0
>>> asyn.onset = 1
>>> asyn.tau

>>> v_vec=h.Vector()
>>> t_vec=h.Vector()
>>> v_vec.record(soma(0.5)._ref_v)
1.0
>>> t_vec.record(h._ref_t)
1.0
>>> h.top=40
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
LookupError: 'top' is not a defined hoc variable name.
>>> h.stop=40
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: not assignable
>>> h.tstop=40
>>> h.run()
0.0
>>> asyn.gmax=1
>>> asyn.onset = 20
>>> import pylab as plt
>>> plt.ion()
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018e6da650>]
>>> asyn.onset = 1000
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018ee8e290>]
>>> esyn = h.ExpSyn(soma(0.5))
>>> esyn.e
0.0
>>> esyn.tau
0.1
>>> stim=h.NetStim()
>>> stim.number=3
>>> stim.start=0
>>> stim.interval=5
>>> stim.start=9
>>> ncstim = h.NetCon(stim,esyn)
>>> ncstim.delay=1
>>> ncstim.weight[0]=0.04
>>> ncstim.threshold
10.0
>>> h.tstop
40.0
>>> h.run()
0.0
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018e586f50>]
>>> ncstim.weight[0]=0.8
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018e586c10>]
>>> ncstim.weight[0]=8
>>> h.run()
0.0
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018bd1c950>]
>>> h.psection()
soma { nseg=1  L=100  Ra=35.4
	/*location 0 attached to cell 0*/
	/* First segment only */
	insert morphology { diam=500}
	insert capacitance { cm=1}
	insert pas { g_pas=0.001 e_pas=-70}
	insert AlphaSynapse { onset=1000 tau=0.1 gmax=1 e=0}
	insert ExpSyn { tau=0.1 e=0}
}
1.0
>>> ncstim.weight[0]=80
>>> h.run()
0.0
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018ee77a10>]
>>> ncstim.weight[0]=800
>>> h.run()
0.0
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018bd3d790>]
>>> ncstim.weight[0]=800000
>>> h.run()
0.0
>>> plt.plot(t_vec,v_vec)
[<matplotlib.lines.Line2D object at 0x7f018e6daa10>]


