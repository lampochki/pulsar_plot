import struct
import os
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import pylab
from mpl_toolkits.mplot3d import Axes3D


#create('DIFX_56341_040423.s0000.b0000', 258, 'RR', 0, 'filname',)
'''
256*1+1=257
256*1+2=258
256*1+3=259
256*2+1=513
256*2+2=514
256*2+3=515
256*3+1=769
256*3+2=770
256*3+3=771
'''

def create(fn='', bn='', pp='', fi='', outfile='', per=''):
	size=os.path.getsize(fn)
	headnumber=0
	headcount=0
	count=0
	x=0
	with open(fn, 'rb') as sw:
		sync=int(struct.unpack('i', sw.read(4))[0])
		sync_search=0 
		period=0
		while sync_search!=sync:
			period+=1
			sw.seek(period)
			sync_search=int(struct.unpack('i', sw.read(4))[0])
		if per==0:
			while headcount<size/period:
				sw.seek(period*headnumber)
				sync=str(sw.read(4))
				head_version,=struct.unpack('i', sw.read(4))
				baseline_number,=struct.unpack('i', sw.read(4))
				MJD,=struct.unpack('i', sw.read(4))
				seconds,=struct.unpack('d', sw.read(8))
				config_index,=struct.unpack('i', sw.read(4))
				source_index,=struct.unpack('i', sw.read(4))
				freq_index,=struct.unpack('i', sw.read(4))
				polarisation_pair=str(sw.read(2))
				pulsar_bin,=struct.unpack('i', sw.read(4))
				data_weight,=struct.unpack('d', sw.read(8))
				U,=struct.unpack('d', sw.read(8))
				V,=struct.unpack('d', sw.read(8))
				W,=struct.unpack('d', sw.read(8))
				if baseline_number==bn and polarisation_pair==pp and freq_index==fi:
					count+=1
				headnumber+=1
				headcount+=1

			headcount=0
			headnumber=0
			data=np.zeros((count,64), dtype=np.complex)
			while headcount<size/period:
				sw.seek(period*headnumber)
				sync=str(sw.read(4))
				head_version,=struct.unpack('i', sw.read(4))
				baseline_number,=struct.unpack('i', sw.read(4))
				MJD,=struct.unpack('i', sw.read(4))
				seconds,=struct.unpack('d', sw.read(8))
				config_index,=struct.unpack('i', sw.read(4))
				source_index,=struct.unpack('i', sw.read(4))
				freq_index,=struct.unpack('i', sw.read(4))
				polarisation_pair=str(sw.read(2))
				pulsar_bin,=struct.unpack('i', sw.read(4))
				data_weight,=struct.unpack('d', sw.read(8))
				U,=struct.unpack('d', sw.read(8))
				V,=struct.unpack('d', sw.read(8))
				W,=struct.unpack('d', sw.read(8))
				if baseline_number==bn and polarisation_pair==pp and freq_index==fi:
					for y in range (64):
						data[x,y]= complex(struct.unpack('f', sw.read(4))[0], 
							struct.unpack('f', sw.read(4))[0])
					x+=1
				headcount+=1
				headnumber+=1
				print count
			
		else:
			data=np.zeros((per,64), dtype=np.complex)
			while x<per:
				sw.seek(period*headnumber)
				sync=str(sw.read(4))
				head_version,=struct.unpack('i', sw.read(4))
				baseline_number,=struct.unpack('i', sw.read(4))
				MJD,=struct.unpack('i', sw.read(4))
				seconds,=struct.unpack('d', sw.read(8))
				config_index,=struct.unpack('i', sw.read(4))
				source_index,=struct.unpack('i', sw.read(4))
				freq_index,=struct.unpack('i', sw.read(4))
				polarisation_pair=str(sw.read(2))
				pulsar_bin,=struct.unpack('i', sw.read(4))
				data_weight,=struct.unpack('d', sw.read(8))
				U,=struct.unpack('d', sw.read(8))
				V,=struct.unpack('d', sw.read(8))
				W,=struct.unpack('d', sw.read(8))
				if baseline_number==bn and polarisation_pair==pp and freq_index==fi:
					print '\niteration\n', x
					for y in range (64):
						data[x,y]= complex(struct.unpack('f', sw.read(4))[0], 
							struct.unpack('f', sw.read(4))[0])
					x+=1
				headcount+=1
				headnumber+=1

	np.save(outfile, data)
	print data.shape
	print data

def plot(fn='', razr=''):
	data=np.load(fn)
	average=data.shape[0]/razr
	new_data=np.zeros((razr,64), dtype=np.complex)
	for x1 in range (razr):
			for x in range (average):
				new_data[x1,0:]+=data[x1+razr*x,0:]
	t=abs(fft.fftshift(fft.ifft(new_data)))
	np.save('outfile', t)
	x = np.arange (0, 64, 1)
	y = np.arange (0, razr, 1)
	xgrid, ygrid = np.meshgrid(x, y)
	zgrid = t
	x, y, z = xgrid,ygrid,zgrid
	fig = pylab.figure()
	axes = Axes3D(fig)
	axes.plot_surface(x, y, z)
	plt.xlabel('delay')
	plt.ylabel('accumulation periods')
	pylab.show()
	

def plot2D(s=''):
	t=np.load('outfile.npy')
	srez=t[s,0:]
	print srez
	plt.plot(srez)
	plt.show()

def snr_u(s=''):
	t=np.load('outfile.npy')
	srez=t[s,0:]
	srez=srez.tolist()
	A=max(srez)
	index_from=srez.index(A)-5
	index_to=srez.index(A)+5
	del srez[index_from:index_to]
	print len(srez)
	average=np.mean(srez)
	sko=np.std(srez)
	SNR=(A-average)/sko
	return SNR, average, sko



	'''
A=max(srez)
	index_a = np.where(np.in1d(srez, [A]))[0]
	index_from=index_a[0]-5
	index_to=index_a[0]+5
	new_srez=
	for elem_number in range (-5,6): #from a to b-1
		index=index_a[0]+elem_number
		print index
		del srez[index]
	average=np.mean(srez)
	sko=std(srez)
	SNR=(A-average)/sko
	return SNR

		headcount=0
		headnumber=0
		x=0
		while headcount<size/period:
			sw.seek(period*headnumber)
			sync=str(sw.read(4))
			head_version,=struct.unpack('i', sw.read(4))
			baseline_number,=struct.unpack('i', sw.read(4))
			MJD,=struct.unpack('i', sw.read(4))
			seconds,=struct.unpack('d', sw.read(8))
			config_index,=struct.unpack('i', sw.read(4))
			source_index,=struct.unpack('i', sw.read(4))
			freq_index,=struct.unpack('i', sw.read(4))
			polarisation_pair=str(sw.read(2))
			pulsar_bin,=struct.unpack('i', sw.read(4))
			data_weight,=struct.unpack('d', sw.read(8))
			U,=struct.unpack('d', sw.read(8))
			V,=struct.unpack('d', sw.read(8))
			W,=struct.unpack('d', sw.read(8))
			if baseline_number==bn and polarisation_pair==pp and freq_index==fi:
				for y in range (64):
					data[x,y]= complex(struct.unpack('f', sw.read(4))[0], 
						struct.unpack('f', sw.read(4))[0])
				x+=1
			headcount+=1
			headnumber+=1
	for elem in data:
		if elem.any==0:
			print 'error'
			'''