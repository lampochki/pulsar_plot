import struct
import os
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import pylab
from mpl_toolkits.mplot3d import Axes3D


#pulsar_plot.create('DIFX_56341_040423.s0000.b0000', 258, 'RR', 0, 'filname',0)

def create(fn='', bn='', pp='', fi='', outfile='', per=''):
	size=os.path.getsize(fn)
	current_head=0
	head_count=0
	accum_periods=0
	x=0
	#measure the distance between sync words
	with open(fn, 'rb') as sw:
		sync=int(struct.unpack('i', sw.read(4))[0])
		sync_search=0 
		sync_distance=0
		while sync_search!=sync:
			sync_distance+=1
			sw.seek(sync_distance)
			sync_search=int(struct.unpack('i', sw.read(4))[0])
		#all accumulation periods, count for array shape	
		if per==0:
			while head_count<size/sync_distance:
				sw.seek(sync_distance*current_head)
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
					accum_periods+=1
				current_head+=1
				head_count+=1
			head_count=0
			current_head=0
			data=np.zeros((accum_periods,64), dtype=np.complex)
			while head_count<size/sync_distance:
				sw.seek(sync_distance*current_head)
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
				head_count+=1
				current_head+=1
			print 'all accumulation periods number:', accum_periods
		#specified number of periods	
		if per>0:
			data=np.zeros((per,64), dtype=np.complex)
			while x<per:
				sw.seek(sync_distance*current_head)
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
				head_count+=1
				current_head+=1
	np.save(outfile, data)
	print data.shape
	print data

def snr_1D(fn='', per='', aver='', sect=''):
	data=np.load(fn)[per,0:]
	average=data.shape[0]/aver
	aver_data=np.zeros((aver,64), dtype=np.complex)
	for x1 in range (aver):
			for x in range (average):
				aver_data[x1,0:]+=data[x1+aver*x,0:]
	ifft_data=abs(fft.fftshift(fft.ifft(aver_data)))
	section=ifft_data[sect,0:]
	section=section.tolist()
	A=max(section)
	index_from=section.index(A)-5
	index_to=section.index(A)+5
	del section[index_from:index_to]
	average=np.mean(section)
	st_deviation=np.std(section)
	SNR=(A-average)/st_deviation
	return SNR

	

def plot(fn='', razr=''):
	np.save('outfile', ifft_data)
	x = np.arange (0, 64, 1)
	y = np.arange (0, razr, 1)
	xgrid, ygrid = np.meshgrid(x, y)
	zgrid = ifft_data
	x, y, z = xgrid,ygrid,zgrid
	fig = pylab.figure()
	axes = Axes3D(fig)
	axes.plot_surface(x, y, z)
	plt.xlabel('delay')
	plt.ylabel('accumulation periods')
	pylab.show()
	

def plot2D(s=''):
	ifft_data=np.load('outfile.npy')
	section=ifft_data[s,0:]
	print section
	plt.plot(section)
	plt.show()


	

