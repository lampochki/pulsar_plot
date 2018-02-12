import struct
import os
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import pylab
from mpl_toolkits.mplot3d import Axes3D

#pulsar_plot.create('DIFX_56341_040423.s0000.b0000', 258, 'RR', 0, 'filname',0)

#pulsar_plot.create('DIFX_56341_040423.', 258, 'RR', 0, '2sec258ch0all',0)

#pulsar_plot.create('DIFX_57844_018000.s0000.b0000', 258, 'RR', 0, 'qaz")

def create(fn='', bn='', pp='', fi='', outfile='', per=''):
	size=os.path.getsize(fn)
	current_head=0
	head_count=0
	accum_periods=0
	spec_chan=2048
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
		print 'distance between syncwords', sync_distance
		spec_chan=(sync_distance-74)/8
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
			data=np.zeros((accum_periods,spec_chan), dtype=np.complex)
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
					for y in range (spec_chan):
						data[x,y]= complex(struct.unpack('f', sw.read(4))[0], 
							struct.unpack('f', sw.read(4))[0])
					x+=1
				head_count+=1
				current_head+=1
			print 'all accumulation periods number:', accum_periods
		#specified number of periods	
		if per>0:
			data=np.zeros((per,spec_chan), dtype=np.complex)
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
					for y in range (spec_chan):
						data[x,y]= complex(struct.unpack('f', sw.read(4))[0], 
							struct.unpack('f', sw.read(4))[0])
					x+=1
				head_count+=1
				current_head+=1
	np.save(outfile, data)
	print data.shape
	print data


def plot(fn='', aver='',periods='', outfile=''):
	periods=periods*aver
 	data=np.load(fn)[0:periods]
 	spec_chan=all_data.shape[1]
 	average=data.shape[0]/aver
 	new_data=np.zeros((aver,spec_chan), dtype=np.complex)
 	for x1 in range (aver):
 			for x in range (average):
 				new_data[x1,0:]+=data[x1+aver*x,0:]
 	ifft_data=abs(fft.fftshift(fft.ifft(new_data)))
 	x = np.arange (0, spec_chan, 1)
 	y = np.arange (0, aver, 1)
 	xgrid, ygrid = np.meshgrid(x, y)
 	zgrid = ifft_data
 	x, y, z = xgrid,ygrid,zgrid
 	fig = pylab.figure()
 	axes = Axes3D(fig)
 	axes.plot_surface(x, y, z)
 	plt.xlabel('delay')
 	plt.ylabel('accumulation periods')
 	axes.view_init(0, 90)
 	pylab.show()
 	np.save(outfile, ifft_data)


def plot2D(fn='',s='', outfile=''):
	data=np.load(fn)
	spec_chan=data.shape[1]
	section=data[0:,s]
	print section
	np.save(outfile, section)
	plt.plot(section)
	plt.show()


def plots(fn=''):
	data=np.load(fn)
	spec_chan=data.shape[1]
	window_number=1
	fig=plt.figure()
	for aver in range (120,124):
		average=data.shape[0]/aver
		aver_data=np.zeros((aver,spec_chan), dtype=np.complex)
		for x1 in range (aver):
				for x in range (average):
					aver_data[x1,0:]+=data[x1+aver*x,0:]
		ifft_data=abs(fft.fftshift(fft.ifft(aver_data)))
		ax=fig.add_subplot(2,2,window_number,projection='3d')
		x = np.arange (0, spec_chan, 1)
 		y = np.arange (0, aver, 1)
 		xgrid, ygrid = np.meshgrid(x, y)
 		zgrid = ifft_data
 		x, y, z = xgrid,ygrid,zgrid
 		surf=ax.plot_surface(x,y,z)
 		plt.xlabel('delay')
 		plt.ylabel('accumulation periods')
 		window_number+=1
 		ax.view_init(0, 90)
 	pylab.show()



def snr(fn=''):
	section=np.load(fn)
	section=section.tolist()
	A=max(section)
	index_from=section.index(A)-5
	index_to=section.index(A)+5
	del section[index_from:index_to]
	average=np.mean(section)
	st_deviation=np.std(section)
	SNR=(A-average)/st_deviation
	print SNR, A


#pulsar_plot.snr_1D('filname.npy', 63, 42)
def snr_1D(fn='', aver='', sect=''):
	spec_chan=160
	periods=aver
	SNR_new=1
	SNR_old=0
	SNR=0
	while SNR_old<SNR_new:
		data=np.load(fn)[0:periods,0:]
		average=data.shape[0]/aver
		aver_data=np.zeros((aver,spec_chan), dtype=np.complex)
		for x1 in range (aver):
				for x in range (average):
					aver_data[x1,0:]+=data[x1+aver*x,0:]
		ifft_data=abs(fft.fftshift(fft.ifft(aver_data)))
		section=ifft_data[0:,sect]
		section=section.tolist()
		A=max(section)
		index_from=section.index(A)-5
		index_to=section.index(A)+5
		del section[index_from:index_to]
		average=np.mean(section)
		st_deviation=np.std(section)
		SNR_old=SNR
		SNR=(A-average)/st_deviation
		SNR_new=SNR
		periods+=aver
		print SNR, 
	return SNR_old


def snr_all(fn='', aver='', sect=''):
	periods=aver
	all_data=np.load(fn)
	spec_chan=all_data.shape[1]
	while all_data.shape[0]-periods>=aver:
		data=np.load(fn)[0:periods,0:]
		average=data.shape[0]/aver
		aver_data=np.zeros((aver,spec_chan), dtype=np.complex)
		for x1 in range (aver):
				for x in range (average):
					aver_data[x1,0:]+=data[x1+aver*x,0:]
		ifft_data=abs(fft.fftshift(fft.ifft(aver_data)))
		section=ifft_data[0:,sect]
		section=section.tolist()
		A=max(section)
		index_from=section.index(A)-5
		index_to=section.index(A)+5
		del section[index_from:index_to]
		average=np.mean(section)
		st_deviation=np.std(section)
		SNR=(A-average)/st_deviation
		per=periods/63
		print'periods=',per 
		print 'snr=' ,SNR
		periods+=aver
		

def find_max(fn=''):
	peaks=[]
	data=np.load(fn)
	spec_chan=data.shape[1]
	datan=data.tolist()
	for aver in range (60,70):
		average=data.shape[0]/aver
		aver_data=np.zeros((aver,spec_chan), dtype=np.complex)
		for x1 in range (aver):
				for x in range (average):
					aver_data[x1,0:]+=data[x1+aver*x,0:]
		ifft_data=abs(fft.fftshift(fft.ifft(aver_data)))
		ifft_data=ifft_data.tolist()
		peak=max(max(ifft_data))
		peaks.append(peak)
	peak=max(peaks)
	peak_index=peaks.index(peak)+60
	print peak, peak_index

