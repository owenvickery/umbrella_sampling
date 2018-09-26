#!/usr/bin/env python

import os, sys 
import numpy as np
from pylab import *
from matplotlib.font_manager import FontProperties
import argparse
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
from shutil import copyfile

def ask_integer(question):
	while True:
		try:
			integer = int(input("\n"+question))
			break
		except KeyboardInterrupt:
			sys.exit('\nInterrupted')
		except:
			print("Oops!  That was not a number.  Try again...")
	return integer

def ask_yes_no(question):
	additional=False
	while True:
		q_result = input("\n"+question)
		try:
			if q_result.lower() in ['yes','y']:
				additional=True
				break
			elif q_result.lower() in ['n', 'no']:
				break
			else:
				print("\nplease enter yes or no.")
		except KeyboardInterrupt:
			sys.exit('\nInterrupted')
	return additional

def get_pull():
	try:
		file_out=np.genfromtxt(args.pull, autostrip=True, comments='@',skip_header=13)
	except:
		sys.exit("Cannot find pull file")
	pull=[[],[]]
	for i in range(len(file_out[:,0])):
		pull[0].append(file_out[:,0][i])
		pull[1].append(file_out[:,1][i])
	return pull

def make_min():
	if not os.path.exists(location+'em.mdp'):
		em = open(location+'/em.mdp','w')
		em.write('integrator = cg\nnsteps     = 5000\nemtol      = 1000\nemstep     = 0.001')

def folders():
	for i in range(len(directories)): 
		try: 
			os.makedirs(directories[i])
		except:
			print(directories[i]+' folder exists')
def backup():
	if args.tpr == True:
		copyfile(args.p, location+'/setup_files_'+timestamp+'/topol.top')
		copyfile(args.mdp, location+'/setup_files_'+timestamp+'/md.mdp')
		copyfile(args.n, location+'/setup_files_'+timestamp+'/index.ndx')


def gromacs(cmd):
	print('\nrunning gromacs: \n '+cmd)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err, out = output.communicate()
	exitcode = output.returncode
	out=out.decode("utf-8")
	checks = open(location+'/setup_files_'+timestamp+'/gromacs_outputs'+'_'+timestamp, 'a')
	checks.write(out)

def final(cmd):
	print('\ncollecting final reaction coordinate')
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = output.communicate()
	out=out.decode("utf-8")
	out.splitlines()
	return out.splitlines()

def setup():
	print('\ninitialising setup')
	folders()


	react_coord_proposed, react_coord_init=[],[]
	pull=get_pull()
	end, proposed, actual = get_conformation(args.start, args.end, args.int, args.offset, pull)
	for i in range(len(proposed)):
		react_coord_proposed.append(proposed[i])
		react_coord_init.append(actual[i])
	react_coord_final=equilibrate(args.offset, end)
	return react_coord_proposed, react_coord_init, react_coord_final 

def equilibrate(offset, offset_end):
	if args.tpr ==True:
		parameters=['md.mdp', 'tpr file','xtc file', 'index file', 'topology file']
		arguments=[args.mdp, args.s, args.f, args.n, args.p]
		make_min()
		if None in arguments:
			for i in range(len(parameters)):
				print('-'+parameters[i]+'\t'+str(arguments[i]))
		else:
			print('\nmaking minimised windows')
			for i in range(offset+1, offset_end+1):
				try: 
					os.makedirs(directories[2]+'/window_'+sign+str(i))
				except:
					print('minise folder exists')
				gromacs(gmx+' grompp -po '+location+'/setup_files_'+timestamp+'/em_out-'+str(i)+' -f '+location+'/em.mdp -p '+args.p+' -n '+args.n+' -maxwarn 2 -c '+directories[1]+'/window_'+sign+str(i)+'.pdb -o '+directories[2]+'/window_'+sign+str(i)+'/window_'+sign+str(i))
			cwd=os.getcwd()
			for i in range(offset+1, offset_end+1):
				os.chdir(directories[2]+'/window_'+sign+str(i))
				gromacs(gmx+' mdrun -v -deffnm window_'+sign+str(i))
				os.chdir(cwd)
			print('\nmaking umbrellas windows')
			for i in range(offset+1, offset_end+1):
				try: 
					os.makedirs(directories[3]+'/window_'+sign+str(i))
				except:
					print('windows folder exists')
				gromacs(gmx+' grompp -po '+location+'/setup_files_'+timestamp+'/md_out-'+str(i)+' -f '+args.mdp+' -p '+args.p+' -n '+args.n+' -maxwarn 2 -c '+directories[2]+'/window_'+sign+str(i)+'/window_'+sign+str(i)+'.gro -o '+directories[3]+'/window_'+sign+str(i)+'/window_'+sign+str(i))
	return final('awk \'/Pull group  natoms  pbc atom/{nr[NR+2]}; NR in nr\' '+location+'/setup_files_'+timestamp+'/gromacs_outputs_'+timestamp+' | awk \'{print $4}\'')

def get_conformation(start, end, interval, offset, pull): 
	print('\nsetting up umbrella window coordinates')
	frametime, distance, dist =[], [], []
	drange=np.arange(start,end,interval)
	if start==end:
		drange=[start]
	for j in range(len(drange)):
		distance.append(min(pull[1], key=lambda x:abs(x-drange[j])))
		frametime.append(pull[0][pull[1].index(min(pull[1], key=lambda x:abs(x-drange[j])))])
		dist.append(drange[j])
	offsetnew=offset
	for x in range(len(frametime)):
		gromacs('echo 0 | '+gmx+' trjconv -f '+args.f+' -s '+args.s+' -b '+str(frametime[x])+' -e '+str(frametime[x])+' -o umbrella_windows/frames/window_'+sign+str(x+1+offset)+'.pdb')
		offsetnew=x+1+offset
	return offsetnew, np.around(dist, decimals=3), np.around(distance, decimals=3)

def get_histograms():
	print('\ngetting histograms from: '+args.hist)
	try:
		histograms=np.genfromtxt(args.hist, autostrip=True, comments='@',skip_header=13)
	except:
		sys.exit("Cannot find Histogram file")
	hist_sum=(histograms[:,1:-1].sum(axis=1))/np.max(histograms[:,1:-1].sum(axis=1))
	hist_rel=histograms[:,1:-1]/np.max(histograms[:,1:-1])
	overlap_cutoff=np.mean(np.max(hist_rel))*0.1 ################################################### 10% of total hisgram height discarded
	overlap=[]
	for i in range(len(histograms[0:,1:-1])):
		overlap.append(np.count_nonzero(hist_rel[i] > overlap_cutoff))
	return histograms[:,0], histograms[:,1:-1],hist_sum, overlap_cutoff, overlap, hist_rel

def fill_gaps():
	print('\nfilling in gaps in PMF')
	folders()
	coord, histograms,histogram_sum, overlap_cutoff, overlap, histogram_rel=get_histograms()
	pull=get_pull()
	react_coord_proposed, react_coord_init=[],[]
	count=0
	start, end=0,0
	done, pull_check, initial, check=False, True,True,True
	initial_offset, offset=args.offset, args.offset
	for i in range(0, len(coord)):
		if overlap[i] < 3 or histogram_sum[i] <= np.mean(histogram_sum)*0.25:
			if args.dir==True:
				colvar=coord[i]*-1
			else:
				colvar=coord[i]
			if colvar >= 0:
				count+=1
				if check == True : 
					start=colvar
					initial=i
					check=False
				else:
					if i != initial+1 and count>=2:
						done=True					
					if i == initial+1 and count>=2:
						end=colvar
						initial=i
					if i != initial+1 and count==1:
						initial=i
					if done==True or i == len(coord)-1:
						if count>=3:
							if args.dir==True:
								offset, proposed, actual = get_conformation( end, start, args.int, offset, pull )
							else:
								offset, proposed, actual = get_conformation(start, end, args.int, offset, pull )
							
							for i in range(len(proposed)):
								react_coord_proposed.append(proposed[i])
								react_coord_init.append(actual[i])
						elif count ==2:
							offset, proposed, actual = get_conformation(start, start, args.int, offset, pull )
							for i in range(len(proposed)):
								react_coord_proposed.append(proposed[i])
								react_coord_init.append(actual[i])
						count=1
						done=False
						initial=i
						start=colvar
	react_coord_final=equilibrate(args.offset, offset)
	return react_coord_proposed, react_coord_init, react_coord_final

def set_to_zero(energy):
	if energy[-1] < 0:
		energy = energy+(0-energy[-1])
	else:
		energy = energy-energy[-1]	
	return energy
def plot_pmf():
	print('\nplotting PMF data')
	# Fonts
	alignment = {'horizontalalignment': 'center', 'verticalalignment': 'baseline'}
	families = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
	styles = ['normal', 'italic', 'oblique']
	font = FontProperties()

	#  tick fonts
	font1 = font.copy()
	font1.set_family('sans-serif')
	font1.set_style('normal')

	# label fonts
	font2 = font.copy()
	font2.set_family('sans-serif')
	font2.set_style('normal')
	font2.set_weight('normal')
	rcParams['mathtext.default'] = 'regular'

	### histogram start###
	histt, histograms,hist_sum, overlap_cutoff, overlap, hist_rel=get_histograms()
	arb_cutoff=np.mean(hist_sum)*0.25
	### histogram end###
	### pmf start ###
	try:
		pmf_current=np.genfromtxt(args.pmf, autostrip=True, comments='@',skip_header=13)
	except:
		sys.exit("Cannot find pmf file")

	if np.isfinite(pmf_current[:,1]).any()==False:
		sys.exit('Your data is NaN, you most likely have a empty histogram')
	pmf_current[:,1]=set_to_zero(pmf_current[:,1])
	other_pmf=False
	other_pmf = ask_yes_no('Do you wish to overlay another pmf:  ')
	if other_pmf==True:
		extra_pmf_file = str(input("\nwhich pmf do wish to plot:  "))
		try:
			extra_pmf=np.genfromtxt(extra_pmf_file, autostrip=True, comments='@',skip_header=13)
		except:
			sys.exit("Cannot find pmf file")
		if np.isfinite(pmf_current[:,1]).any()==False:
			sys.exit('Your data is NaN, you most likely have a empty histogram in your set')
		extra_pmf[:,1]=set_to_zero(extra_pmf[:,1])
	### pmf end ###

	step_y= ask_integer('PMF tick interval length on the Y axis [eg 10]:  ') ## y axis tick interval

	### min max start ###
	min_x, max_x, step_x=np.round(min(pmf_current[:,0]), 2)-0.25, np.round(max(pmf_current[:,0]), 2)+0.25, 1
	try:
		min_y, max_y = [int(x) for x in input("\nmin and max Y (press enter to use defaults) :  ").split()]
		print(min_y, max_y)
	except:
		print("\ndefault values used")
		min_y, max_y= np.round(min(pmf_current[:,1]), 2)-10, np.round(max(pmf_current[:,1]), 2)+10
		print(min_y, max_y)
	### min max end ###

	rcParams['axes.linewidth']=3
	figure(1, figsize=(10,10))
	# #energy landscape
	subplot(4,1,1)
	title('PMF',  fontproperties=font1, fontsize=15,y=1)
	if other_pmf==True:
		plot(extra_pmf[:,0],extra_pmf[:,1], linewidth=3, color='blue')
		if len(extra_pmf[0]) == 3:
			fill_between(extra_pmf[:,0],extra_pmf[:,1]-extra_pmf[:,2], extra_pmf[:,1]-extra_pmf[:,2], alpha=0.3, facecolor='green')	
	plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, color='red')	
	# print(len(pmf_current[0]))
	if len(pmf_current[0]) == 3:
		fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);yticks(np.arange(-500,500,step_y), fontproperties=font1, fontsize=15)
	ylim(min_y, max_y);xlim(min_x,max_x)
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)
	xlabel('Distance (nm)', fontproperties=font2,fontsize=15);ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=15) 

	# Sum of histograms
	subplot(4,1,2)
	title('Normalised histogram sum',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,hist_sum, linewidth=3, color='black')
	plot([-100,100],[arb_cutoff,arb_cutoff], linewidth=3, color='red')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=15)#np.arange(0,np.max(hist_sum), 1000000), 
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)
	xlabel('Distance (nm)', fontproperties=font2,fontsize=15);ylabel('Sum of counts', fontproperties=font2,fontsize=15) 

	# Histogram overlap
	subplot(4,1,3)
	title('Histogram overlap',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,overlap, linewidth=3, color='black')
	plot([-100,100],[3,3], linewidth=3, color='red')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,np.max(overlap), 2), fontproperties=font1, fontsize=15)
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)
	xlabel('Distance (nm)', fontproperties=font2,fontsize=15);ylabel('Overlap', fontproperties=font2,fontsize=15) 

	# Histograms
	subplot(4,1,4)
	title('Histogram windows',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,hist_rel)
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=15)
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)
	xlabel('Distance (nm)', fontproperties=font2,fontsize=15);ylabel('Counts', fontproperties=font2,fontsize=15) 

	subplots_adjust(left=0.15, wspace=0.4, hspace=0.8, top=0.95, bottom=0.1)

	show()
def results(miscs):
	react_coord_proposed, react_coord_init,react_coord_final=miscs
	result_write = open(location+'/setup_files_'+timestamp+'/collective_variable_position'+'_'+timestamp, 'w')
	if args.tpr==True:
		if len(react_coord_final)==len(react_coord_proposed):
			line=str('proposed       selected         final           S-P              F-S         window')
			print(line)
			result_write.write(line+'\n')
		else:
			sys.exit('mismatch between tpr output and expected output, check gromacs_outputs')
	elif args.tpr==False:
		line=str('proposed       selected           S-P         window')      
		print(line)
		result_write.write(line+'\n')
	for i in range(len(react_coord_proposed)):
		if args.tpr==True:
			try:
				line=str('  '+str(react_coord_proposed[i])+'\t\t'+str(react_coord_init[i])+'\t\t'+str(react_coord_final[i])+'\t\t'+str(np.around(float(react_coord_init[i])-float(react_coord_proposed[i]), decimals=3 ))+ \
				'\t\t'+str(np.around(float(react_coord_final[i])-float(react_coord_init[i]), decimals=3 ))+'\t\t'+str(args.offset+1+i))
				print(line)
			except ValueError:
				line=str(args.offset+1+i)+'   error in final grompp most likely'
			result_write.write(line+'\n')
		elif args.tpr==False:
			try:
				line=str('  '+str(react_coord_proposed[i])+'\t\t'+str(react_coord_init[i])+'\t\t'+str(np.around(float(react_coord_init[i])-float(react_coord_proposed[i]), decimals=3 ))+'\t\t'+str(args.offset+1+i))	
				print(line)
			except ValueError:
				line=str(args.offset+1+i)+'   error in final grompp most likely'
			result_write.write(line+'\n')
	backup()
			
# def run_wham():
	# core = str(ask_integer('what core would you like to run this on 0-11: '))
	# gromacs('taskset --cpu-list '+core+' gmx wham -if '+args.en+' -it '+args.tpr+' -hist '+args.hist+' -bsres '+args.pmf+' -nBootstrap '+ \
	# str(args.boot)+' -b '+str(args.b)+' -o '+args.profile+' -bsprof '+args.bsprof+' -temp '+str(args.temp)+' '+str(args.extra[0]))

parser = argparse.ArgumentParser()

parser.add_argument('-mdp', help='umbrella mdp file',metavar='md.mdp', type=str)
parser.add_argument('-func', help='what to do initial setup, plot or fill',metavar='[setup, plot, fill]',type=str, choices= ['setup', 'plot', 'fill'])
parser.add_argument('-f', help='xtc file for setup',metavar='pull.xtc',type=str)
parser.add_argument('-s', help='structure file for setup',metavar='pull.tpr',type=str)
parser.add_argument('-n', help='index file',metavar='index.ndx', type=str)
parser.add_argument('-p', help='topology file',metavar='topol.top', type=str)
parser.add_argument('-pull', help='pull file for setup',metavar='pullx.xvg',type=str)
parser.add_argument('-offset', help='window offset',metavar='5',type=int)
parser.add_argument('-tpr', help='make tpr files default (True)', action='store_true')
parser.add_argument('-int', help='interval for umbrella windows (nm)',metavar='0.05', type=float)
parser.add_argument('-start', help='where to start on reaction coordinate',metavar='0',type=float)
parser.add_argument('-end', help='where to end on reaction coordinate',metavar='5', type=float)
parser.add_argument('-pmf', help='location of pmf ',metavar='bsres.xvg',type=str)
parser.add_argument('-hist', help='location of histogram and name if used with wham',metavar='histo.xvg',type=str)
parser.add_argument('-dir', help='direction default (positive)', action='store_true')
args = parser.parse_args()

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())


#flags to change if needed

gmx='gmx'
location=os.getcwd()+'/umbrella_windows'  #  use line below if you want to change location to a absolute path
#location = 'xxx/xxx/xxx/xxx'
directories=[location,location+'/frames', location+'/minimised', location+'/windows',location+'/analysis',location+'/setup_files_'+timestamp]



if args.dir==True:
	sign='-'
else:
	sign=''

if args.func == 'setup':
	parameters=['pull file', 'start', 'end', 'interval','offset']
	arguments=[args.pull, args.start, args.end, args.int, args.offset]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		misc = setup()
		results(misc)
elif args.func == 'fill':
	parameters=['pull file', 'interval','offset']
	arguments=[args.pull, args.int, args.offset]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		misc = fill_gaps()
		results(misc)
elif args.func== 'plot':
	parameters=['pmf', 'hist']
	arguments=[args.pmf, args.hist]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		plot_pmf()

# elif args.f == 'wham':
# 	parameters=['en.dat', 'tpr.dat','pmf', 'hist', 'bootstrap', 'time', 'profile', 'bsprof', 'temperature']
# 	arguments=[args.en, args.tpr, args.pmf, args.hist, args.boot, args.b, args.profile, args.bsprof, args.temp]
# 	if None in arguments:
# 		for i in range(len(parameters)):
# 			print('-'+parameters[i]+'\t'+str(arguments[i]))
# 	else:
# 		run_wham()
# 		plot_results=ask_yes_no('Do you wish to plot the results? (yes/no)')
# 		if plot_results==True:
# 			plot_pmf()
else:
	sys.exit('Try again, enter  [-f setup, plot, fill]')



