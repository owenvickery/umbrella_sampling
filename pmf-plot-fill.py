#!/usr/bin/env python

import os, sys 
import numpy as np
from pylab import *
from matplotlib.font_manager import FontProperties
import argparse
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime


def get_file(file):
	try:
		file_out=np.genfromtxt(file, autostrip=True, comments='@',skip_header=13)
	except:
		sys.exit("Cannot find pull file")
	print('\ngetting file from: '+file)
	return file_out
def ask_direction():
	order=['+','-', '']
	while True:
		direction = input("\nPlease enter direction (+, - or enter for blank):  ")
		try:
			if direction in order:
				break
		except KeyboardInterrupt:
			sys.exit('\nInterrupted')
		except:
			print("Oops!  That was not + or -  Try again...")
	return direction
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
def pull_question(pull_file):
	collective_variable = str(input("\nDo you wish to use a different pull file (press enter to use default):  "))
	if len(collective_variable)>0:
		pull=get_file(collective_variable)
	else:
		pull=get_file(str(args.pull)+pull_file)
		print("\ndefault pull file used")	
	return pull
def filename(location):
	print('\ngetting file names from pull directory')
	pull_multi, tpr_multi, xtc_multi=False, False, False
	for file in os.listdir(location):
		if file.endswith('_pullx.xvg'):
			if pull_multi==False:
				pull_file=file
				pull_multi=True
			else:
				sys.exit('too many pullx files')
		if file.endswith('tpr'):
			if tpr_multi==False:
				tpr_file=file
				tpr_multi=True
			else:
				sys.exit('too many tpr files')
		if file.endswith('xtc'):
			if xtc_multi==False:
				xtc_file=file
				xtc_multi=True
			else:
				sys.exit('too many xtc files')
	return pull_file, xtc_file, tpr_file
def gromacs(cmd):
	print('\nrunning gromacs: \n '+cmd)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err, out = output.communicate()
	exitcode = output.returncode
	out=out.decode("utf-8")
	checks = open('gromacs_outputs'+'_'+timeStamp, 'a')
	checks.write(out)
def final(cmd):
	print('\ncollecting final reaction coordinate')
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err, out = output.communicate()
	err=err.decode("utf-8")
	err.splitlines()
	return err.splitlines()
def create_additional(additional, initial_offset, offset, direction):
	if additional ==True:
		parameters=['md.mdp', 'output directory', 'topology file', 'index file', 'min.mdp']
		arguments=[args.mdp, args.o, args.p, args.n, args.min]
		if None in arguments:
			for i in range(len(parameters)):
				print('-'+parameters[i]+'\t'+str(arguments[i]))
		else:
			minimise(initial_offset, offset, direction)
			make_umbrellas(initial_offset, offset, direction)
	return final('awk \'/Pull group  natoms  pbc atom/{nr[NR+2]}; NR in nr\' gromacs_outputs_'+timeStamp+' | awk \'{print $4}\'')
def setup():
	print('\ninitialising setup')
	react_coord_proposed, react_coord_init=[],[]
	direction=ask_direction()
	offset = ask_integer('Please enter window offset [eg 12]:  ')
	initial_offset=offset
	offset, proposed, actual, pull_check = get_conformation(args.s, args.e, args.int, offset, direction, True)
	for i in range(len(proposed)):
		react_coord_proposed.append(proposed[i])
		react_coord_init.append(actual[i])
	additional = ask_yes_no('Do you wish to create umbrella tpr files (yes/no):  ')  
	react_coord_final=create_additional(additional, initial_offset, offset, direction)
	return react_coord_final, react_coord_proposed, react_coord_init, additional, direction, initial_offset
def get_histograms():
	print('\ngetting histograms from: '+args.hist)
	try:
		histo=np.genfromtxt(args.hist, autostrip=True, comments='@',skip_header=13) #read in histograms sys.argv[2]
	except:
		sys.exit("Cannot find Histogram file")

	histt, histograms=histo[:,0],histo[:,1:-1] # simplify variables
	hist_sum=(histograms.sum(axis=1))/np.max(histograms.sum(axis=1))
	hist_rel=histograms/np.max(histograms)#.sum(axis=1)
	overlap_cutoff=np.mean(np.max(histo[:,1:-1]))*0.1
	overlap=[]
	for i in range(len(histo[0:,1:-1])):
		overlap.append(np.count_nonzero(histo[0:,1:-1][i] > overlap_cutoff))

	return histt, histograms,hist_sum, overlap_cutoff, overlap, hist_rel
def get_conformation(start, end, interval, offset, direction, pull_check, pull): 
	print('\nsetting up umbrella window coordinates')
	pull_file, xtc_file, tpr_file=filename(args.pull)
	if pull_check==True:
		pull=pull_question(pull_file)
		pull_check=False

	pulltime, pulldist=pull[:,0],pull[:,1]
	frametime, distance, dist =[], [], []
	frame=0
	drange=np.arange(start,end,interval)
	if start==end:
		drange=[start]
	for j in range(len(drange)):
		run=True
		for i in range(frame, len(pulltime)):
			if float(drange[j])-0.01 <= float(pulldist[i]) <= float(drange[j])+0.01:
				if run==True:
					frametime.append(pulltime[i])
					frame=i
					distance.append(pulldist[i])
					dist.append(drange[j])
					run=False
	offsetnew=offset
	drange=drange[:len(frametime)]
	for x in range(len(frametime)):
		gromacs('echo 0 | gmx trjconv -f '+args.pull+'/'+xtc_file+' -s '+args.pull+'/'+tpr_file+' -b '+str(frametime[x])+' -e '+ \
		str(frametime[x])+' -o '+args.window+'/conf_'+direction+str(x+1+offset)+'.pdb'+str(args.extra))
		offsetnew=x+1+offset
	return offsetnew, np.around(drange, decimals=3), np.around(distance, decimals=3), pull_check, pull 
def fill_gaps():
	react_coord_proposed, react_coord_init=[],[]
	print('\nfilling in gaps in PMF')
	histt, histograms,hist_sum, overlap_cutoff, overlap,hist_rel=get_histograms()
	count=0
	start, end, pull=0,0,[]
	done, pull_check, initial, check=False, True,True,True
	direction=ask_direction()
	offset = ask_integer('Please enter window offset [eg 12]:  ')
	initial_offset=offset
	for i in range(0, len(histt)):
		if overlap[i] < 3 or hist_sum[i] <= np.mean(hist_sum)*0.25:
			if direction=='-':
				colvar=histt[i]*-1
			if direction=='+' or direction == '':
				colvar=histt[i]
			if colvar >= 0:
				count+=1
				if check == True :   ### for initialisation of script
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
					if done==True or i == len(histt)-1:
						if count>=3:
							if direction == '+' or direction == '':
								offset, proposed, actual, pull_check, pull =get_conformation(start, end, args.int, offset, direction, pull_check, pull )
								for i in range(len(proposed)):
									react_coord_proposed.append(proposed[i])
									react_coord_init.append(actual[i])
							else:
								offset, proposed, actual, pull_check, pull =get_conformation(end,start, args.int, offset, direction, pull_check, pull )
								for i in range(len(proposed)):
									react_coord_proposed.append(proposed[i])
									react_coord_init.append(actual[i])
						elif count ==2:
							offset, proposed, actual, pull_check, pull =get_conformation(start, start, args.int, offset, direction, pull_check, pull )
							for i in range(len(proposed)):
								react_coord_proposed.append(proposed[i])
								react_coord_init.append(actual[i])
						count=1
						done=False
						initial=i
						start=colvar
	additional = ask_yes_no('Do you wish to create umbrella tpr files (yes/no):  ')
	react_coord_final=create_additional(additional, initial_offset, offset, direction)
	return react_coord_final, react_coord_proposed, react_coord_init, additional, direction, initial_offset
def minimise(offset_initial, offset, direction):
	print('\nmaking minimised windows')
	for i in range(offset_initial+1, offset+1):
		direct=False
		try: 
			os.makedirs(args.window+'/min/r'+direction+str(i))
			direct=True
		except:
			additional=ask_yes_no(args.window+'/min/r'+str(i)+'\t already exists \nDo you wish add to folder anyway (yes/no):  ')
			if additional ==True:
				gromacs('gmx grompp -f '+args.min+' -p '+args.p+' -n '+args.n+' -maxwarn 1 -c '+args.window+'/conf_'+direction+str(i)+'.pdb -o '+args.window+'/min/r'+direction+str(i)+'/window_'+direction+str(i))
		if direct==True:
			gromacs('gmx grompp -f '+args.min+' -p '+args.p+' -n '+args.n+' -maxwarn 1 -c '+args.window+'/conf_'+direction+str(i)+'.pdb -o '+args.window+'/min/r'+direction+str(i)+'/window_'+direction+str(i))
	cwd=os.getcwd()
	for i in range(offset_initial+1, offset+1):
		os.chdir(args.window+'/min/r'+direction+str(i))
		gromacs('gmx mdrun -v -deffnm window_'+direction+str(i))
		os.chdir(cwd)
		
def make_umbrellas(offset_initial, offset, direction):
	print('\nmaking umbrellas windows')
	for i in range(offset_initial+1, offset+1):
		direct=False
		try: 
			os.makedirs(args.o+'/r'+direction+str(i))
			direct=True
		except:
			additional=ask_yes_no(args.o+'/r'+str(i)+'\t already exists \nDo you wish overwrite (yes/no):  ')
			if additional ==True:
				gromacs('gmx grompp -f '+args.mdp+' -p '+args.p+' -n '+args.n+' -maxwarn 1 -c '+args.window+'/min/r'+direction+str(i)+'/window_'+direction+str(i)+'.gro -o '+args.o+'/r'+direction+str(i)+'/window_'+direction+str(i))
		if direct==True:
			gromacs('gmx grompp -f '+args.mdp+' -p '+args.p+' -n '+args.n+' -maxwarn 1 -c '+args.window+'/min/r'+direction+str(i)+'/window_'+direction+str(i)+'.gro -o '+args.o+'/r'+direction+str(i)+'/window_'+direction+str(i))
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
	pmf_current=get_file(args.pmf)

	if np.isfinite(pmf_current[:,1]).any()==False:
		sys.exit('Your data is NaN, you most likely have a empty histogram')
	pmf_current[:,1]=set_to_zero(pmf_current[:,1])

	other_pmf = ask_yes_no('Do you wish to plot another pmf:  ')
	if other_pmf==True:
		extra_pmf_file = str(input("\nwhich pmf do wish to plot extra:  "))
		extra_pmf=get_file(extra_pmf_file)
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
	#energy landscape
	subplot(4,1,1)
	title('PMF',  fontproperties=font1, fontsize=15,y=1)
	if other_pmf==True:
		plot(extra_pmf[:,0],extra_pmf[:,1], linewidth=3, color='blue')
		fill_between(extra_pmf[:,0],extra_pmf[:,1]-extra_pmf[:,2], extra_pmf[:,1]-extra_pmf[:,2], alpha=0.3, facecolor='green')	
	plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, color='red')	
	fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);yticks(np.arange(-500,500,step_y), fontproperties=font1, fontsize=15)
	ylim(min_y, max_y);xlim(min_x,max_x)
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right='off', top='off')
	xlabel('Distance(nm)', fontproperties=font2,fontsize=15);ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=15) 

	# Sum of histograms
	subplot(4,1,2)
	title('Normalised histogram sum',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,hist_sum, linewidth=3, color='black')
	plot([-100,100],[arb_cutoff,arb_cutoff], linewidth=3, color='red')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=15)#np.arange(0,np.max(hist_sum), 1000000), 
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right='off', top='off')
	xlabel('Distance(nm)', fontproperties=font2,fontsize=15);ylabel('Sum of counts', fontproperties=font2,fontsize=15) 

	# Histogram overlap
	subplot(4,1,3)
	title('Histogram overlap',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,overlap, linewidth=3, color='black')
	plot([-100,100],[3,3], linewidth=3, color='red')
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,np.max(overlap), 2), fontproperties=font1, fontsize=15)
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right='off', top='off')
	xlabel('Distance(nm)', fontproperties=font2,fontsize=15);ylabel('Overlap', fontproperties=font2,fontsize=15) 

	# Histograms
	subplot(4,1,4)
	title('Histogram windows',  fontproperties=font1, fontsize=15,y=1)
	plot(histt,hist_rel)
	xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);xlim(min_x,max_x)
	yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=15)#np.arange(0,np.max(hist_sum), 1000000), 
	tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right='off', top='off')
	xlabel('Distance(nm)', fontproperties=font2,fontsize=15);ylabel('Counts', fontproperties=font2,fontsize=15) 

	subplots_adjust(left=0.15, wspace=0.4, hspace=0.8, top=0.95, bottom=0.1)

	show()
def results(miscs):
	react_coord_final, react_coord_proposed, react_coord_init,tpr, direction, initial_offset=miscs
	result_write = open('collective_variable_position'+'_'+timeStamp, 'w')
	if tpr==True:
		if len(react_coord_final)==len(react_coord_proposed):
			line=str('proposed       selected         final         window')
			print(line)
			result_write.write(line+'\n')
		else:
			sys.exit('mismatch between tpr output and expected output, check gromacs_outputs')
	elif tpr==False:
		line=str('proposed       selected         window')      
		print(line)
		result_write.write(line+'\n')
	for i in range(len(react_coord_proposed)):
		if tpr==True:
			line=str('  '+str(react_coord_proposed[i])+'\t\t'+str(react_coord_init[i])+'\t\t'+str(react_coord_final[i])+'\t\t'+direction+str(initial_offset+1+i))
			print(line)
			result_write.write(line+'\n')
		elif tpr==False:
			line=str('  '+str(react_coord_proposed[i])+'\t\t'+str(react_coord_init[i])+'\t\t'+direction+str(initial_offset+1+i))	
			print(line)
			result_write.write(line+'\n')
def run_wham():
	core = str(ask_integer('what core would you like to run this on 0-11: '))
	gromacs('taskset --cpu-list '+core+' gmx wham -if '+args.en+' -it '+args.tpr+' -hist '+args.hist+' -bsres '+args.pmf+' -nBootstrap '+ \
	str(args.boot)+' -b '+str(args.b)+' -o '+args.profile+' -bsprof '+args.bsprof+' -temp '+str(args.temp)+' '+str(args.extra[0]))

parser = argparse.ArgumentParser()
parser.add_argument('-pmf', help='location of pmf with bootstrap and name if used with wham',metavar='bsres.xvg',type=str)
parser.add_argument('-hist', help='location of histogram and name if used with wham',metavar='histo.xvg',type=str)
parser.add_argument('-pull', help='location of pull directory',metavar='/*/*/*/', type=str)
parser.add_argument('-f', help='what to do initial setup, plot, fill or wham',metavar='[setup, plot, fill, wham]',type=str, choices= ['setup', 'plot', 'fill', 'wham'])
parser.add_argument('-s', help='where to start on reaction coordinate',metavar='0',type=float)
parser.add_argument('-e', help='where to end on reaction coordinate',metavar='5', type=float)
parser.add_argument('-mdp', help='umbrella mdp file',metavar='md.mdp', type=str)
parser.add_argument('-o', help='umbrella output location',metavar='/*/*/*/', type=str)
parser.add_argument('-window', help='umbrella structure file output/location',metavar='/*/*/*/', type=str)
parser.add_argument('-p', help='topology file',metavar='topol.top', type=str)
parser.add_argument('-n', help='index file',metavar='index.ndx', type=str)
parser.add_argument('-int', help='interval for umbrella windows (nm)',metavar='0.05', type=float)
parser.add_argument('-en', help='energy.dat file',metavar='en.dat', type=str)
parser.add_argument('-tpr', help='tpr.dat file',metavar='tpr.dat', type=str)
parser.add_argument('-boot', help='number of bootstraps to do',metavar='200', type=int)
parser.add_argument('-b', help='First time to analyse (ps)',metavar='200', type=int)
parser.add_argument('-profile', help='location of pmf profile',metavar='profile.xvg',type=str)
parser.add_argument('-bsprof', help='location of bootstrap profile',metavar='bsprof.xvg',type=str)
parser.add_argument('-temp', help='simulation temperature',metavar='310',type=int)
parser.add_argument('-extra', help='any extra commands for gromacs, each command should be in \'\' e.g. \'-min 5 -max 6\' ',metavar='-max',type=str, nargs='*')
parser.add_argument('-min', help='umbrella minimise mdp file',metavar='min.mdp', type=str)

args = parser.parse_args()

timeStamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

if args.extra== None:
	args.extra=''
if args.f == 'setup':
	parameters=['pull', 's', 'e', 'window']
	arguments=[args.pull, args.s, args.e , args.window]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		misc = setup()
		results(misc)
elif args.f== 'plot':
	parameters=['pmf', 'hist']
	arguments=[args.pmf, args.hist]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		plot_pmf()
elif args.f == 'fill':
	parameters=['pmf', 'hist', 'pull', 'windows']
	arguments=[args.pmf, args.hist, args.pull, args.window]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		misc = fill_gaps()
		results(misc)
elif args.f == 'wham':
	parameters=['en.dat', 'tpr.dat','pmf', 'hist', 'bootstrap', 'time', 'profile', 'bsprof', 'temperature']
	arguments=[args.en, args.tpr, args.pmf, args.hist, args.boot, args.b, args.profile, args.bsprof, args.temp]
	if None in arguments:
		for i in range(len(parameters)):
			print('-'+parameters[i]+'\t'+str(arguments[i]))
	else:
		run_wham()
		plot_results=ask_yes_no('Do you wish to plot the results? (yes/no)')
		if plot_results==True:
			plot_pmf()
else:
	sys.exit('Try again, enter  [-f setup, plot, fill or wham]')



