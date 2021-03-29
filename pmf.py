#!/usr/bin/env python3

import os, sys 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import argparse
import subprocess
from time import gmtime, strftime
import time
from shutil import copyfile
import distutils.spawn
import multiprocessing as mp
import scipy.signal as signal
from shutil import copyfile


###############################
os.environ['GMX_SUPPRESS_DUMP'] = '1'

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

#flags to change if needed

#location=os.getcwd()+'/umbrella_windows'  #  use line below if you want to change location to a absolute path

#location = 'xxx/xxx/xxx/xxx'

#directories=[location,location+'/frames', location+'/minimised', location+'/windows',location+'/analysis',location+'/setup_files_'+timestamp]

###############################
def ask_number(question):
    while True:
        try:
            number = float(input("\n"+question))
            break
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except:
            print("Oops!  That was not a number.  Try again...")
    return number

### checks whether all arguements are present ###
def check_arguments(arguments):

    for argue in arguments:
        if argue not in options or options[argue] == None:
            print('\narguement: \''+argue+'\' is missing')
            failed = True
    if 'failed' in locals():
        print('\n')
        parser.print_help()
        sys.exit()


### finds gromacs installation
def find_gromacs():
    if args.gmx != None:
        args.gmx=distutils.spawn.find_executable(args.gmx)
    else:
        args.gmx=distutils.spawn.find_executable('gmx')
    if args.gmx is None or type(args.gmx) != str:
        if os.environ.get("GMXBIN") != None:
            for root, dirs, files in os.walk(os.environ.get("GMXBIN")):
                for file_name in files:
                    if file_name.startswith('gmx') and file_name.islower() and '.' not in file_name:
                        args.gmx=distutils.spawn.find_executable(file_name)
                        if type(args.gmx) == str and args.gmx != None :
                            break
                        else:
                            args.gmx=None
                break
        if args.gmx is None:
            sys.exit('Cannot find gromacs installation')

def file_copy_and_check(file_in,file_out):
    if not os.path.exists(file_out) and os.path.exists(file_in):
        copyfile(file_in, file_out)


## gets window offset
def find_offset():
    if os.path.exists(frames) and args.offset == 0:
        pdb_frame = []
        for file in os.listdir(frames):
            if file.endswith('.pdb'):
                pdb_frame.append(int(file[7:-4]))
        if len(pdb_frame) == 0:
            args.offset=0
        else:
            args.offset=max(pdb_frame)

### gets CV from pull file 
def get_pull():
    file_out = [[],[]]
    if os.path.exists(args.pull):
        with open(args.pull, 'r') as pull_input:
            for line_nr, line in enumerate(pull_input.readlines()):
                if len(line)>0:
                    if line[0] not in ['@', '#']:
                        if len(line.split())==2:
                            file_out[0].append(float(line.split()[0]))
                            file_out[1].append(float(line.split()[1]))
        return file_out
    else:
        sys.exit("cannot find pull file")
### make minimise mdp file
def make_min():
    if not os.path.exists(running_dir+'em.mdp'):
        em = open(running_dir+'/em.mdp','w')
        em.write('integrator = steep\nnsteps     = 5000\nemtol      = 1000\nemstep     = 0.001')

### make file directories for umbrella sampling
def folders():
    for directory in [running_dir, frames, minimised_frames, umbrella_windows, analysis, setup_files]:
        if not os.path.exists(directory):
            os.mkdir(directory)
        else:
            print(directory+' already exists')

### copies input files to backup
def backup():
    if args.tpr:
        copyfile(args.p, setup_files+'/topol.top')
        copyfile(args.mdp, setup_files+'/md.mdp')
        copyfile(args.n, setup_files+'/index.ndx')

### runs gromacs commands
def gromacs(gro):
    cmd = gro[0]
    if args.v >= 1:
        print('\nrunning gromacs: \n '+cmd+'\n')
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err, out = output.communicate()
    out=out.decode("utf-8")
    if args.func != 'wham':
        checks = open(setup_files+'/gromacs_outputs'+'_'+timestamp, 'a')
    else:
        checks = open('gromacs_outputs'+'_'+timestamp, 'a')
    checks.write(out)
    #### standard catch for failed gromacs commands
    if 'File input/output error:' in out:
        sys.exit('\n'+out)
    elif 'Error in user input:' in out:
        sys.exit('\n'+out)
    elif 'did not converge to Fmax ' in out:
        sys.exit('\n'+out)
    elif 'Segmentation fault (core dumped):' in out:
        sys.exit('\n'+out)
    elif 'Fatal error:' in out:
        sys.exit('\n'+out)
    elif 'but did not reach the requested Fmax' in out:
        sys.exit('\n'+out)
    elif 'number of atoms in the topology (' in out:
        sys.exit('\n'+out+'\n\nIf it is only 2 atoms out check cysteine distances, and increase -cys cutoff')

    if len(gro) == 3: 
        gro[2].put(gro[1])
        return gro[1]
### collects final CV position from grompp
def final(cmd):
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = output.communicate()
    out=out.decode("utf-8")
    out.splitlines()
    return out.splitlines()

### sets up pmf windows
def setup():
    print('\ninitialising setup')
    pull=get_pull()                     ##  gets CV coordinates
    start_pull = min(pull[1]) if args.start is None else args.start              ##  selection of start point
    end_pull = max(pull[1]) if args.end is None else args.end

    print('\ncreating umbrella windows between\n')
    print('{0:^10}{1:^10}{2:^10}\n{3:^10}{4:^10}{5:^10}\n'.format('start', 'end','interval', np.round(start_pull, 2), end_pull, args.int))   
    end, react_coord_proposed, react_coord_init = get_conformation(start_pull, end_pull, args.offset, pull,[], []) ## gets conformations from pull trajectory
    react_coord_final=equilibrate(end)         ##  equilibrates umbrella windows 
    results([react_coord_proposed, react_coord_init, react_coord_final])        ##  returns CVs

def report_complete(func, size, resid):
    complete = np.round((size/resid)*100,2)
    if complete.is_integer():
        print('Running '+func+' on '+str(resid)+' molecules. Percentage complete: ',complete,'%', end='\r')
    time.sleep(0.1)


def equilibrate(offset_end):
    pool = mp.Pool(mp.cpu_count())
    m = mp.Manager()
    q = m.Queue()
    if not args.tpronly:
        if args.min:
            make_min()
            for i in range(args.offset+1, offset_end+1):     ##  makes individual window folders
                if not os.path.exists(minimised_frames+'window_'+str(i)):
                    os.makedirs(minimised_frames+'window_'+str(i))
            pool_process = pool.map_async(gromacs, [(args.gmx+' grompp -po '+setup_files+'/em_out-'+str(i)+
                                    ' -f '+running_dir+'/em.mdp -p '+args.p+' -n '+args.n+' -maxwarn 2 -c '+
                                    frames+'/window_'+str(i)+'.pdb -o '+minimised_frames+'window_'+str(i)+'/window_'+str(i), i, q) 
                                    for i in range(args.offset+1, offset_end+1)])         ## minimisation grompp parallised
            while not pool_process.ready():
                report_complete('GROMPP', q.qsize(), offset_end-args.offset)
            print('{:<130}'.format(''), end='\r')
            print('GROMPP completed\n')
            pool.close
            cwd=os.getcwd()
            for i in range(args.offset+1, offset_end+1):             ##  changes to minimised directory and runs mdrun 
                os.chdir(minimised_frames+'/window_'+str(i))
                gromacs([args.gmx+' mdrun -v -nt 10 -deffnm window_'+str(i)])    
                print('Running mdrun on '+str(args.offset+1)+'-'+str(offset_end+1)+' frames: percentage complete: ',np.round(((i-args.offset)/(offset_end-args.offset))*100,2),'%', end='\r')      
                os.chdir(cwd)
            print('{:<130}'.format(''), end='\r')
            print('Minimisation completed\n')
    if args.tpr:
        make_umbrella_windows(offset_end)

    ### sorts out frame ordering
    window = final('grep windows/window_ '+setup_files+'/gromacs_outputs_'+timestamp+' | grep grompp | awk \'{print substr($0,length($0),1)} \'')
    CV = final('awk \'/Pull group  natoms  pbc atom/{nr[NR+2]}; NR in nr\' '+setup_files+'/gromacs_outputs_'+timestamp+' | awk \'{print $4}\'')
    final_CV=np.stack((window, CV),axis=-1)
    final_CV=np.sort(final_CV, axis=0)
    return final_CV[:,1]

def make_umbrella_windows(offset_end):
    for i in range(args.offset+1, offset_end+1):
        if not os.path.exists(umbrella_windows+'/window_'+str(i)):
            os.makedirs(umbrella_windows+'/window_'+str(i))
        else:
            print('windows folder exists') 

    pool = mp.Pool(mp.cpu_count()) 
    m = mp.Manager()
    q = m.Queue()       
    if args.min:
        pool_process = pool.map_async(gromacs, [(args.gmx+' grompp -po '+setup_files+'/md_out-'+str(i)+' -f '+args.mdp+' -p '+args.p+' -n '+args.n\
        +' -maxwarn 2 -c '+minimised_frames+'window_'+str(i)+'/window_'+str(i)+'.gro -r '+minimised_frames+'/window_'+str(i)+'/window_'+str(i)+'.gro -o '\
        +umbrella_windows+'window_'+str(i)+'/window_'+str(i), i, q) for i in range(args.offset+1, offset_end+1)])           ## makes umbrella windows from minimised frames
    else:   
        pool_process = pool.map_async(gromacs, [(args.gmx+' grompp -po '+setup_files+'/md_out-'+str(i)+' -f '+args.mdp+' -p '+args.p+' -n '+args.n\
        +' -maxwarn 2 -c '+frames+'window_'+str(i)+'.pdb -r '+frames+'/window_'+str(i)+'.pdb  -o '+umbrella_windows+'/window_'\
        +str(i)+'window_'+str(i), i, q) for i in range(args.offset+1, offset_end+1)])                 ##  ## makes umbrella windows from non-minimised frames
    while not pool_process.ready():
        report_complete('GROMPP', q.qsize(), offset_end-args.offset)
        print('                                                                       ', end='\r')
    pool.close()

### gets frames from pull 
def get_conformation(start, end, offset, pull,react_coord_proposed, react_coord_init): 
    pool = mp.Pool(mp.cpu_count())
    m = mp.Manager()
    q = m.Queue()
    # print('\nsetting up umbrella window coordinates')
    frametime, distance, dist =[], [], []
    if start==end:                          ##  if only 1 frame
        drange=[start]
    else:
        drange = np.arange(start,end,args.int) if start<=end else np.arange(end,start,args.int)

    for CV in drange:                   ##  gets frametime and CV closest to suggested value
            if min(pull[1]) <= CV <= max(pull[1]): 
                react_coord_init.append(np.around(min(pull[1], key=lambda x:abs(x-CV)), decimals=3))
                frametime.append(pull[0][pull[1].index(min(pull[1], key=lambda x:abs(x-CV)))])
                react_coord_proposed.append(np.around(CV, decimals=3))
    if not args.tpronly:   ### runs trjconv to get frames uses multiprocessing
        pool_process = pool.map_async(gromacs, [('echo 0 | '+args.gmx+' trjconv -f '+args.f+' -s '+args.s+' -b '+str(ftime)+' -e '+str(ftime)+' -o '+frames+'window_'\
            +str(x+1+offset)+'.pdb', x, q) for x,ftime in enumerate(frametime)])
        while not pool_process.ready(): 
            report_complete('trjconv', q.qsize(), len(frametime))
        pool.close()
    offsetnew=len(frametime)+offset
    return offsetnew, react_coord_proposed, react_coord_init

### reads in histograms
def get_histograms():
    print('\ngetting histograms from: '+args.hist+'\n')
    if os.path.exists(args.hist):
        histograms=np.genfromtxt(args.hist, autostrip=True, comments='@',skip_header=13)
    else:
        sys.exit("Cannot find Histogram file")
    hist_sum=(histograms[:,1:].sum(axis=1))/np.max(histograms[:,1:].sum(axis=1))        ##  gives the sum of histograms relative to max sum
    hist_rel=histograms[:,1:]/np.max(histograms[:,1:])         
    ##  makes histograms relative to maximum
    overlap_cutoff=np.mean(np.max(hist_rel[:,1:], axis=0))*0.05
    overlap=[]
    for i in range(len(histograms[0:,1:])):  
        overlap.append(np.count_nonzero(hist_rel[i] > overlap_cutoff))  ##  provides overlap count if there is a overlap above cutoff
    skewed_data = check_histogram_normality(histograms)
    return histograms[:,0], histograms[:,1:],hist_sum, overlap_cutoff, np.array(overlap), hist_rel, skewed_data

def check_histogram_normality(data):
    skewed_data = []
    for histogram in range(1,len(data[0,:])):
        non_zero = np.where(data[:,histogram] > 0)
        series = data[:,histogram][non_zero]
        indexes, _ = signal.find_peaks(series)#, distance=(np.max(series)-np.min(series))/4)
        if len(indexes) > 1 and np.all(series[indexes]>np.max(series)*0.1):
            if 'header' not in locals():
                print('These windows may need to be repeated with a higher force constant\n')
                header = True
            print('multiple peaks found in window:',histogram)
            skewed_data.append(histogram)
    print()
    return np.array(skewed_data)

def preprocess_arrays(pull, coord, overlap, histogram_sum):
    available_coord = np.where(np.logical_and(coord>np.min(pull[1]),coord<np.max(pull[1])))
    coord_cut = coord[available_coord]
    overlap_cut = overlap[available_coord]
    histogram_sum_cut = histogram_sum[available_coord]

    coord_fill = coord_cut[np.where(np.logical_or(overlap_cut < args.cutoff, histogram_sum_cut <= np.mean(histogram_sum)*0.25))]
    if args.start != 0:
        coord_fill = coord_fill[np.where(coord_fill > args.start)]
    if args.end != None:
        coord_fill = coord_fill[np.where(coord_fill < args.end)]
    return coord_fill

### fills in gaps from the histograms
def fill_gaps():
    print('\nfilling in gaps in PMF')
    coord, histograms,histogram_sum, overlap_cutoff, overlap, histogram_rel, skewed_data=get_histograms()        ##  gets histogram data
    pull=get_pull()                                                                                 ##  gets pull information
    react_coord_proposed, react_coord_init=[],[]
    initial_offset, offset=args.offset, args.offset
    coord_cut = preprocess_arrays(pull, coord, overlap, histogram_sum)
    cv_start=True
    for i, cv in enumerate(coord_cut):
        if cv_start==True:
            cv_initial=cv 
            cv_end=cv
            cv_start=False
        if cv_initial-args.int <= cv <= cv_initial+args.int:
            cv_end=cv
        else:
            cv_start=True
            offset, react_coord_proposed, react_coord_init = get_conformation(cv_initial, cv_end, offset, pull, react_coord_proposed, react_coord_init)
    if cv_start==False:
        offset, react_coord_proposed, react_coord_init = get_conformation(cv_initial, cv_end, offset, pull, react_coord_proposed, react_coord_init)
    results([react_coord_proposed, react_coord_init,equilibrate(offset)])

def set_to_zero(energy):
    if energy[-1] < 0:
        energy = energy+(0-energy[-1])
    else:
        energy = energy-energy[-1]  
    return energy

def get_min_max(pmf_current, min_x=0, max_x=0, min_y=0, max_y=0):
    if min_x > np.round(min(pmf_current[:,0]), 2)-0.25 or min_x ==0:
        min_x = np.round(min(pmf_current[:,0]), 2)-0.25
    if max_x < np.round(max(pmf_current[:,0]), 2)+0.25 or max_x ==0:
        max_x = np.round(max(pmf_current[:,0]), 2)+0.25

    if args.y_range != None:
        min_y, max_y = args.y_range[0], args.y_range[1]
    else:
        if min_y > np.round(min(pmf_current[:,1]), 2)-10 or min_y ==0:
            min_y = np.round(min(pmf_current[:,1]), 2)-10
        if max_y < np.round(max(pmf_current[:,1]), 2)+10 or max_y ==0:
            max_y = np.round(max(pmf_current[:,1]), 2)+10
    return min_x, max_x, min_y, max_y

def plot_pmf():
    print('\nplotting PMF data')
    # Fonts
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
    plt.rcParams['mathtext.default'] = 'regular'

    ### histogram start###
    histt, histograms,hist_sum, overlap_cutoff, overlap, hist_rel, skewed_data=get_histograms()
    ### histogram end###


    energy_time=[]              ## collects energy minima for timecourse
### plotting header
    plt.rcParams['axes.linewidth']=3
    if len(skewed_data) > 0:
        plt.figure(1, figsize=(15,25))
        add_to_plot = 1
    else:
        plt.figure(1, figsize=(15,20))
        add_to_plot = 0
    # #energy landscape
    plt.subplot(4+add_to_plot,1,1)
    plt.title('PMF',  fontproperties=font1, fontsize=25,y=1.18)

    for pmf_number, pmf in enumerate(args.pmf):
        ### pmf plot start ###
        try:
            pmf_current=np.genfromtxt(pmf, autostrip=True, comments='@',skip_header=13)
        except:
            sys.exit("Cannot find pmf file")

        if np.isfinite(pmf_current[:,1]).any()==False:
            sys.exit('Your data is NaN, you most likely have a empty histogram')
        pmf_current[:,1]=set_to_zero(pmf_current[:,1])

    ### pmf plot end ###
        ### min max start ###
        if pmf_number == 0:
            min_x, max_x, min_y, max_y = get_min_max(pmf_current)

        ### min max end ###
            energy_minima=np.round(min(pmf_current[:,1]), 2)
            energy_minima_error=0
            print('pmf '+str(pmf_number+1)+': energy minima = '+str(energy_minima)+' +/- '+ 
                  str(np.round( np.mean(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)))
            plt.plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, color='red', label='PMF 1')    
            if len(pmf_current[0]) == 3:
                plt.fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')
                energy_minima_error=np.round( np.mean(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)
            energy_time.append([energy_minima, energy_minima_error])
        else:
            min_x, max_x, min_y, max_y = get_min_max(pmf_current, min_x, max_x, min_y, max_y)
            energy_minima=np.round(min(pmf_current[:,1]), 2)
            energy_minima_error=0
            print('pmf '+str(pmf_number+1)+': energy minima = '+str(energy_minima)+' +/- '+ 
                  str(np.round( np.mean(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)))
            plt.plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, label='PMF '+str(pmf_number+1))#, color='blue')
            if len(pmf_current[0]) == 3:
                plt.fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')  
                energy_minima_error=np.round( np.mean(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)
            energy_time.append([energy_minima, energy_minima_error])
    plt.xticks(np.arange(-500,500,args.xtick), fontproperties=font1, fontsize=25)
    if args.ytick != None:
        plt.yticks(np.arange(-500,500,args.ytick), fontproperties=font1, fontsize=25)
    else:
        plt.yticks(fontproperties=font1, fontsize=25)

    plt.ylim(min_y, max_y);plt.xlim(min_x,max_x)
    plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)
    plt.xlabel('Distance (nm)', fontproperties=font2,fontsize=25);plt.ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=25) 
    plt.legend(prop={'size': 25}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0.)#len(args.pmf)
    # Sum of histograms
    plt.subplot(4+add_to_plot,1,2)
    plt.title('Normalised histogram sum',  fontproperties=font1, fontsize=25,y=1)
    plt.plot(histt,hist_sum, linewidth=3, color='black')
    plt.plot([-100,100],[0.25, 0.25], linewidth=3, color='red')
    plt.xticks(np.arange(-500,500,args.xtick), fontproperties=font1, fontsize=25);plt.xlim(min_x,max_x)
    plt.yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=25)#np.arange(0,np.max(hist_sum), 1000000), 
    plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)
    plt.xlabel('Distance (nm)', fontproperties=font2,fontsize=25);plt.ylabel('Sum of counts', fontproperties=font2,fontsize=25) 

    # Histogram overlap
    plt.subplot(4+add_to_plot,1,3)
    plt.title('Histogram overlap',  fontproperties=font1, fontsize=25,y=1)
    plt.plot(histt,overlap, linewidth=3, color='black')
    plt.plot([-100,100],[args.cutoff,args.cutoff], linewidth=3, color='red')
    plt.xticks(np.arange(-500,500,args.xtick), fontproperties=font1, fontsize=25);plt.xlim(min_x,max_x)
    plt.yticks(np.arange(0,np.max(overlap)+2, 2), fontproperties=font1, fontsize=25)
    plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)
    plt.xlabel('Distance (nm)', fontproperties=font2,fontsize=25);plt.ylabel('Overlap', fontproperties=font2,fontsize=25) 

    # Histograms
    plt.gca().set_prop_cycle(None)
    plt.subplot(4+add_to_plot,1,4)
    plt.gca().set_prop_cycle(None)
    plt.title('Histogram windows',  fontproperties=font1, fontsize=25,y=1)
    plt.plot(histt,hist_rel)
    plt.plot([-100,100],[overlap_cutoff, overlap_cutoff], linewidth=3, color='red')

    plt.xticks(np.arange(-500,500,args.xtick), fontproperties=font1, fontsize=25);plt.xlim(min_x,max_x)
    plt.yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=25)
    plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)
    plt.xlabel('Distance (nm)', fontproperties=font2,fontsize=25);plt.ylabel('Counts', fontproperties=font2,fontsize=25) 

    if len(skewed_data) > 0:
        plt.gca().set_prop_cycle(None)
        plt.subplot(4+add_to_plot,1,5)
        
        plt.title('Skewed Histogram windows',  fontproperties=font1, fontsize=25,y=1)
        plt.plot(histt,hist_rel, alpha=0.15, linewidth=2)
        plt.plot(histt,hist_rel[:,skewed_data-1], linewidth=3)
        plt.xticks(np.arange(-500,500,args.xtick), fontproperties=font1, fontsize=25);plt.xlim(min_x,max_x)
        plt.yticks(np.arange(0,1.01,0.25),fontproperties=font1, fontsize=25)
        plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)
        plt.xlabel('Distance (nm)', fontproperties=font2,fontsize=25);plt.ylabel('Counts', fontproperties=font2,fontsize=25)         
        plt.subplots_adjust(left=0.15, wspace=0.4, hspace=0.4, top=0.95, bottom=0.1)
    else:
        plt.subplots_adjust(left=0.15, wspace=0.4, hspace=0.4, top=0.95, bottom=0.1)
    if args.save_plot == None:
        plt.savefig('energy_landscape_'+timestamp+'.png', dpi=500)
    else:
        plt.savefig(args.save_plot+'.png', dpi=500)


    if len(energy_time)>1:
        plt.figure(2, figsize=(10,10))
        energy_time=np.array(energy_time)
        # if len(energy_time) > 1 and args.timestep != 1:
        #     args.timestep=ask_number('what is the timestep? ')
        # #energy landscape
        plt.title('PMF timecourse',  fontproperties=font1, fontsize=25,y=1.1)
        plt.errorbar(np.arange(0+args.timestep,len(energy_time[:,0])*args.timestep+args.timestep,args.timestep),energy_time[:,0], yerr=energy_time[:,1], color='red', zorder=1)
        plt.plot(np.arange(0+args.timestep,len(energy_time[:,0])*args.timestep+args.timestep,args.timestep),energy_time[:,0],'o',color='k') 
        plt.ylim(np.min(energy_time[:,0])-5, np.max(energy_time[:,0])+5);plt.xlim(0,len(energy_time[:,0])*args.timestep+(2*args.timestep))
        plt.tick_params(axis='both', which='major', width=3, length=5, labelsize=25, direction='in', pad=10, right=False, top=False)

        plt.xlabel('time (ns)', fontproperties=font2,fontsize=25);plt.ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=25) 
        plt.subplots_adjust(left=0.18, right=0.98, wspace=0.4, hspace=0.4, top=0.8, bottom=0.1)

        if args.save_conv == None:
            plt.savefig('energy_landscape_convergence_'+timestamp+'.png', dpi=500)
        else:
            plt.savefig(args.save_conv+'.png', dpi=500)
    if args.show:
        plt.show()

def results(miscs):
    react_coord_proposed, react_coord_init,react_coord_final=miscs
    result_write = open(setup_files+'/collective_variable_position'+'_'+timestamp, 'w')
    if args.tpr:
        if len(react_coord_final)==len(react_coord_proposed):
            line='\n{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}'.format('proposed', 'selected', 'final', 'S-P', 'F-S', 'window') 
            print(line)
            result_write.write(line+'\n')
        else:
            sys.exit('mismatch between tpr output and expected output, check gromacs_outputs')
    else:
        line='{0:^10}{1:^10}{2:^10}{3:^10}'.format('proposed','selected','S-P','window')   
        print(line)
        result_write.write(line+'\n')
    for i in range(len(react_coord_proposed)):
        if args.tpr:
            try:
                line='{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}'.format(react_coord_proposed[i], react_coord_init[i], react_coord_final[i], 
                            np.around(float(react_coord_init[i])-float(react_coord_proposed[i]), decimals=3), 
                            np.around(float(react_coord_final[i])-float(react_coord_init[i]), decimals=3 ), args.offset+1+i)
                print(line)
            except ValueError:
                line=str(args.offset+1+i)+'   error in final grompp most likely'
            result_write.write(line+'\n')
        else:
            try:
                line='{0:^10}{1:^10}{2:^10}{3:^10}'.format(react_coord_proposed[i],react_coord_init[i], 
                            np.around(float(react_coord_init[i])-float(react_coord_proposed[i]), decimals=3), 
                            args.offset+1+i)
                print(line)
            except ValueError:
                line=str(args.offset+1+i)+'   error in final grompp most likely'
            result_write.write(line+'\n')
    backup()


def pull_concat(window):
        if args.current:
            xvg_loc = args.win_dir
        else:
            xvg_loc = args.win_dir+'/window_'+str(window)+'/'

        if os.path.exists(xvg_loc):
            files_range = get_file_range(xvg_loc, window)
            if len(files_range) == 1:
                os.system('cp '+xvg_loc+files_range[0]+' window_'+str(window)+'_pullf_com.xvg')     
                return [window, len(files_range)]
            elif len(files_range) >= 2: 
                time, pull = read_in_pullf_files(xvg_loc,files_range)
                if len(time) > 0:
                    write_concat_pullf(window, time, pull)
                os.system('cp '+xvg_loc+'window_'+str(window)+'.tpr'+' window_'+str(window)+'.tpr')
                return [window, len(files_range)]                
            else:
                return [window, 'MISSING']
        else:
            return [window, 'MISSING']

def get_file_range(xvg_loc, window):
    files_range = []
    for filename in os.listdir(xvg_loc):
        if filename.endswith('.xvg') and 'pullf' in filename:
            if 'window_'+str(window)+'.' in filename or 'window_'+str(window)+'_' in filename and not '_com' in filename:
                files_range.append(filename)
    return files_range

def read_in_pullf_files(xvg_loc,files_range):
    time, pull = [], []
    for part in range(len(files_range)):
        for line in open(xvg_loc+files_range[part], 'r').readlines():
            if not line[0] in ['#', '@']:
                if len(line) > 0 and len(line.split()) == 2:
                    try:
                        time.append(float(line.split()[0]))
                        pull.append(float(line.split()[1]))
                    except:
                        break
    return time, pull 

def write_concat_pullf(window, time, pull):
    time_ord, pull_ord = (list(t) for t in zip(*sorted(zip(time, pull))))
    em = open('window_'+str(window)+'_pullf_com.xvg','w')
    for j in range(len(time_ord)):
        em.write(str(time_ord[j])+'\t'+str(pull_ord[j])+'\n')

def run_wham():
    print('Now running wham\n')
    if int(args.core) > mp.cpu_count():
        args.core = str(int(ask_number('what core would you like to run this on 0-11: ')))
    if len(args.pmf) == 1:
        args.pmf = args.pmf[0]
    else:
        sys.exit('Please input a single name for flag -pmf.')
    if args.te != '0':
        trim = ' -e '+args.te
    else:
        trim = ''
    gromacs(['taskset --cpu-list '+args.core+' '+args.gmx+' wham -if '+args.en_list+' -it '+args.tpr_list+' -temp '+args.temp+' -hist '+args.hist+' -bsres '+args.pmf+' -nBootstrap '+args.boot+' -bins '+str(args.bins)+' -b '+args.ts+trim])


if __name__ == '__main__':

    start = time.time()

    parser = argparse.ArgumentParser(description='Sets up and analyses PMF calculations', prog='pmf', epilog='Enjoy the program and best of luck!\n')

    
    parser.add_argument('-func', help='what to do initial setup, plot, concat, wham or fill',metavar='[setup, plot, concat, wham, fill]',type=str, choices= ['setup', 'plot', 'concat','wham', 'fill'], required=True)
    parser.add_argument('-f', help='xtc file for setup (eg your pull xtc)',metavar='pull.xtc',type=str)
    parser.add_argument('-mdp', help='umbrella mdp file',metavar='md.mdp', type=str)
    parser.add_argument('-s', help='structure file for setup (use the pull.tpr)',metavar='pull.tpr',type=str)
    parser.add_argument('-n', help='index file for the system',metavar='index.ndx', type=str)
    parser.add_argument('-p', help='topology file',metavar='topol.top', type=str)
    parser.add_argument('-pull', help='pull file for setup',metavar='pullx.xvg',type=str)
    parser.add_argument('-offset', help='window offset',metavar='5',type=int, default=0)
    parser.add_argument('-tpr', help='do not make tpr files', action='store_false')
    parser.add_argument('-min', help='switch off minisation', action='store_false')
    parser.add_argument('-int', help='interval for umbrella windows (nm)',metavar='0.05', type=float, default=0.05)
    parser.add_argument('-start', help='where to start on reaction coordinate',metavar='0',type=float, default=0)
    parser.add_argument('-end', help='where to end on reaction coordinate',metavar='5', type=float)
    parser.add_argument('-ws', help='windows start',metavar='1',type=int)
    parser.add_argument('-we', help='windows end',metavar='5', type=int)
    parser.add_argument('-ts', help='time start (equilibration)',metavar='5', type=str, default='0')
    parser.add_argument('-te', help='time end (equilibration)',metavar='5', type=str, default='0')
    parser.add_argument('-boot', help='number of bootstraps to run',metavar='5', type=str, default='1')
    parser.add_argument('-pmf', help='location of pmf ',metavar='bsres.xvg',type=str, nargs='*', default=['bsres.xvg'])
    parser.add_argument('-hist', help='location of histogram and name if used with wham',metavar='histo.xvg',type=str, default='histo.xvg')
    parser.add_argument('-en_list', help='energy file list',metavar='en.dat',type=str, default='en.dat')
    parser.add_argument('-tpr_list', help='tpr file list',metavar='tpr.dat',type=str, default='tpr.dat')
    parser.add_argument('-temp', help='temperature',metavar='310',type=str, default='310')
    parser.add_argument('-tpronly', help='only makes tpr files default (False) requires energy minimised files', action='store_true')
    parser.add_argument('-current', help='to concat in current directory', action='store_true')
    parser.add_argument('-cutoff', help='histogram overlap cutoff for filling and plotting', default=3, type=int)
    parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels) (Optional)")
    parser.add_argument('-loc', help='output folder name, (default = umbrella_windows)',metavar='umbrella_windows',type=str, default='umbrella_windows')
    parser.add_argument('-gmx', help='gromacs executable name (Optional)',metavar='gmx_avx',type=str)
    parser.add_argument('-core', help='which core to run wham on',metavar='1', type=str, default='0')
    parser.add_argument('-bins', help='number of bins for wham',metavar='1', type=str, default='200')
    parser.add_argument('-save_plot', help='plot save name',metavar='pmf_plot', type=str)
    parser.add_argument('-save_conv', help='convergence save name',metavar='pmf_conv', type=str)
    parser.add_argument('-ytick', help='pmf Y-axis tick interval', type=int)
    parser.add_argument('-xtick', help='pmf Y-axis tick interval', type=int, default=1)
    parser.add_argument('-y_range', help='range of Y-axis ',type=float, nargs=2)
    parser.add_argument('-show', help='show interactive plots', action='store_true')
    parser.add_argument('-timestep', help='Convergence timestep', type=int, default=1 )
    parser.add_argument('-win_dir', help='Window directory',metavar='../windows/', type=str, default='../windows/')


    args = parser.parse_args()
    options = vars(args)
    start_dir        = os.getcwd()+'/'  ### initial working directory
    running_dir      = os.getcwd()+'/'+args.loc
    frames           = running_dir+'/frames/'   ### working directory 
    minimised_frames = running_dir+'/minimised/'  ### final directory for run files
    umbrella_windows = running_dir+'/windows/'  ### contains input run files
    analysis         = running_dir+'/analysis/'  ### contains run files
    setup_files      = running_dir+'/setup_files_'+timestamp




    if args.func not in ['plot','concat']:
        find_gromacs()
    find_offset()
    


    tpr_arg =[]
    if args.func == 'setup':
        if args.tpr:
            tpr_arg = ['mdp', 'tpr','f', 'n', 'p', 's']
        check_arguments(['pull', 'int']+tpr_arg)
        folders()
        setup()
    elif args.func == 'fill':
        if args.tpr:
            tpr_arg = ['mdp', 'tpr','f', 'n', 'p', 's']
        check_arguments(['pull', 'int', 'offset']+tpr_arg)
        folders()
        fill_gaps()
    elif args.func== 'plot':
        check_arguments(['pmf', 'hist'])
        plot_pmf()
    elif args.func== 'concat':
        print('Concatonating files in the range of: ', args.ws,'-', args.we)
        print('Within the location: ', args.win_dir)
        pool = mp.Pool(mp.cpu_count())
        check_arguments(['ws', 'we'])
        concatonated = pool.map(pull_concat, [(window) for window in range(args.ws, args.we+1)])          ## makes umbrella windows from minimised frames
        print('\nwindow\ttotal part numbers')
        np.array(concatonated).sort(axis=0)
        for line in concatonated:
            if line[1] != 0:
                print(line[0], '\t', line[1])
            else:
                print(line[0], '\t', 'SKIPPED')

    elif args.func== 'wham':
        check_arguments(['pmf', 'boot', 'ts', 'te', 'core'])
        run_wham()

    end = time.time()

    print('\nThis script took: '+str(np.round(end-start,1))+'s to run')