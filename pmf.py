#!/usr/bin/env python3

import os, sys 
import numpy as np
from pylab import *
from matplotlib.font_manager import FontProperties
import argparse
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
from shutil import copyfile
import multiprocessing as mp

###############################

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

#flags to change if needed

gmx='gmx'
location=os.getcwd()+'/umbrella_windows'  #  use line below if you want to change location to a absolute path

#location = 'xxx/xxx/xxx/xxx'

directories=[location,location+'/frames', location+'/minimised', location+'/windows',location+'/analysis',location+'/setup_files_'+timestamp]

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
    correct=True
    for argue in arguments:
        if not bool(options[argue]):
            print('\narguement: \''+argue+'\' is missing')
            correct=False
    if not correct:
        print('\n')
        parser.print_help()
    return correct


### gets CV from pull file 
def get_pull():
    try:
        file_out=np.genfromtxt(args.pull, autostrip=True, comments='@',skip_header=13)
    except:
        print('Cannot find pull file or something is wrong with it')
        print('removing final line to fix (hopefully)')
        try:
            file_out=np.genfromtxt(args.pull, autostrip=True, comments='@',skip_header=13, skip_footer=1)
        except:
            sys.exit("Cannot find pull file or something is wrong with it: "+args.pull)
    pull=[[],[]]
    for i in range(len(file_out[:,0])):
        pull[0].append(file_out[:,0][i]) ## time
        pull[1].append(file_out[:,1][i]) ## CV
    return pull

### make minimise mdp file
def make_min():
    if not os.path.exists(location+'em.mdp'):
        em = open(location+'/em.mdp','w')
        em.write('integrator = steep\nnsteps     = 5000\nemtol      = 1000\nemstep     = 0.001')

### make file directories for umbrella sampling
def folders():
    for i in range(len(directories)): 
        try: 
            os.makedirs(directories[i])
        except:
            print(directories[i]+' folder exists')

### copies input files to backup
def backup():
    if args.tpr:
        copyfile(args.p, location+'/setup_files_'+timestamp+'/topol.top')
        copyfile(args.mdp, location+'/setup_files_'+timestamp+'/md.mdp')
        copyfile(args.n, location+'/setup_files_'+timestamp+'/index.ndx')

### runs gromacs commands
def gromacs(gro):
    cmd = gro[0]
    if args.v >= 1:
        print('\nrunning gromacs: \n '+cmd+'\n')
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err, out = output.communicate()
    exitcode = output.returncode
    out=out.decode("utf-8")
    if args.func != 'wham':
        checks = open(location+'/setup_files_'+timestamp+'/gromacs_outputs'+'_'+timestamp, 'a')
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
    folders()                           ##  makes folder directories
    pull=get_pull()                     ##  gets CV coordinates
    if args.start == None:              ##  selection of start point
        start_pull=min(pull[1])
    else:
        start_pull=args.start
    if args.end == None:                ## selection of end point
        end_pull=max(pull[1])
    else:
        end_pull=args.end
    print('\ncreating umbrella windows between\n')
    print('{0:^10}{1:^10}\n{2:^10}{3:^10}'.format('start', 'end', np.round(start_pull, 2), end_pull))   
    react_coord_proposed, react_coord_init=[],[]
    end, react_coord_proposed, react_coord_init = get_conformation(start_pull, end_pull, args.offset, 
                                                                    pull,react_coord_proposed, react_coord_init) ## gets conformations from pull trajectory
    react_coord_final=equilibrate(args.offset, end)         ##  equilibrates umbrella windows 
    return react_coord_proposed, react_coord_init, react_coord_final        ##  returns CVs

def report_complete(func, size, resid):
        print('Running '+func+' on '+str(resid)+' frames: percentage complete: ',np.round((size/resid)*100,2),'%', end='\r')
        time.sleep(0.1)

def equilibrate(offset, offset_end):
    pool = mp.Pool(mp.cpu_count())
    m = mp.Manager()
    q = m.Queue()
    if not args.tpronly:
        if args.min:
            make_min()
            print('\nmaking minimised windows')
            for i in range(offset+1, offset_end+1):     ##  makes individual window folders
                try: 
                    os.makedirs(directories[2]+'/window_'+str(i))
                except:
                    pass
            pool_process = pool.map_async(gromacs, [(gmx+' grompp -po '+location+'/setup_files_'+timestamp+'/em_out-'+str(i)+
                                    ' -f '+location+'/em.mdp -p '+args.p+' -n '+args.n+' -maxwarn 2 -c '+
                                    directories[1]+'/window_'+str(i)+'.pdb -o '+directories[2]+'/window_'+str(i)+'/window_'+str(i), i, q) 
                                    for i in range(offset+1, offset_end+1)])         ## minimisation grompp parallised
            while not pool_process.ready():
                report_complete('GROMPP', q.qsize(), offset_end-offset)
                print('                                                                       ', end='\r')
                # print('GROMPP completed') 
            pool.close
            cwd=os.getcwd()
            for i in range(offset+1, offset_end+1):             ##  changes to minimised directory and runs mdrun 
                os.chdir(directories[2]+'/window_'+str(i))
                gromacs([gmx+' mdrun -v -nt 10 -deffnm window_'+str(i)])    
                print('Running mdrun on '+str(offset+1)+'-'+str(offset_end+1)+' frames: percentage complete: ',np.round((i/(offset_end-offset))*100,2),'%', end='\r')      
                os.chdir(cwd)
            print('\n')
    if args.tpr:
        make_umbrella_windows(offset, offset_end)

    ### sorts out frame ordering
    window = final('grep /umbrella_windows/windows/window_ '+location+'/setup_files_'+timestamp+'/gromacs_outputs_'+timestamp+' | grep grompp | awk \'{print substr($0,length($0),1)} \'')
    CV = final('awk \'/Pull group  natoms  pbc atom/{nr[NR+2]}; NR in nr\' '+location+'/setup_files_'+timestamp+'/gromacs_outputs_'+timestamp+' | awk \'{print $4}\'')
    final_CV=np.stack((window, CV),axis=-1)
    final_CV=np.sort(final_CV, axis=0)
    return final_CV[:,1]

def make_umbrella_windows(offset, offset_end):
    print('\nmaking umbrellas windows')             ##  makes individual window folders
    for i in range(offset+1, offset_end+1):
        try: 
            os.makedirs(directories[3]+'/window_'+str(i))
        except:
            print('windows folder exists') 
    pool = mp.Pool(mp.cpu_count()) 
    m = mp.Manager()
    q = m.Queue()       
    if args.min:
        pool_process = pool.map_async(gromacs, [(gmx+' grompp -po '+location+'/setup_files_'+timestamp+'/md_out-'+str(i)+' -f '+args.mdp+' -p '+args.p+' -n '+args.n\
        +' -maxwarn 2 -c '+directories[2]+'/window_'+str(i)+'/window_'+str(i)+'.gro -r '+directories[2]+'/window_'+str(i)+'/window_'+str(i)+'.gro -o '\
        +directories[3]+'/window_'+str(i)+'/window_'+str(i), i, q) for i in range(offset+1, offset_end+1)])           ## makes umbrella windows from minimised frames
    else:   
        pool_process = pool.map_async(gromacs, [(gmx+' grompp -po '+location+'/setup_files_'+timestamp+'/md_out-'+str(i)+' -f '+args.mdp+' -p '+args.p+' -n '+args.n\
        +' -maxwarn 2 -c '+directories[1]+'/window_'+str(i)+'.pdb -r '+directories[1]+'/window_'+str(i)+'.pdb  -o '+directories[3]+'/window_'\
        +str(i)+'/window_'+str(i), i, q) for i in range(offset+1, offset_end+1)])                 ##  ## makes umbrella windows from non-minimised frames
    while not pool_process.ready():
        report_complete('GROMPP', q.qsize(), offset_end-offset)
        print('                                                                       ', end='\r')
    pool.close()

### gets frames from pull 
def get_conformation(start, end, offset, pull,react_coord_proposed, react_coord_init): 
    pool = mp.Pool(mp.cpu_count())
    m = mp.Manager()
    q = m.Queue()
    # print('\nsetting up umbrella window coordinates')
    frametime, distance, dist =[], [], []
    if start<=end:
        drange=np.arange(start,end,args.int)    ##  gets umbrella range
    else:
        drange=np.arange(end,start,args.int)
    if start==end:                          ##  if only 1 frame
        drange=[start]
    for CV in drange:                   ##  gets frametime and CV closest to suggested value
            if min(pull[1]) <= CV <= max(pull[1]): 
                react_coord_init.append(np.around(min(pull[1], key=lambda x:abs(x-CV)), decimals=3))
                frametime.append(pull[0][pull[1].index(min(pull[1], key=lambda x:abs(x-CV)))])
                react_coord_proposed.append(np.around(CV, decimals=3))
    if not args.tpronly:   ### runs trjconv to get frames uses multiprocessing
        pool_process = pool.map_async(gromacs, [('echo 0 | '+gmx+' trjconv -f '+args.f+' -s '+args.s+' -b '+str(ftime)+' -e '+str(ftime)+' -o umbrella_windows/frames/window_'\
            +str(x+1+offset)+'.pdb', x, q) for x,ftime in enumerate(frametime)])
        while not pool_process.ready():
            report_complete('trjconv', q.qsize(), len(frametime))
            print('                                                                       ', end='\r')
            # print('trjconv completed')
        pool.close()
    offsetnew=len(frametime)+offset
    return offsetnew, react_coord_proposed, react_coord_init

### reads in histograms
def get_histograms():
    print('\ngetting histograms from: '+args.hist)
    try:
        histograms=np.genfromtxt(args.hist, autostrip=True, comments='@',skip_header=13)
    except:
        sys.exit("Cannot find Histogram file")
    hist_sum=(histograms[:,1:].sum(axis=1))/np.max(histograms[:,1:].sum(axis=1))        ##  gives the sum of histograms relative to max sum
    hist_rel=histograms[:,1:]/np.max(histograms[:,1:])                                  ##  makes histograms relative to maximum
    overlap_cutoff=np.mean(np.max(hist_rel))*0.1                                        ##  10% of total hisgram height discarded
    overlap=[]
    for i in range(len(histograms[0:,1:])):                                             ##  provides overlap count if there is a overlap above cutoff
        overlap.append(np.count_nonzero(hist_rel[i] > overlap_cutoff))
    return histograms[:,0], histograms[:,1:],hist_sum, overlap_cutoff, np.array(overlap), hist_rel

### fills in gaps from the histograms
def fill_gaps():
    print('\nfilling in gaps in PMF')
    folders()                                                                                       ##  makes_folders
    coord, histograms,histogram_sum, overlap_cutoff, overlap, histogram_rel=get_histograms()        ##  gets histogram data
    pull=get_pull()                                                                                 ##  gets pull information
    react_coord_proposed, react_coord_init=[],[]
    initial_offset, offset=args.offset, args.offset
    coord_cut = coord[np.where(np.logical_and(coord>np.min(pull[1]),coord<np.max(pull[1])))]
    coord_fill = np.where(np.logical_or(overlap < args.cutoff, histogram_sum <= np.mean(histogram_sum)*0.25))

    cv_start=True
    for i, cv in enumerate(coord_cut[coord_fill]):
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
        offset, react_coord_proposed, react_coord_init = get_conformation(cv_initial, cv_end, offset, pull,react_coord_proposed, react_coord_init)
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

    step_y= ask_number('PMF tick interval length on the Y axis [eg 10]:  ') ## y axis tick interval

    energy_time=[]              ## collects energy minima for timecourse
### plotting header
    rcParams['axes.linewidth']=3
    figure(1, figsize=(10,20))
    # #energy landscape
    subplot(4,1,1)
    title('PMF',  fontproperties=font1, fontsize=15,y=1.18)

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
            
            min_x, max_x, step_x=np.round(min(pmf_current[:,0]), 2)-0.25, np.round(max(pmf_current[:,0]), 2)+0.25, 1
            try:
                min_y, max_y = [float(x) for x in input("\nmin and max Y (press enter to use defaults) :  ").split()]
                print(min_y, max_y)
            except:
                print("\ndefault values used")
                min_y, max_y= np.round(min(pmf_current[:,1]), 2)-10, np.round(max(pmf_current[:,1]), 2)+10
                print(min_y, max_y)
        ### min max end ###
            energy_minima=np.round(min(pmf_current[:,1]), 0)
            energy_minima_error=0

            print('pmf '+str(pmf_number+1)+': energy minima = '+str(np.round(min(pmf_current[:,1]), 2))+' +/- '+ str(np.round( float(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)))
            plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, color='red', label='PMF 1')    
            if len(pmf_current[0]) == 3:
                fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')
                energy_minima_error=np.round( float(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)
            energy_time.append([energy_minima, energy_minima_error])
        else:
            energy_minima=np.round(min(pmf_current[:,1]), 0)
            energy_minima_error=0
            print('pmf '+str(pmf_number+1)+': energy minima = '+str(np.round(min(pmf_current[:,1]), 2))+' +/- '+ str(np.round( float(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)))
            plot(pmf_current[:,0],pmf_current[:,1], linewidth=3, label='PMF '+str(pmf_number+1))#, color='blue')
            if len(pmf_current[0]) == 3:
                fill_between(pmf_current[:,0], pmf_current[:,1]-pmf_current[:,2], pmf_current[:,1]+pmf_current[:,2], alpha=0.3, facecolor='black')  
                energy_minima_error=np.round( float(pmf_current[:,2][np.where(pmf_current[:,1] == min(pmf_current[:,1]))]), 2)
            energy_time.append([energy_minima, energy_minima_error])
    xticks(np.arange(-500,500,step_x), fontproperties=font1, fontsize=15);yticks(np.arange(-500,500,step_y), fontproperties=font1, fontsize=15)
    ylim(min_y, max_y);xlim(min_x,max_x)
    tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)
    xlabel('Distance (nm)', fontproperties=font2,fontsize=15);ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=15) 
    legend(prop={'size': 10}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(args.pmf), mode="expand", borderaxespad=0.)
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
    plot([-100,100],[args.cutoff,args.cutoff], linewidth=3, color='red')
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
    savefig('energy_landscape_'+timestamp+'.png', dpi=300)

    figure(2, figsize=(10,10))
    energy_time=np.array(energy_time)
    timestep=1
    if len(energy_time) > 1:
        timestep=ask_number('what is the timestep? ')
    # #energy landscape
    title('PMF timecourse',  fontproperties=font1, fontsize=15,y=1.18)
    # plot(np.arange(0+timestep,len(energy_time[:,0])*timestep+timestep,timestep),energy_time[:,0], linewidth=3, color='blue', zorder=1)
    errorbar(np.arange(0+timestep,len(energy_time[:,0])*timestep+timestep,timestep),energy_time[:,0], yerr=energy_time[:,1], color='red', zorder=1)
    plot(np.arange(0+timestep,len(energy_time[:,0])*timestep+timestep,timestep),energy_time[:,0],'o',color='k') 

    ylim(np.min(energy_time[:,0])-5, np.max(energy_time[:,0])+5);xlim(0,len(energy_time[:,0])*timestep+(2*timestep))
    # ylim(np.min(energy_time[:,0])-5,0);xlim(0,len(energy_time[:,0])*timestep+(2*timestep))

    tick_params(axis='both', which='major', width=3, length=5, labelsize=15, direction='in', pad=10, right=False, top=False)

    xlabel('time (ns)', fontproperties=font2,fontsize=15);ylabel('Energy (kJ mol$^{-1}$)', fontproperties=font2,fontsize=15) 

    savefig('energy_landscape_error_'+timestamp+'.png', dpi=300)

    show()

def results(miscs):
    react_coord_proposed, react_coord_init,react_coord_final=miscs
    result_write = open(location+'/setup_files_'+timestamp+'/collective_variable_position'+'_'+timestamp, 'w')
    if args.tpr==True:
        if len(react_coord_final)==len(react_coord_proposed):
            line='\n{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}'.format('proposed', 'selected', 'final', 'S-P', 'F-S', 'window') 
            print(line)
            result_write.write(line+'\n')
        else:
            sys.exit('mismatch between tpr output and expected output, check gromacs_outputs')
    elif args.tpr==False:
        line='{0:^10}{1:^10}{2:^10}{3:^10}'.format('proposed','selected','S-P','window')   
        print(line)
        result_write.write(line+'\n')
    for i in range(len(react_coord_proposed)):
        if args.tpr==True:
            try:
                line='{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}'.format(react_coord_proposed[i], react_coord_init[i], react_coord_final[i], 
                            np.around(float(react_coord_init[i])-float(react_coord_proposed[i]), decimals=3), 
                            np.around(float(react_coord_final[i])-float(react_coord_init[i]), decimals=3 ), args.offset+1+i)
                print(line)
            except ValueError:
                line=str(args.offset+1+i)+'   error in final grompp most likely'
            result_write.write(line+'\n')
        elif args.tpr==False:
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
        files_range, time, pull = [], [], []
        if args.current:
            xvg_loc = './'
        else:
            xvg_loc = '../windows/window_'+str(window)+'/'
        
        for root, dirs, files in os.walk(xvg_loc):
            for filename in files:
                if 'window_'+str(window)+'.' in filename or 'window_'+str(window)+'_' in filename and not '_com' in filename:
                    if 'pullf' in filename:
                        files_range.append(filename)
            break

        if len(files_range) == 1:
            xvgs=len(files_range)
            os.system('cp '+xvg_loc+files_range[0]+' window_'+str(window)+'_pullf_com.xvg')     
        elif len(files_range) >= 2: 
            xvgs=len(files_range)
            for x in range(len(files_range)):
                for line in open(xvg_loc+files_range[x], 'r').readlines():
                    if not line[0] in ['#', '@']:
                        if len(line) > 0 and len(line.split()) ==2:
                            try:
                                time.append(float(line.split()[0]))
                                pull.append(float(line.split()[1]))
                            except:
                                break
            if len(time) > 0:
                time_ord, pull_ord = (list(t) for t in zip(*sorted(zip(time, pull))))
                em = open('window_'+str(window)+'_pullf_com.xvg','w')
                for j in range(len(time_ord)):
                    em.write(str(time_ord[j])+'\t'+str(pull_ord[j])+'\n')
        else:
            xvgs=len(files_range)
        return [window, xvgs]

def run_wham():
    core = str(ask_integer('what core would you like to run this on 0-11: '))
    gromacs(['taskset --cpu-list '+core+' gmx wham -if en.dat -it tpr.dat -bsres '+args.pmf+' -temp 310 -nBootstrap '+str(args.boot)+' -b '+str(args.start)])

start = time.time()

parser = argparse.ArgumentParser()

parser.add_argument('-mdp', help='umbrella mdp file',metavar='md.mdp', type=str)
parser.add_argument('-func', help='what to do initial setup, plot, concat, wham or fill',metavar='[setup, plot, concat, wham, fill]',type=str, choices= ['setup', 'plot', 'concat','wham', 'fill'])
parser.add_argument('-f', help='xtc file for setup (eg your pull xtc)',metavar='pull.xtc',type=str)
parser.add_argument('-s', help='structure file for setup (use the pull.tpr)',metavar='pull.tpr',type=str)
parser.add_argument('-n', help='index file for the system',metavar='index.ndx', type=str)
parser.add_argument('-p', help='topology file',metavar='topol.top', type=str)
parser.add_argument('-pull', help='pull file for setup',metavar='pullx.xvg',type=str)
parser.add_argument('-offset', help='window offset',metavar='5',type=int, default=0)
parser.add_argument('-tpr', help='do not make tpr files', action='store_false')
parser.add_argument('-min', help='switch off minisation', action='store_false')
parser.add_argument('-int', help='interval for umbrella windows (nm)',metavar='0.05', type=float, default=0.05)
parser.add_argument('-start', help='where to start on reaction coordinate',metavar='0',type=float)
parser.add_argument('-end', help='where to end on reaction coordinate',metavar='5', type=float)
parser.add_argument('-boot', help='number of bootstraps to run',metavar='5', type=int)
parser.add_argument('-pmf', help='location of pmf ',metavar='bsres.xvg',type=str, nargs='*')
parser.add_argument('-hist', help='location of histogram and name if used with wham',metavar='histo.xvg',type=str)
parser.add_argument('-tpronly', help='only makes tpr files default (False) requires energy minimised files', action='store_true')
parser.add_argument('-current', help='to concat in current directory', action='store_true')
parser.add_argument('-cutoff', help='histogram overlap cutoff for filling and plotting', default=3, type=int)
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels) (Optional)")


args = parser.parse_args()
options = vars(args)

if args.func == 'setup':
    correct = check_arguments(['pull', 'int'])
    if args.tpr:
        correct = check_arguments(['mdp', 'tpr','f', 'n', 'p'])
    if correct:
        misc = setup()
        results(misc)
elif args.func == 'fill':
    correct = check_arguments(['pull', 'int', 'offset'])
    if args.tpr:
        correct = check_arguments(['mdp', 'tpr','f', 'n', 'p'])
    if correct:
        misc = fill_gaps()
        results(misc)
elif args.func== 'plot':
    correct = check_arguments(['pmf', 'hist'])
    if correct:
        plot_pmf()
elif args.func== 'concat':
    pool = mp.Pool(mp.cpu_count())
    correct = check_arguments(['start', 'end'])
    if correct:
        concatonated = pool.map(pull_concat, [(window) for window in range(int(args.start), int(args.end)+1)])          ## makes umbrella windows from minimised frames
        pool.join
        np.array(concatonated).sort(axis=0)
        print('\nwindow\ttotal part numbers')
        for line in concatonated:
            if line[1] != 0:
                print(line[0], '\t', line[1])
            else:
                print(line[0], '\t', 'SKIPPED')
elif args.func== 'wham':
    correct = check_arguments(['pmf', 'boot', 'start'])
    if correct:
        run_wham()
else:
    sys.exit('Try again, enter  [-f setup, plot, fill or concat]')


end = time.time()

print('\nThis script took: '+str(np.round(end-start,1))+'s to run')