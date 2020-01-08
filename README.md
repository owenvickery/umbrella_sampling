UMBRELLA SAMPLING README

If you are using this script please acknowledge me (Dr Owen Vickery) and cite the following DOI please.

DOI: 10.5281/zenodo.3592318

---------------------------------------------------------------------------------------

                                        REQUIREMENTS

Python v3 or higher

Numpy

---------------------------------------------------------------------------------------

                                        FLAGS

This script sets up and analyses umbrella sampling trajectories.

Flags to use

      -h, --help         
      -mdp (.mdp)        
      -func (setup, concat, wham, plot, fill)             
      -f (.xtc)         
      -s (.tpr)          
      -n (.ndx)          
      -p (.top)          
      -pull (.xvg)       
      -offset (int)      
      -tpr              
      -min               
      -int (float)      
      -start (float)       
      -end (float)   
      -boot (int)     
      -pmf (list)     
      -hist (.xvg)      
      -tpronly    
      -current     
      -cutoff (int) 
      -v


This script is designed to follow on from your initial pull simulation.

---------------------------------------------------------------------------------------

                                        SETUP

To setup a new PMF use the flag (-func setup). As this sets up all the umbrella windows the script requires all the files for running the grompp step.

eg.
index, 
topology,
production mdp, 
structure, 
pull tracjectory,
pull reaction coordinate,
window interval, 
pull start,
pull end.

for example if you have all the build files in the folder 'BUILD' and all the pull files in the folder 'PULL' and wish to run the umbrella sampling between 1 and 3 nm with a window spacing of 0.05 nm.

The pull reaction coordinate file has to be 2 columns time and CV. 

the command line would look like:

    python pmf.py -func setup -n BUILD/index.ndx -p BUILD/topol.top -mdp BUILD/production.mdp -s PULL/pull.tpr -f PULL/pull.xtc -pull PULL/pullx.xvg -start 1 -end 3 -int 0.05


If you wish to extend the PMF from 3 to 4, you can use the offset flag (e.g. if you have 40 windows already, enter 40. the script will write windows from 41 onwards (flag -offset 40)).

e.g.

    python pmf.py -func setup -n BUILD/index.ndx -p BUILD/topol.top -mdp BUILD/production.mdp -s PULL/pull.tpr -f PULL/pull.xtc -pull PULL/pullx.xvg -start 1 -end 3 -int 0.05 -offset 40

The scipt will provide you with a overview of your PMF windows.

Column 1 (proposed) is the CV windows you choose.

Column 2 (selected) is the CV windows selected by the script ( closest value to Col 1 )

Column 3 (final) is the CV in the final production window

Column 4 (S-P) is the selected CV minus the proposed CV

Column 5 (F-S) is the final CV minus the selected CV

Columm 6 is the window number for each CV used in the script

    proposed       selected         final           S-P              F-S         window
      0.5           0.5             0.510           0.0             0.01            40
      0.575         0.575           0.573           0.0             -0.002          41
      0.65          0.65            0.662           0.0             0.012           42
      0.725         0.726           0.736           0.001           0.01            43
      0.8           0.799           0.806           -0.001          0.007           44

The output from this script is in the following format.

    | --    umbrella_windows
                | --    setup_files_(timestamp)
                                - collective variables, gromacs outputs, topology, mdp, index, em_out-[1..40]
                | --    frames
                                - window_[1..40].pdb
                | --    minimised
                            | --   window_[1..40]
                                        | --    tpr, gro, edr, log
                | --    windows
                            | --   window_[1..40]
                                        | --    production tpr
                | --    analysis

There are a couple of other flags for for more odd cases.

If you do not wish to make the tpr files (flag -tpr), however you can make the tpr files only but this requires energy minimised structures (flag -tpronly).
if you wish to skip energy minimisation (flag -min).

---------------------------------------------------------------------------------------

                                        WHAM ANALYSIS

Once you have finished your initial PMF, you need to analyse it. Here the script can run a rudimentary analysis on your PMF. 

1st you need to fetch all the pullf files from the windows directory. The script can do this for you.

To concatonate part files if you do your umbrella sampling with the noappend flag, it will read the windows directory and copy them to the analysis folder.

If they are already in your current directory you can use the flag -current.

use the flags -start/-end to specify window range to concatonate.

If you change to the analysis directory and run the following command.



e.g. from windows directory

    python pmf.py -func concat -start 1 -end 40 

e.g. in current directory

    python pmf.py -func concat -start 1 -end 40 -current



Your analysis folder should now contain the pullf.xvg and .tpr files from the windows directory.

To analysis your create your landscape with gmx wham you need to files (please use these names): tpr.dat and en.dat 

These files contain a single column of either the names of the tpr files of the pullf files. Note that they have to be in the same order.

e.g.

      tpr.dat               en.dat
    window_1.tpr    window_1_pullf_com.xvg
    window_2.tpr    window_2_pullf_com.xvg
    window_3.tpr    window_3_pullf_com.xvg
    window_4.tpr    window_4_pullf_com.xvg
    window_5.tpr    window_5_pullf_com.xvg


The script can run wham for you, however it only uses the basic setting (This wham is run at 310K) so I would advise you to run it separately.

The flag -pmf provides the output file name of your PMF.
The flag -boot is the number of bootstraps to do.
The flag -start is the amount of time to dicard as equilibration.

e.g. 

    python pmf.py -func wham -pmf bsres.xvg -boot 200 -start 5000

For some reason gmx wham runs equally well on 1 core as it does on all the cores. This allows you to run multiple whams simultaneously.
Therefore the script will ask you which core to run on (note they cores start from 0)  .


---------------------------------------------------------------------------------------

                                        PLOTTING

Once you have your intial PMF and you need to run some basic quality control checks.

Here this script will provide:

Energy landscape and error
Normalised histogram sum
Histogram overlap
Histograms

The red lines denote my level of quality I want in the PMF. 
The normalised histogram sum has a cutoff of 20 %.
The Histgram overlap has a default cutoff of 3 (changed with flag -cutoff). 

The flag -pmf provides the name of your energy landscape.
The flag -hist provides the name of your histograms.
The flag -cutoff allows the change of the overlap cutoff.

e.g. 

    python pmf.py -func plot -pmf bsres.xvg -hist histo.xvg

you be asked various questions, reply with a numerical value or return.

    PMF tick interval length on the Y axis [eg 10]: 10
    min and max Y (press enter to use defaults) : 0 100      (min and max separated by a space)

The plot will be saved as energy_landscape_(time_stamp).png
Also the energy minima is also plotted as a line graph.


If you wish to check the convergence of the pmf with increasing simulation time.

e.g. 
    5-10 ns 
    5-15 ns
    5-20 ns

you can provide multiple bsres.xvg files to the -pmf flag. (each file should be separated by a space).

e.g. 

    python pmf.py -func plot -pmf bsres_5-10.xvg bsres_5-15.xvg bsres_5-20.xvg -hist histo.xvg

In this case you be asked a additional question. 

    what is the timestep? 5

This sets the x axis of the line graph showing the energy minima over time.

---------------------------------------------------------------------------------------

                                        GAPS

Unless you have a simple system, you will get gaps between your umbrella windows.

Therefore this script will fill in the gaps using the same criteria from the analysis plot.

The normalised histogram sum has a cutoff of 20 %.
The Histgram overlap has a default cutoff of 3 (changed with flag -cutoff). 

To fill in the gaps in the umbrella windows coverage we will use a similar command as in the setup.  

However you need to provide the histogram file as well and the offset flag otherwise it will overwrite your existing windows.

e.g.

    python pmf.py -func fill -n BUILD/index.ndx -p BUILD/topol.top -mdp BUILD/production.mdp -s PULL/pull.tpr -f PULL/pull.xtc -pull PULL/pullx.xvg -start 1 -end 3 -int 0.05 -offset 40 -hist umbrella_sampling/analysis/histo.xvg

Iterate over these steps until your have a coverged PMF.

Good luck    