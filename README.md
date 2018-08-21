# umbrella_sampling
setting up and analysing PMFs 

To setup a new pmf or fill a specific place (flag -f setup):
files needed:
specify location of the pull directory (flag -pull)
start of collective variable you wish to sample (flag -s)
end of collective variable you wish to sample (flag -e)
location where you want the configuration windows eg pdb input files (flag -window)
interval between umbrella windows in nm eg 0.05 spaces a window every 0.5 Angstroms (flag -int)

you will be asked for a direction to setup +/- or just press enter to not incorporate a direction.
you will be asked for the window offset, if you have 100 windows already enter 100. the script will write windows from 101 onwards.
you will be asked whether you want to use the pullx file from the pull directory or a user defined one.
You will finally be asked if you want to create the tpr files for the configuration files.

If you wish to make the tpr files in addition (need to be specified at start of script):
the mdp file needs to be specified (flag -mdp)
the topology file needs to be specified (flag -p)
the index file needs to be specified (flag -n)
the location where you want the windows to be outputed (flag -o) 

you get 2 outputs:
the gromacs outputs if you wish to check for errors or prosterity.
the collective variable: proposed, found from xvg file, the final coordinate from the grompp and the corresponding window

To fill in missing gaps (flag -f fill):
files needed:
specify location of the pull directory (flag -pull)
location where you want the configuration windows eg pdb input files (flag -window)
interval between umbrella windows in nm eg 0.05 spaces a window every 0.5 Angstroms (flag -int)

you will be asked for a direction to setup +/- or just press enter to not incorporate a direction.
you will be asked for the window offset, if you have 100 windows already enter 100. the script will write windows from 101 onwards.
you will be asked whether you want to use the pullx file from the pull directory or a user defined one.
You will finally be asked if you want to create the tpr files for the configuration files.

If you wish to make the tpr files in addition (need to be specified at start of script):
the mdp file needs to be specified (flag -mdp)
the topology file needs to be specified (flag -p)
the index file needs to be specified (flag -n)
the location where you want the windows to be outputed (flag -o) 

To plot (flag -f plot):
files needed:
The energy landscape with error, ie bsres.xvg from gmx wham (flag -pmf)
The histogram file from gmx wham (flag -hist)

To run wham (flag -f wham):
files needed:
location of the tpr.dat file contains list of tprs (flag -tpr)
location of the en.dat file contains pullx file list (flag -en)
The energy landscape output directory  (flag -pmf)
The histogram file  output directory (flag -hist)
The profile landscape output directory  (flag -profile)
The bootstrap profile file  output directory (flag -bsprof)
the time to skip from start of pmf (flag -time)
the temperature of the system (flag -temp)
the number of bootstraps (flag -boot)
any extra commands should be inputed as 'min 5' 'max 6' eg "-extra '-min 5 -max 6'"  

You will be asked if you want to plot



