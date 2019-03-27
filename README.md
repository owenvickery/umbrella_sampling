# umbrella_sampling
setting up and analysing PMFs 

To setup a new pmf or fill a specific place (flag -func setup):

specify the pullx file (has to be 2 columns; time and CV) (flag -pull)

specify a structure file, eg your pull tpr (flag -s)

start of collective variable you wish to sample (flag -start)

end of collective variable you wish to sample (flag -end)

interval between umbrella windows in nm eg 0.05 spaces a window every 0.5 Angstroms (flag -int)

window offset, if you have 100 windows already, enter 100. the script will write windows from 101 onwards (flag -offset

If you do notwish to make the tpr files (flag -tpr):
the mdp file needs to be specified (flag -mdp)
the topology file needs to be specified (flag -p)
the index file needs to be specified (flag -n)

if you wish to skip energy minimisation use -min flag


To fill in missing gaps (flag -func fill) use the above commands and it will fill in any gaps (remember to use -offset flag otherwise it will overwrite you existing windows)

To plot (flag -f plot):
files needed:
The energy landscape with error, ie bsres.xvg from gmx wham (flag -pmf)
The histogram file from gmx wham (flag -hist)

you be asked various questions, reply with Y/N or a numerical value.

To concatonate part files if you do your umbrella sampling with the noappend flag, it will read the windows directory.
use the flags -start/-end to specify window range to concatonate.

to run wham (I would do manually) but you can use this code (flag -wham).
-pmf specifies pmf output name
-boot specifies number of bootstraps
-start specifies when to run the wham from in picoseconds (eg for equilibration)

2 files are hardcoded, you need a file with your tprs and one with your pullf.xvg files listed. these are called tpr.dat and en.dat 
This wham is run at 310K 




