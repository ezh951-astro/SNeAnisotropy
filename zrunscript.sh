#!/bin/bash -l
#SBATCH -J PyAnalys
#SBATCH -p saleslab
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --time=6:00:00
#SBATCH -o zPyOutLog/pysc.out
#SBATCH -e zPyOutLog/pysc.err
#SBATCH --mail-user=ezhan039@ucr.edu
#SBATCH --mail-type=ALL

#cd /rhome/ezhan039/bigdata/analysis/FLAGSHIP-DIRWGT-M19
cd /rhome/ezhan039/bigdata/analysis/zAnisotropyPaper24

python F6Gif.py



#python CalculateOutflowData.py
#python HsmlDensity.py
#python morphology.py -r 6789 -v 0221
#python transverse_velocity.py -r 01234567 -v 1127
#python zGifMaker.py
#python velocity_projection.py -r 01234567 -v 1127
#python newstars.py -r 012345 -v 1208
#python rvphase.py -r 13 -v 1103
#python CustomDensity.py -r 0123 -v 1107
#python zGifMaker.py
