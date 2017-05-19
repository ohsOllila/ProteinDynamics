## CALCULATE INERTIA TENSORS AND ROOT MEAN SQUARE ANGLE DEVIATIONS

# Source for Python script to calculate inertia tensors and root mean square angle deviations 
calcINERTIA=/wrk/einstein/samuli/ProteinDynamics/scripts/calcINERTIA2.py

# Source for trajectory and *gro file with protein centered
trajPROTcent=/wrk/einstein/pauline/Calmodulin-felicitas-BK/02-CaM_Amber99sb-ildn/09-CaM_4MK_scaled_resid/03-total_center_CaM.xtc
groPROTcent=/wrk/einstein/pauline/Calmodulin-felicitas-BK/02-CaM_Amber99sb-ildn/09-CaM_4MK_scaled_resid/03-total_center_CaM.gro

#Output files for calcINERTIA.
INERTIAoutput=./

python3  $calcINERTIA $trajPROTcent $groPROTcent $INERTIAoutput


## CALCULATE OVERALL CORRELATION FUNCTIONS

# Source for correlation function folders
# This folder should contain subfolders oriented and unoriented with calculated correlation functions
# and empty subfolders overall, overallFIT and scaledrotation.
correlationFUNCfolder=./

# Amount of correlation functions in subfolders
NumberOfCorrFs=144

for((i=0;i<=$NumberOfCorrFs;i++))
do
    paste $correlationFUNCfolder/oriented/tmp/NHrotaCF_"$i".xvg $correlationFUNCfolder/unoriented/tmp/NHrotaCF_"$i".xvg | awk '{if($0!~"#" && $0!~"@" && $0!~"&")print $1" "$4/$2}' > $correlationFUNCfolder/overall/NHrotaCF_"$i".xvg
done

## CALCULATE OVERALL ROTATION DIFFUSION CONSTANTS AND SCALED CORRELATION FUNCTIONS

# Source for Python script to calculate rotational diffusion constants and
# correlation functions with scaled diffusion constant
calcROTdiff=/wrk/einstein/samuli/ProteinDynamics/scripts/calcROTdiff2.py

# File to save the rotational diffusion constants
ROTdiffs=./rotationalDIFFUSIONAmber99SCALED.dat

# Scaling factor for new rotational diffusion constants
scalingF=2.9

# Lag time (ns) used for linear fit for diffusion coefficient.
# I usually use one hundredth of total length of the analyzed trajectory
LagTimeForDiffusionConstant=10

python3 $calcROTdiff $INERTIAoutput"RMASD.dat" $correlationFUNCfolder $scalingF $LagTimeForDiffusionConstant $NumberOfCorrFs > $ROTdiffs
