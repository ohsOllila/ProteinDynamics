## CALCULATE INERTIA TENSORS AND ROOT MEAN SQUARE ANGLE DEVIATIONS

# Source for Python script to calculate inertia tensors and root mean square angle deviations 
calcINERTIA=~/Dropbox/ProteinDynamics/scripts/calcINERTIA2.py

# Source for trajectory and *gro file with protein centered
trajPROTcent=/Users/osollila/Calmodulin/withAMBERchargesSCALED/03-total_center_CaM_400-1000ns.xtc
groPROTcent=/Users/osollila/Calmodulin/withAMBERchargesSCALED/03-total_center_CaM_400-1000ns.gro

#Output files for calcINERTIA.
INERTIAoutput=../Data/amber400-1000nsSCALED

python3  $calcINERTIA $trajPROTcent $gorPROTcent $INERTIAoutput


## CALCULATE OVERALL CORRELATION FUNCTIONS

# Source for correlation function folders
# This folder should contain subfolders oriented and unoriented with calculated correlation functions
# and empty subfolders overall and scaledrotation.
correlationFUNCfolder=../Data/NHcorrelationFunctions/withCHARMMscaled/

for((i=0;i<=144;i++))
do
    paste $correlationFUNCfolder/oriented/NHrotaCF_"$i".xvg $correlationFUNCfolder/unoriented/NHrotaCF_"$i".xvg | awk '{if($0!~"#" && $0!~"@" && $0!~"&")print $1" "$4/$2}' > $correlationFUNCfolder/overall/NHrotaCF_"$i".xvg
done

## CALCULATE OVERALL ROTATION DIFFUSION CONSTANTS AND SCALED CORRELATION FUNCTIONS

# Source for Python script to calculate rotational diffusion constants and
# correlation functions with scaled diffusion constant
calcROTdiff=~/Dropbox/ProteinDynamics/scripts/calcROTdiff2.py

# File to save the rotational diffusion constants
ROTdiffs=../Data/rotationalDIFFUSION.dat

python3 $calcROTdiff $INERTIAoutput"RMASD.dat" $correlationFUNCfolder > $ROTdiffs
