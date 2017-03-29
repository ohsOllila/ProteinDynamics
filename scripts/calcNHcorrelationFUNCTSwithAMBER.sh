ndxFILE=/home/samuli/Data/Calmodulin/02-CaM_Amber99sb-ildn/02-CaM-4MK_noscaled/03-md_noPOSRES_1000ns_good_itp/HN.ndx
trajectory=/home/samuli/Data/Calmodulin/02-CaM_Amber99sb-ildn/02-CaM-4MK_noscaled/03-md_noPOSRES_1000ns_good_itp/02-total_center.xtc
trajectoryORIENTED=/home/samuli/Data/Calmodulin/02-CaM_Amber99sb-ildn/02-CaM-4MK_noscaled/03-md_noPOSRES_1000ns_good_itp/02-total_centerFIT.xtc
tpr=/home/samuli/Data/Calmodulin/02-CaM_Amber99sb-ildn/02-CaM-4MK_noscaled/03-md_noPOSRES_1000ns_good_itp/01-1CLL_md_noPOSRES_1000ns.tpr
gro=/home/samuli/Data/Calmodulin/02-CaM_Amber99sb-ildn/02-CaM-4MK_noscaled/03-md_noPOSRES_1000ns_good_itp/03-total_center_CaM.gro
correlationFUNCTfolder=~/Dropbox/Calmodulin/Data/NHcorrelationFunctions/withAMBER/unoriented/
correlationFUNCTfolderORIENTED=~/Dropbox/Calmodulin/Data/NHcorrelationFunctions/withAMBER/oriented/
awk -f ~/Dropbox/CBD64/scripts/makeNHindex.awk $gro > $ndxFILE
numberOFfuncs=$(grep "\[" $ndxFILE | tail -n 1 | awk '{print $2}')
for((i=1;i<=$numberOFfuncs;i++))
do
    echo $i | gmx rotacf -f "$trajectory" -s $tpr -n $ndxFILE -o $correlationFUNCTfolder/NHrotaCF_$i.xvg -P 2 -d -xvg none  #-nice 20 &
    echo $i | gmx rotacf -f "$trajectoryORIENTED" -s $tpr -n $ndxFILE -o $correlationFUNCTfolderORIENTED/NHrotaCF_$i.xvg -P 2 -d -xvg none  #-nice 20 &
done
