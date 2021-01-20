#
#$1 input is site
/bin/rm rotate.macro list_of_files
ls *.Z.$1.SAC | awk '{print substr($1,1,12)}' > list_of_files
foreach i (`cat list_of_files`)
echo "r " $i".E.$1.SAC" $i".N.$1.SAC" >> rotate.macro.$1
echo "rmean;rtrend;rmean" >> rotate.macro.$1
echo "highpass corner 0.02 npoles 4" >> rotate.macro.$1
echo "rot to gcp reversed" >> rotate.macro.$1
echo "w " $i".$1.r" $i".$1.t" >> rotate.macro.$1
echo "r " $i".Z.$1.SAC" >> rotate.macro.$1
echo "rmean;rtrend;rmean" >> rotate.macro.$1
echo "highpass corner 0.02 npoles 4" >> rotate.macro.$1
echo "w " $i".$1.z" >> rotate.macro.$1
echo "  " >> rotate.macro.$1
end

#produce a list with station names in it
awk '{print $1 "." A}' A=$1 list_of_files > pick-list.$1

echo "run rotate.macro.$1 then execute make-pick-macro pick-list.$1"