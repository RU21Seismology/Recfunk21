# after going through the data with SAC make-pick-macro
# this shell will look through the data files, and collect those that 
# have flags T1 or T2 set into a list for use with recfunk09_pick
# !!!to run this you need to have a code called sachdr on your machine!!!
#Input - same list of events as for make-sac-macro
# Output - a file to be used as input for recfunk09_input,
# and a list of events where both T1 and T2 got marked (erroneously) - if it has
# any entries, examine the files, and unset one of the flags

/bin/rm $1.recpick_input $1.two_marks
set DIR = `pwd`

set t1 = " "
set t2 = " "
set t3 = " "

foreach i ( `cat $1` )
set t1 = `sachdr -l $i.r | grep t1 `
set t2 = `sachdr -l $i.r | grep t2 `
set t3 = `sachdr -l $i.r | grep t3 `
echo $DIR $i $t1 $t2 $t3 | awk ' NF > 2  {print $1 "/" $2 ".?"}' >> $1.recpick_input
echo $i $t1 $t2 $t3 | awk ' NF > 4' >> $1.two_marks
end
