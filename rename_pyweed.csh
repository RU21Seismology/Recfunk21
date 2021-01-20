# this renames the name of SAC files downloaded using pyweed from whatever junk to YYYYMMDDHHMM.C.SAC. There are many events so if there are events occuring on the same day and causing complications just abandon them. This is for both BH? and HH? channels. Also, pyweed sometimes miss a channel for one event, for weird reasons. 
# also takes the station name as the input

/bin/rm do-rename

foreach i (E N Z)

ls *.HH$i* | awk '{print "mv " $1, substr($1,14,4) substr($1,19,2) substr($1,22,2) substr($1,25,2) substr($1,28,2) "." A "." B ".SAC"}' A=$i B=$1 >> do-rename

#ls *.BH$i* | awk '{print "mv " $1, substr($1,14,4) substr($1,19,2) substr($1,22,2) substr($1,25,2) substr($1,28,2) "." A "." B ".SAC"}' A=$i B=$1 >> do-rename

end

