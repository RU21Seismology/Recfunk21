# shell to prepare data 
# for RF analysis with Park-Levin MTC codes, version 2021
# assumes data to be in SAC format
# the shell produces two files - a list of events and 
# Input - $1 station name
# YEAR.JDAY.HH.MM.SITE.C, where C stands for z, r or t, an entry in the list will be YEAR.JDAY.HH.MM.SITE
#
# Output a SAC macro that will open each event and let you put time marks in
#
# this shell assumes that data rotation and necessary 
# cleanup, like demeaning and detrending and maybe a high-pass filter,
# have all been done)
#

/bin/rm $1.recpick_input *.sac.macro

#below it is assumed that data were cut with 300 sec pre-event window.
# if the window is different - change the value 
# TAKES STATION NAME AS INPUT

foreach i (`ls *.z | awk '{print substr($1,1,12)}'`)
echo "r $i.$1.?; chnhdr a 300; listhdr gcarc baz;  plotpk markall on; writehdr over" >> $1.sac.macro
end

echo RUN SAC MACRO
	

