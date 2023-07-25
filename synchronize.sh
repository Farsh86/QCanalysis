#!/bin/bash
# Read Password
echo -n "Password for itsshift:" 
read -s password

while true
do
	echo -n "$password" | kinit itsshift
	rsync -avruh Plots/ShiftReport24h_THR_* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/THR
	rsync -avruh Plots/ShiftReport24h_FHR_* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/FHR
	rsync -avruh Data/Output_all-IB-layers_FHRMAPS* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/Output/FHR
        rsync -avruh Data/Output_all-IB-layers_THRMAPS* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/Output/THR
	rsync -avruh logs/logTHR* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/logs/THR
	rsync -avruh logs/logFHR* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/logs/FHR
	echo "Sleeping 1h"
	sleep 1h
done

