#!/bin/bash
for file in `ls -1 *.minimap2.tab`; do
	while read line; do
		middle=`echo -e "$line" | cut -f 2-5`
		end=`echo -e "$line" | cut -f 7-`
		sstart=`echo -e "$line" | cut -f 3`
		send=`echo -e "$line" | cut -f 4`
		let "slength = $send - $sstart"
		slength=${slength#-}
		if [[ $slength -gt 10000 ]]; then
			echo -e "$line" >> ${file}.10kb.tab
		fi
	done < $file
done
