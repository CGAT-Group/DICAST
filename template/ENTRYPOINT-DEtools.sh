#!/bin/bash
source config.sh
source func.sh

#see func.sh for function help
readsamfiles
for filename in $wd/tmp/$tool-samlist do;
makebamfromsam $filename
done;

#run tool
# example: CASH java -jar -Xmx10g cash.jar \ --Case:Case $bamCase \ --Control:Control $bamControl \ --GTF $wd/$gtf --Output $wd/$out 
#wait for all processes to end
 wait 
#make output accessible and clean up temp files
trap cleaner EXIT
