#!/bin/bash
/Users/darren/Development/sarcastic/scripts/fixHeaders.py $1 > /tmp/fhtmp 
if [ $? -eq 0 ]
then
	mv /tmp/fhtmp $1
else
	cat /tmp/fhtmp
fi
