#!/bin/bash
for number in {44..213}
do
/Applications/Mitsuba.app/Contents/MacOS/mitsuba -DframeNum=$number -o water_$number.exr water.xml
done
exit 0