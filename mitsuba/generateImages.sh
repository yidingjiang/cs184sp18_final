#!/bin/bash
for number in {0..0}
do
C:\\Users\\cs199-btx\\Desktop\\Mitsuba0.5.0\\mitsuba -DframeNum=$number -o water_$number.exr water.xml
done
exit 0
