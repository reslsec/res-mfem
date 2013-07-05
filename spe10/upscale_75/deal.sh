#!/bin/bash
grep Producer1 "$1">p1.txt
sed -i 's/,/ /g' p1.txt
sed -i 's/Producer1:/ /g' p1.txt
grep Producer2 "$1">p2.txt
sed -i 's/Producer2:/ /g' p2.txt
sed -i 's/,/ /g' p2.txt
grep Producer3 "$1">p3.txt
sed -i 's/Producer3:/ /g' p3.txt
sed -i 's/,/ /g' p3.txt
grep Producer4 "$1">p4.txt
sed -i 's/Producer4:/ /g' p4.txt
sed -i 's/,/ /g' p4.txt
grep average_press "$1">avp.txt
sed -i 's/Time:/ /g' avp.txt
sed -i 's/average_press:/ /g' avp.txt
sed -i 's/day,/ /g' avp.txt
sed -i 's/Psi/ /g' avp.txt
grep inje_press "$1">inje.txt
sed -i 's/Time:/ /g' inje.txt
sed -i 's/inje_press:/ /g' inje.txt
sed -i 's/day,/ /g' inje.txt
sed -i 's/Psi/ /g' inje.txt
