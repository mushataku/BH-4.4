
#!/usr/bin/bash

while IFS= read row; do
  No=`echo ${row} | cut -d , -f 1`
  r_ring=`echo ${row} | cut -d , -f 2`
  beta=`echo ${row} | cut -d , -f 3`
  spin_a=`echo ${row} | cut -d , -f 4`
  
  if [ "${No}" != "Run Number" ]; then
    ./a.out $r_ring $beta $spin_a
  fi
done  < calculation_condition.csv

