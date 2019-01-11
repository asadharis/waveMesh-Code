
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

while ( $s <= 100 )
	qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  dataAdd.csh $s
  @ s++
end




