
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

while ( $s <= 100 )
  foreach dist ( unif ) 
    foreach n ( 100 500 128 256 512 1024 )
		  # par1: $n sample size, 100 -> 1000
		  # par2: Dist of x, unif, exp, norm, grid, mix
		  # par3: $s seed value

		  qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulationsAdd.csh $n $dist $s
     end
  end
  @ s++
end




