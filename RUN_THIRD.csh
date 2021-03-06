
#!/bin/csh

#Run a loop to do our analysis for each seed

foreach f ( f.poly f.sin f.ppoly f.heavi f.doppler f.bumps )
  foreach dist ( unif norm exp mix ) 
    foreach n ( 75 100 300 500 )
		  # par1: $n sample size, 100 -> 1000
		  # par2: Dist of x, unif, exp, norm, grid, mix
		  # par3: Regression fucntion, f.poly, f.exp, f.sin

		  qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulation3.csh $n $dist $f 
     end
  end
end




