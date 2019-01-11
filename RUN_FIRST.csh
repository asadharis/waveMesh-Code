
#!/bin/csh

#Run a loop to do our analysis for each seed

foreach n ( 64 128 256 512 )
		# par1: $n sample size: 64, 123, ...,  514
		# par2: Dist of x: unif, norm
		# par3: Regression fucntion: f.poly, f.sin, ..., f.bumps

		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.poly
		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.sin
		
		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.heavi
		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.ppoly

		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.doppler
		qsub  -cwd -q "shojaie*"  -o /dev/null -e /dev/null  simulations.csh $n unif f.bumps

end



