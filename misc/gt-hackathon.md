# GT HPC Hackathon Prep!

## Dates

* Mentor meeting: Jan 11 (Tues)
* Day 1: Jan 18 (Tues)
* Day 2-4: Jan 24-26 (Mon-Wed)

## Things to have done before we begin (Jan 18)

**Focus on Summit!**

* Mauro: 
	* Figuring out what's going on with big jobs, parallel **pre-processing**, running out of RAM on each node
		* Hypotheses: (1) Just had a parallel problem, (2) Running out of integers (this also doesn't matter for 2D/3D)
			* Seems like for (1) the memory problem is the main thing, 
			* Fix: Not pushing all domain data to each core? (I think)
			* Can we fix this so we can do bigger scaling runs! (> 768 GPUs)
	* Merging phase-change with hypoelasticity branches

* Anand
	* Implement Order-3 RK Time stepper
	- Done (RK2, RK3)
	* Parallel IO (with Nsight Systems), say every 100 time steps or so
	- Done (Parallel IO takes less than 2 timesteps)
	* Do some preliminary scaling runs on Summit (starting with weak scaling)
		* Do a 1D and a 3D case, two alphas/components if possible
	* Note: Can't merge GPU->master branch without time-steppers and multi-components verified working
	- Done (Verified timesteppers as well as multi-component)

* Jean
	* Post to Slack the results of scaling on XSEDE machines?
	* Compare CPU performance of hypoelastic vs regular cases (around 2x slowdown expected?)
	* Merge GPU-3D-unmanaged -> Hypoelastic; Hypoelastic -> GPU (merge)
	* Mention to Jose that he can join and 'submerge' himself into how this works (no key thrusts for him, though)

* Henry
	* Build, test, verify that new MFC master branch build system works on Summit
	* Getting new build system working with GPUs
	* Watch OpenACC YouTube video series, might be an Nvidia (separate) series on YT as well?

* Esteban
	1. Get Ascent access
	2. Build MFC GPU-3D-unmanaged branch on Ascent/Summit
	3. Make sure it works and you know how to run it properly
	4. Watch OpenACC YouTube video series, might be an Nvidia (separate) series on YT as well?
	4. Start adding the simplest possible chemistry kernel into MFC


## Key thrusts (Jan 18-26)

* Scaling (Mauro + Henry + Anand)
	* If scaling isn't perfect, why and how to improve

* Optimizing performance (Anand + Henry)
	* Ensuring the current kernel implementations are "as good as we can reasonably do" (Nsight Compute tests)
	* MPI halo transfers (CUDA-aware MPI?), can these be faster?
	* Parallel I/O: How does it impact scaling, should we do asynchronous ACC here?, etc.

* Adding features
	* Sub-grid bubbles. Start with classes, then QBMM. (Anand)
	* Hypoelasticity (Jean)
	* Phase change (Mauro)
	* Chemistry (Esteban)

