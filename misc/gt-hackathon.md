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
		* This was fixed by Anand 
	* **Merging phase-change with hypoelasticity branches**
	* Create regression tests under `/tests` 
	* Merge master -> phase-change

* Anand
	* Implement Order-3 RK Time stepper
		- Done (RK2, RK3)
	* Parallel IO (with Nsight Systems), say every 100 time steps or so
		- Done (Parallel IO takes less than 2 timesteps)
	* Do some preliminary scaling runs on Summit (starting with weak scaling)
		* Do a 1D and a 3D case, two alphas/components if possible
		- Done (Scales well in 3D upto 500 GPUs)
	* Note: Can't merge GPU->master branch without time-steppers and multi-components verified working
		- Done (Verified timesteppers as well as multi-component)

* Jean
	* Post to Slack the results of scaling on XSEDE machines
		* Done
	* Mention to Jose that he can join and 'submerge' himself into how this works (no key thrusts for him, though)
		* Done
	* Create regression tests under `/tests` for hypoelastic cases
	* Merge Master -> Hypoelastic

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

## In-between thrusts

Anand
* We need regression tests under `/tests` that fully exercise the code piece by piece
	* Generate these "golden" results for comparison using best known working code. Can probably just come from the MFC-develop master branch. Can double check WENO3 against a much older version of the code just in case.
	* Cases we should have: WENO1, WENO3, WENO5, mp_weno, mapped_weno, HLL, HLLC, 1 component, 2 comp., 3 comp., kdivu (aka altsoundspeed), 1D, 2D, 3D, CPU parallel (4 mpi ranks), CPU serial, etc.
	* Be careful if some of these tests like kdivu are not available on GPU yet. Perhaps skip that test in such cases?

Henry
* Finish GPU build system and merge into master


## Key thrusts (Jan 18 + 24-26)

* Optimizing performance (Anand + Henry)
	1. Ensuring the current kernel implementations are "as good as we can reasonably do" (Nsight Compute tests)
		* Another thing to consider: Does assigning more than 1 MPI rank per GPU help speed-up? We found in a previous hackathon that doing this and enabling `gpumps` gave an additional bump in speed-up over the CPU-only case. Can ask Jean/Mauro about this as well. You enable `gpumps` via a `bsub` or `jsrun` command line flag. It should be in the Summit user guide.  
		* Should we be using MACROS!! for certain arrays and/or branch statements in kernels like Riemann/WENO? That are expanded by cpp c pre processor at compile time??
	3. MPI halo transfers (CUDA-aware MPI?), can these be faster?
		* Conclusion from preliminary analysis: This is unimportant because more physics means less important communication time
	4. Parallel I/O: How does it impact scaling, should we do asynchronous ACC here?, etc.
		* Conclusion from preliminary analysis: This is unimportant

* Adding features/physics
	* Sub-grid bubbles. Start with classes, then QBMM. (Anand)
	* Hypoelasticity (Jean)
	* Phase change (Mauro)
	* Chemistry (Esteban)

