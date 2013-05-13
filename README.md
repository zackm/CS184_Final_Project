CS184 Final Project
===================
###Fluid Simulation using Smoothed Particle Hydrodynamics
###Tyler Brabham and Zack Mayeda



[Procedure and Results Writeup](https://github.com/zackm/CS184_Final_Project/blob/master/final%20images/final%20paper.pdf)

[Demo Video](http://youtu.be/oogV7GYl630)

[Temporary Demo Website](http://www-inst.eecs.berkeley.edu/~cs184-bg/final.html)


****Instructions****

*Command-Line Arguments:*

* -help : List all options
* -s # : Load scene number #
	*  1 - Dam Break
	*  2 - Drop
	*  3 - Evenly Distributed Across Scene, with Horizontal Velocity
	*  4 - No Gravity
	*  5 - Collision
	*  Default - Horizontal Velocity
* -m : Export Still Frames for movie into Images/
* -rm : Export Raytraced Still Frames for movie into Multi_Trace/output_pics/
* -i : Render Isosurface
* -p : Render Particles
* -v : Increase viscosity and surface tension
* -h : Hide Display Window
* -delay # : delay movie output by # frames

*Live Keyboard Commands:*

* 'i' : toggle isosurface/particles
* 'r' : output single raytraced image and exit
* 'f' : output current FPS to command line
* 'p' : save current image, not raytrace, and continue
* spacebar : pause
* 'q' : quit