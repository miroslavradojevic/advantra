package advantra.trace;

public interface TinyBranch {

	/*
	 * TRACE PARAMETERS (will probably be the same for each BranchTrace inherited class instance)
	 */
	
	static int 		N 				= 100;					// maximum number of estimations used for tracing one branch
	
	static int 		N_orientations 	= 50;
	static double 	radius_init 	= 2.0;
	static double 	radius_step		= 0.5;
	static double 	radius_limit	= 5.0;
	
	public static double 	k				= 2.0; 			// k*					radius is the total radius of the hypothesis
	public static double	jump_ahead		= 1.2; 			// jump_ahead*			radius is the jump ahead when tracing
	
	static double	radius_std 		= 1.0; // voxel
	static double	direction_std	= 0.5; // rad
	
	public static double 	check_bifurcations = 1.0; 		
	// check_bifurcations*hypothesis_radius marks the spherical distance at which bifurcations are checked
	//static double 			check_bifurcations_start = 0.95;// take those from the shell 0.95*3.0*radius - 1.0*3.0*radius
	
	// likelihood estimation
	static int 		likelihood_cylinder_samples = 150;
	
	// mean-shift
	static int		number_of_convergence_points 	= 100;
	static int 		threshold_convergence_points 	= 20;
	
	// mean-shift detection 3d
	static int 		extract_sphere_resolution 		= 32;
		
	// mean-shift in 2d 
	static int 		extract_plane_resolution 		= 32;
	
}
