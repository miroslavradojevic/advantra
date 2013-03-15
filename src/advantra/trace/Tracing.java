package advantra.trace;

public interface Tracing {

	static int 		N 				= 100;					// maximum number of estimations used for tracing one branch
	static int		CPU_NR 			= 4;
	static int 		N_orientations 	= CPU_NR*20; 
	static int 		MS_PTS			= CPU_NR*40;
	
	static int 		STEPS_BACKWARDS = 5;
	static double	RADIUS_SCALE_BORDER	= 1.5;
	
	static double 	radius_init 	= 1.0;
	static double 	radius_step		= 0.5;
	static double 	radius_limit	= 5.0;

	public static double 	k				= 2.0; 			// k*					radius is the total radius of the hypothesis
	public static double	jump_ahead		= 2.0; 			// jump_ahead*			radius is the jump ahead when tracing
	
	static double	radius_std 		= 2.0; // voxel
	static double	direction_std	= 0.7; // rad
	
	// mean-shift detection 3d
	static int 		extract_sphere_resolution 		= 16;
		
	// mean-shift in 2d 
	static int 		extract_plane_resolution 		= 32;
	
}
