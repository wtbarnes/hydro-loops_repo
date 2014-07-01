#include "hydro-loops.h"

int main(int argc, char *argv[])
{
	//Read in command line arguments
	if(argc != 4)
	{
		printf("Incorrect number of input arguments. Exiting program.\n");
		return 1;
	}
	else
	{
		double L = atof(argv[1]);			//Read in loop half-length
		double Sh = atof(argv[2]);			//Read in heating scale height
		int heat_key = atof(argv[3]);		//Choose (0) uniform or (1) non-uniform heating
	}
	
	//DEBUG--print command line arguments
	printf("The loop half-length is %f Mm",L);
	printf("The heating scale height is %f Mm",Sh);
	printf("The heat key is %d",heat_key);
	
	//Declare variables
	
	//Exit with no errors
	return 0;
}