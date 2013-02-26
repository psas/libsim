/**
 * Main
 *
 * Author: Nathan Bergey
 */

#include <stdio.h>
#include <stdlib.h>
#include "libsim_types.h"
#include "libsim.h"

/**
 * main
 */
int main(int argc, char **argv)
{
	// libsim.stage
	stage sustainer;

	sustainer = Init_Stage();
	
	return 0; //exit
}
