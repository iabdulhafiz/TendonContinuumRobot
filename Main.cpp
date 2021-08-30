
#include "ForwardDynamics.h"


int main(int, char**) {
	
	ContinuumRobotLibrary::ForwardDynamics CR;

	std::vector<double> q = { 0, 0 };
	for (int i = 0; i < 400; i++) {

		CR.Next(q);
		q[0] = -0.01619 * i / 400.0;
		q[1] = -0.00619 * i / 200.0;
	}

	return 1;
}