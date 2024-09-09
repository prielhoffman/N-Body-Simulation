#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mpi.h"

///////////////////// parallel C code with MPI to simulate a naive 2D N-body problem, O(n^2)  /////////////////////

#define N 992 // stars
#define num_steps 8500 // number of steps (that gave us at least 2 min for one processor)
#define timestep 2e15 // rate (that gave us at least 2 min for one processor)
#define G 6.674e-11 // Gravitational constant
#define m_sun 2e30 // mass of sun
#define v_avg 200000 // Average speed
#define domain 9e17 //  domain's size
#define PI 3.14159265358979323846 // for random theta
#define start_f "positions_start.txt" // file with the stars positions in the first iteration  
#define middle_f "positions_middle.txt" // file with the stars positions at the middle (num_of_steps/2)
#define end_f "positions_end.txt" // file with the stars positions at the end (of all iterations)

// this struct was made for represent one star
typedef struct {
    double x, y;     // position
    double vx, vy;   // velocity
    double fx, fy;   // force
} Body;

Body galaxy[N]; // this array represent the whole galaxy (all stars)

// this function gets the array of all the stars positions and an index (that represents the part of the simulation- start/middle/end) and saves the values from the array in a txt file
void save_positions(int index) 
{
	char filename[50];
	if (index==1)// index=1 for represent the start (the first iteration)
		sprintf(filename, start_f);
	if (index==2)// index=2 for represent the middle (num_of_steps/2)
		sprintf(filename, middle_f);
	if (index==3)// index=3 for represent the end (the last iteration)
		sprintf(filename, end_f);
	FILE* file = fopen(filename, "w");
	if (file == NULL) // try to create a file
	{
        	printf("Error opening file for writing: %s\n", filename);
        	return;
    	}
	// fill the file withr the stars psitions
	int i;
	for (i = 0; i < N; i++) 
	{
		fprintf(file, "%lf %lf\n", galaxy[i].x, galaxy[i].y); // x and y coordinate
	}

    fclose(file); // close this file
}

// this function generates a random double between min and max
double randomDouble(double min, double max) 
{
    return min + ((double)rand() / RAND_MAX) * (max - min); // generate rand number using rand() function
}

// this function generate initial random positions and velocities, in a square domain.
// velocities in the range 0.5v<v<1.5v with uniform direction distributions
void initialize() 
{
	int i;
    for (i = 0; i < N; i++) // for all the start (N=992)
{
        galaxy[i].x = randomDouble(0.0, 1)*domain; // choose a random x position
        galaxy[i].y = randomDouble(0.0, 1)*domain; // choose a random y position
		double v = randomDouble(0.5, 1.5); // range of velocity: 0.5v<v<1.5v
        double theta = randomDouble(0.0, 2.0 * PI);
        galaxy[i].vx = v * cos(theta); // choose a random velocity in X axis 
        galaxy[i].vy = v * sin(theta); // choose a random velocity in Y axis 
    }
}

// this function compute force on every star (name: body) 
// this calculation is per one rank (need to compute only num_bodies of stars)
void compute_forces(int num_bodies, int rank) 
{
	double dx, dy, distance, force;
	int i, j;
	// compute only for a spesific rank (we will compute the force of only N(number of stars)/size(number of ranks)
	for (i = rank*num_bodies; i < num_bodies*(rank+1); i ++) 
	{	
		galaxy[i].fx = 0.0;
		galaxy[i].fy = 0.0;
		for (j = 0; j < N; j++) // check which stars influence (in terms of force) the specific i star
		{
			if (j != i) // check that is another star
			{	//compute force on ith body like we saw in class
				dx = galaxy[j].x - galaxy[i].x;
				dy = galaxy[j].y - galaxy[i].y;
				distance = sqrt(dx * dx + dy * dy);
				force = (G * m_sun * m_sun) / (distance * distance);
				galaxy[i].fx += force * (dx / distance);
				galaxy[i].fy += force * (dy / distance);
		    	}
		}
    }
}

// this function compute and update the positions and the velocities of all the stars for specific rank
// it means that every rank get chank of  stars (num_bodies=N/size) and need to compute and update the positions and the velocities only for this chank of stars
void update_positions_velocities(int num_bodies, int rank) 
{
	int i,chank_start,chank_end;
	chank_start=rank*num_bodies; // beginning of the chank
	chank_end=chank_start+num_bodies; // end of the chank
    	for (i = chank_start; i < chank_end; i++) // each rank calculate only this chank of stars
	{
		// update velocity & position like we saw in class
		galaxy[i].x += galaxy[i].x * timestep; // new velocity in X axis
		galaxy[i].y += galaxy[i].y * timestep; // new velocity in Y axis
		galaxy[i].x += galaxy[i].vx * timestep; // new position in X axis
		galaxy[i].y += galaxy[i].vy * timestep; // new position in Y axis
		// check if a body exits the domain we will choose new position (in the domain) -> it was up to us to decide what action to take
		if (galaxy[i].x < 0.0) galaxy[i].x=randomDouble(0.0, 1)*domain;
		if (galaxy[i].x > domain) galaxy[i].x=randomDouble(0.0, 1)*domain;
		if (galaxy[i].y < 0.0) galaxy[i].y=randomDouble(0.0, 1)*domain;
		if (galaxy[i].y > domain) galaxy[i].y=randomDouble(0.0, 1)*domain;
	}
}

int main(int argc, char** argv) 
{
	int rank, size, step;
	double start_time, end_time;
	Body *bodies;
	MPI_Init(&argc, &argv); // initializes the MPI execution environment
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // retrieves the rank of the calling process within the communicator
	MPI_Comm_size(MPI_COMM_WORLD, &size); // retrieves the size of the communicator
	if (rank == 0) // "master"
	{
		initialize(); // init all stars (whole galaxy)
		save_positions(1); // need to save the positions.1 is the index for represent the start part
	}
	MPI_Bcast(galaxy, N*sizeof(Body), MPI_BYTE, 0, MPI_COMM_WORLD); // master sent the (whole galaxy) by MPI_Bcast
	int num_bodies = N/size; // (create a chank) divide the work between all ranks
	bodies=(Body *)malloc(num_bodies * sizeof(Body)); // this is the chank of stars per each rank
	if (rank == 0) // if you are the "master" you need to start the timer
	{
		start_time=MPI_Wtime(); // start timer
	}
	// start simulation
	for (step = 0; step < num_steps; step++) 
	{
		// compute forces only per chank (per rank) with all the stars in the galaxy 
		compute_forces(num_bodies, rank); 
		// compute & update positions & velocities only per chank (per rank) 
		update_positions_velocities(num_bodies, rank);
		// copy the chank of stars from the whole galaxy to the chank that we going to send and gather (by MPI_Allgather)
		memcpy(bodies,(galaxy+(num_bodies*rank)),num_bodies*sizeof(Body));
		// gather final positions by MPI_Allgather (from every rank)
		MPI_Allgather(bodies, num_bodies*sizeof(Body), MPI_BYTE, galaxy, num_bodies*sizeof(Body), MPI_BYTE, MPI_COMM_WORLD);
		if(rank == 0 && step == num_steps/2) // if you are the "master" and we pass half of  the total iterations 
		{
			save_positions(2); // need to save the positions. 2 is the index for represent the middle part
		}
	}
	if (rank == 0) // "master"
	{
		end_time=MPI_Wtime(); // if you are the "master" you need to stop the timer
		save_positions(3); // need to save the positions. 3 is the index for represent the finish part
		double elapsed_time = end_time - start_time; // need to calculate the total time of the simulation
        	printf("Done Simulation!!! Elapsed time: %f seconds\n", elapsed_time); // print the total time of the simulation
    	}
	// clean all by free
	free(bodies); // free the specific chank
	MPI_Finalize();
	return 0;
}