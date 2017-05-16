/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
using std::normal_distribution;
using std::default_random_engine;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	//  Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 200;
  	weights.resize(num_particles, 1.0f);

	//Create a normal Gaussian
	//default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);
	Particle prt;
     for (int i = 0; i < num_particles; i++) {
		double sample_x, sample_y, sample_psi;

		sample_x = dist_x(gen);
		sample_y = dist_y(gen);

		sample_psi =dist_psi(gen);

		//store the inital id, x,y,theta, and weight for each particle
		prt.id = i;

		prt.x = sample_x;

		prt.y = sample_y;

		prt.theta = sample_psi;

		prt.weight = 1.0f;

		particles.push_back(prt);

	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	//  Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	static default_random_engine gen;
	double xf, yf, thetaf;
	for (int i =0; i < num_particles; i++) {
		if (yaw_rate != 0) { 
			xf = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			yf = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			thetaf = particles[i].theta + yaw_rate * delta_t;
		}
		else {
			xf = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			yf = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			thetaf = particles[i].theta;
		}
	
		//Create a normal Gaussian
		normal_distribution<double> dist_x(0, std_pos[0]);
		normal_distribution<double> dist_y(0, std_pos[1]);
		normal_distribution<double> dist_psi(0, std_pos[2]);
		//adding noise
		xf = xf + dist_x(gen);
		yf = yf + dist_y(gen);
		thetaf = thetaf + dist_psi(gen);
		//store the x, y, theta for each particle
		particles[i].x = xf;
		particles[i].y = yf;
		particles[i].theta = thetaf;
	} 
}

//void
std::vector<LandmarkObs> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	//  Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	/* Nearest Neighbor:find the least distance from obs to landmark and choose the closest landmark with the least distance*/

	double dist, min_dist;
	double delta_x, delta_y;
	int min_id;
	std::vector<LandmarkObs> Lmin_dist;

	for (int i = 0; i < observations.size(); i++) {
		//assign min index to 0 at the beginning of each iteration	
		min_id = 0;
		for (int j = 0; j < predicted.size(); j++) {
			delta_x = observations[i].x - predicted[j].x;
			//std::cout << "dataAssociations here2, delta_x:  " << std::endl << delta_x << std::endl;

			delta_y = observations[i].y - predicted[j].y;

			dist = sqrt(delta_x * delta_x + delta_y * delta_y);
			//initialize minimum distance
			if (j == 0) {
				min_dist = dist;
				min_id = 0;
			}
			else {
				if (dist < min_dist) {
					min_dist = dist;
					min_id = j;
				}
			}
		} //endfor j loop


		Lmin_dist.push_back(predicted[min_id]);

	} //endfor i loop
	return Lmin_dist;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	//  Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	//(x,y) measurement transformation for observation
	double OBSx, OBSy;
	double p, pf;
	double sigmax = 0.3;
	double sigmay = 0.3;
	std::vector<LandmarkObs> L; 
		
	for (int i = 0; i < particles.size(); i++) {
		//assign p and pf = 1 for each iteration
		p = 1.0f;
		pf = 1.0f;

     		// find the closest landmark(L)
		//convert all Map type landmark(x,y) into LandmarkObs type landmark(x,y)
		std::vector<LandmarkObs> pred_landmarks;
		for (auto lm : map_landmarks.landmark_list) {
			LandmarkObs lmp;
			lmp.x = lm.x_f;
			lmp.y = lm.y_f;
			lmp.id = lm.id_i;
			
			//check if the sensor is within range
			double diffx = lmp.x - particles[i].x;
			double diffy = lmp.y - particles[i].y;
			if ((diffx * diffx + diffy * diffy) <= (sensor_range * sensor_range)) {
				pred_landmarks.push_back(lmp);
			}
			
		}

		std::vector<LandmarkObs> obs_landmarks;
		//transformation for observations[i]
		for (int j = 0; j < observations.size(); j++) {
			LandmarkObs olm;
			OBSx = particles[i].x + 
				observations[j].x * 
				cos(particles[i].theta) 
				- observations[j].y * 
				sin(particles[i].theta);

			OBSy = particles[i].y + 
				observations[j].x * 
				sin(particles[i].theta) 
				+ observations[j].y * 
				cos(particles[i].theta);

				olm.x = OBSx;

				olm.y = OBSy;
				olm.id = j;
			obs_landmarks.push_back(std::move(olm));

		}

		L = dataAssociation(pred_landmarks, obs_landmarks);
		for (int k = 0; k < L.size(); k++) {
			//observation multivariate Gaussian probability(p)
			p = (1/(2*3.14*sigmax*sigmay))*
				exp(-1*((L[k].x - obs_landmarks[k].x)*
				(L[k].x -obs_landmarks[k].x)/
				(2 * sigmax * sigmax) 
				+ (L[k].y - obs_landmarks[k].y)*
				(L[k].y-obs_landmarks[k].y)/
				(2 * sigmay * sigmay)));

			//initialize final probability(pf) when k = 0
			if (k == 0) {
				pf = p;
				
			}
			else {

				pf *= p; //that is pf * p
				
			}
		}
		//total weight
		particles[i].weight = pf;
		
		weights[i] = pf; 

	} 
	
}

void ParticleFilter::resample() {
	//  Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//static default_random_engine gen;
	std::discrete_distribution<int> disc_dist(weights.begin(), weights.end());
	std::vector<Particle> nparticles;
	for (int i = 0; i < num_particles; i++) {
		int index = disc_dist(gen);
		nparticles.push_back(std::move(particles[index]));
	}
	particles = std::move(nparticles);

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
