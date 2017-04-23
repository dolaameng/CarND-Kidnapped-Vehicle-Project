/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <cassert>

#include "particle_filter.h"

// random generator engine
// Better put it in the class declariation, but the submission
// doesn't allow me to change particle_filter.h
std::default_random_engine random_generator;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	std::normal_distribution<double> random_x(x, std[0]/*=std_x*/);
	std::normal_distribution<double> random_y(y, std[1]/*=std_y*/);
	std::normal_distribution<double> random_theta(theta, std[2]/*=std_theta*/);

	num_particles = 500;
	for (auto i = 0; i < num_particles; ++i) {
		auto x_val = random_x(random_generator),
			 y_val = random_y(random_generator),
			 theta_val = random_theta(random_generator),
			 weight_val = 1.0;
		particles.push_back(Particle{/*id=*/i, x_val, y_val, theta_val, weight_val});
		weights.push_back(weight_val); // BAD CODE: WHAT'S THE POINT OF USING `weights` IF PARTICLE ALREADY HAS A WEIGHT MEMBER??
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	auto & g = random_generator;
	double dt = delta_t,
		   r = (fabs(yaw_rate) > 1e-3) ? (velocity / yaw_rate) : 0.0,
		   vdt = velocity * dt,
		   ydt = yaw_rate * dt;

	transform(
		particles.begin(), particles.end(),
		particles.begin(), 
		[&g, std_pos, dt, r, yaw_rate, ydt, vdt](Particle p){
			auto theta = p.theta + ydt,
				 x = (fabs(yaw_rate) > 1e-3) ? (p.x + r * (sin(theta) - sin(p.theta))) : (p.x + vdt * cos(p.theta)),
				 y = (fabs(yaw_rate) > 1e-3) ? (p.y + r * (cos(p.theta) - cos(theta))) : (p.y + vdt * sin(p.theta));

			std::normal_distribution<double> random_x(x, std_pos[0]/*=std_x*/); 
			std::normal_distribution<double> random_y(y, std_pos[1]/*=std_y*/); 
			std::normal_distribution<double> random_theta(theta, std_pos[2]/*=std_theta*/); 

			return Particle{
				p.id,
				random_x(g),
				random_y(g),
				random_theta(g),
				p.weight
			};}
	);

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// update the ids in observations to make them synchronized with index in predicted
	for (LandmarkObs & ob : observations) {
		auto x = ob.x,
			 y = ob.y;
		std::vector<double> to_landmark_dists(predicted.size());
		transform(
			predicted.begin(), predicted.end(),
			to_landmark_dists.begin(),
			[x, y](const LandmarkObs & landmark) {
				return dist(x, y, landmark.x, landmark.y);
			}
		);
		auto imin = std::min_element(to_landmark_dists.begin(), to_landmark_dists.end()) - to_landmark_dists.begin();
		ob.id = imin;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.) - WHAT???
	//   http://planning.cs.uiuc.edu/node99.html

	// normalizing factor for multivariate gaussian
	double P = 2 * M_PI * std_landmark[0] * std_landmark[1];

	// for each particle, update its weight
	for (Particle & p : particles) {
		auto x = p.x,
			 y = p.y,
			 theta = p.theta;

		// transform observation from p coord to map coord
		std::vector<LandmarkObs> obs_in_map{observations.size()};
		transform(
			observations.begin(), observations.end(),
			obs_in_map.begin(),
			[x, y, theta](LandmarkObs ob) {
				return LandmarkObs{
					/*id=*/ob.id,
					/*x=*/ ob.x * cos(theta) - ob.y * sin(theta) + x,
					/*y=*/ ob.x * sin(theta) + ob.y * cos(theta) + y
				};
			}
		);

		// predicted observations as a subset of landmarks, 
		// which are "observable" from the particle
		std::vector<LandmarkObs> predicted_obs;
		for (const auto & landmark : map_landmarks.landmark_list) {
			bool in_range = (dist(x, y, landmark.x_f, landmark.y_f) <= sensor_range);
			if (!in_range) continue;
			predicted_obs.push_back(LandmarkObs{
				/*id=*/landmark.id_i,
				/*x=*/ landmark.x_f,
				/*y=*/ landmark.y_f
			});
		}

		// associate observations_in_map with observable landmarks,
		// by using nearst neighbor
		// after this, the ids in obs_in_map will be synchronized with predicted_obs.
		dataAssociation(predicted_obs, obs_in_map);

		// assign weights by a multivariate Gaussian
		double weight = 0;
		if (predicted_obs.empty()) { /*cannot find any landmark in range, which means it is a bad particle*/
			weight = 1e-10;
		} else {                     /*find at least one landmark to compare with*/
			weight = 1;
			auto var_x = std_landmark[0] * std_landmark[0],
				 var_y = std_landmark[1] * std_landmark[1];
			for (const auto & ob : obs_in_map) {
				double dx = ob.x - predicted_obs[ob.id].x; // ob.id corresponds to index of predicted in dataAssociation
				double dy = ob.y - predicted_obs[ob.id].y;
				weight *= exp( -dx*dx/(2*var_x) - dy*dy/(2*var_y) ) / P;
			}
		}

		// these are unnormalized weights
		p.weight = weight;
		weights[p.id] = weight; 
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::discrete_distribution<int> sample_dist(weights.begin(), weights.end());
	std::vector<Particle> new_particles(particles.size());
	for (auto i = 0; i < particles.size(); ++i) {
		Particle & p = new_particles[i];
		p.id = i;
		p.weight = 1.;
		int choice = sample_dist(random_generator);
		p.x = particles[choice].x;
		p.y = particles[choice].y;
		p.theta = particles[choice].theta;
	}
	particles = new_particles;
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
