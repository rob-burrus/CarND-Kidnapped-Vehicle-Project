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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  num_particles = 100;
  min_yr = 0.0001;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  default_random_engine gen;
  
  for(int i =0; i < num_particles; i++){
    Particle p = {
      i, //id
      dist_x(gen), //x
      dist_y(gen), //y
      dist_theta(gen), //theta
      1.0 //weight
    };
    particles.push_back(p);
    weights.push_back(1.0);
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  
  //add gaussian noise by sampling froma  gaussian distribution with mean = updated particle position and std = std of measurements
  default_random_engine gen;
  
  for(int i = 0; i < num_particles; i++){
    
    Particle& p = particles[i];
    
    if(fabs(yaw_rate) > min_yr){
      p.x = p.x + (velocity/yaw_rate)*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
      p.y = p.y + (velocity/yaw_rate)*(cos(p.theta) - cos(p.theta  + yaw_rate*delta_t));
      p.theta = p.theta + yaw_rate*delta_t;
    }else{
      p.x = p.x + velocity * delta_t * cos(p.theta);
      p.y = p.y + velocity * delta_t * sin(p.theta);
    }
    
    normal_distribution<double> dist_x(p.x, std_pos[0]);
    normal_distribution<double> dist_y(p.y, std_pos[1]);
    normal_distribution<double> dist_theta(p.theta, std_pos[2]);
    
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    
    
  }
  
  
  
  

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


  

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
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
  
  
  double weights_sum = 0.0;
  for(int i = 0; i < num_particles; i++){
    
    Particle& p = particles[i];
    
    
    
    for(int j = 0; j < observations.size(); j++){
      
      //transform the observation from vehicle coordinates to map coordinates, as if the observation had been recorded by the particle instead of the vehicle
      LandmarkObs& particle_ob = observations[j];
      LandmarkObs map_ob;
      map_ob.x = particle_ob.x * cos(p.theta) - particle_ob.y * sin(p.theta) + p.x;
      map_ob.y = particle_ob.x * sin(p.theta) + particle_ob.y * cos(p.theta) + p.y;
      
      
      //associate each observed landmark with the nearest map landmark (must be within sensor range)
      Map::single_landmark_s nearest_landmark;
      double nearest_distance = sensor_range;
      bool landmark_found = false;
      for(int k = 0; k < map_landmarks.landmark_list.size(); k++){
        
        Map::single_landmark_s lm = map_landmarks.landmark_list[k];
        
        double distance = dist(lm.x_f, lm.y_f, map_ob.x, map_ob.y);
        if(distance < nearest_distance){
          nearest_distance = distance;
          nearest_landmark = lm;
          landmark_found = true;
        }
        
      }
      
      
      //udpate particle weight with the measurement probability for the nearest landmark, if a landmark was found
      if(landmark_found){
        const double x_i = map_ob.x;
        const double y_i = map_ob.y;
        const double x_mu = nearest_landmark.x_f;
        const double y_mu = nearest_landmark.y_f;
        const double x_diff = pow(x_i - x_mu, 2) / pow(std_landmark[0],2);
        const double y_diff = pow(y_i - y_mu, 2) / pow(std_landmark[1],2);
      
        const double measurementProbability = 1 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]) * exp(-(x_diff + y_diff)/2);
        p.weight *= measurementProbability;
      }

      
    }

    weights[i] = p.weight;
    weights_sum += p.weight;
    
  }
  
  //normalize weights
  for( int i = 0; i< num_particles; i++){
    particles[i].weight /= weights_sum;
  }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  default_random_engine gen;
  discrete_distribution<> weight_distribution(weights.begin(), weights.end());
  
  vector<Particle> resampled_particles;
  
  for (int i = 0; i < num_particles; i++) {
    resampled_particles.push_back( particles[ weight_distribution(gen) ] );
  }
  
  particles = resampled_particles;
  

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
