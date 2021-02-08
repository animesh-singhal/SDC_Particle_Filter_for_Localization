/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <map>
#include <string>
#include <vector>

#include <random>  //Need this to access 

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::cout;
using std::endl;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  

  // Setting GPS provided state of the car
  // These will act as the mu for the gaussian distribution in which particles will fall
  double gps_x = x;
  double gps_y = y;
  double gps_theta = theta;
  
  // TODO: Set standard deviations for x, y, and theta
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2]; 

  // Creating normal (Gaussian) distributions 'dist_x' for x, sim for y & theta
  normal_distribution<double> dist_x(gps_x, std_x);
  normal_distribution<double> dist_y(gps_y, std_y);
  normal_distribution<double> dist_theta(gps_theta, std_theta);

  // Initializing random engine "gen"  
  std::default_random_engine gen;

  //Creating particles  
  for (int i =0; i<num_particles; i++){
  
    Particle random_robo;   //Refer Particle struct in particle_filter.h (Line#16)

    random_robo.id = i;
    random_robo.x = dist_x(gen);     // Passing random engine gem in Gaussian dist_x. Gives a value of x within the spread 
    random_robo.y = dist_y(gen);     // We call this process: Sampling from a gaussian distribution
    random_robo.theta = fmod(dist_theta(gen) , 2*M_PI );
    random_robo.weight = 1;
    
    particles.push_back(random_robo);  //Refer 'vector<Particle> particles' in header file (Line#112)
    }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {

  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */


  // Here, std_pos uncertainities are sigmas in 3 dimensions: [x [m], y [m], theta [rad]]
  // We could also get uncertainities in velocities & yaw rate. Here we don't have that info
    
  // Initializing random engine "gen". Will be used to pick a sample fromm our gaussians
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);   // Creating gaussian for x
  normal_distribution<double> dist_y(0, std_pos[1]);   // Creating gaussian for y
  normal_distribution<double> dist_theta(0, std_pos[2]);   // Creating gaussian for theta

  
  for (unsigned int i =0; i<particles.size(); i++){
        //Note: particles[i] is member of 'vector<Particle> particles'. Each item is a struct having multiple features

    double theta = particles[i].theta;
    
    if(fabs(yaw_rate)<0.00001){
    	particles[i].x +=  velocity * delta_t* cos(theta);
        particles[i].y +=  velocity * delta_t* sin(theta);      
	}
    else{
        particles[i].x += velocity/yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
        particles[i].y += velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
        particles[i].theta += yaw_rate * delta_t;
    }
    
	//Add noise
    particles[i].x += dist_x(gen);
    particles[i].y  += dist_y(gen);
    particles[i].theta  += dist_theta(gen);
      
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs>& transformed_obs, 
                                     vector<LandmarkObs> map_landmark) {
  /**
   * TODO: Find the landmark that is closest to each 
   *   transformed_observation and assign the particular landmark to this 
   *   observed measurement.
   *   OBJECTIVE: assign an id of a landmark to the transformed_observation id variable which it maps to 
   */

  for(unsigned int i =0; i< transformed_obs.size(); i++){
    double x_obs = transformed_obs[i].x;
    double y_obs = transformed_obs[i].y;
    
    
    //We'll now make an initial assignment to distance variable 
    // This is basically the Initial assignment of a landmark to our observation: 
    //  Then we'll optimize this in the loop 
    
    int id_lm = map_landmark[0].id;   // Extracting landmark's id
    transformed_obs[i].id = id_lm;     // Assigning it to transformed_observation
    
    double x_lm = map_landmark[0].x;  
    double y_lm = map_landmark[0].y;      
         
    double dx = x_obs - x_lm;
    double dy = y_obs - y_lm;
    double distance_init = sqrt( dx*dx + dy*dy );

    //Now, looping over remainig landmarks for a precise allotment. Starting loop from 2nd element
    for(unsigned int j =1; j< map_landmark.size(); j++){
      
      //We'll reuse the above defined variables
      id_lm = map_landmark[j].id;
      x_lm = map_landmark[j].x;
      y_lm = map_landmark[j].y;      
            
      dx = x_obs - x_lm;
      dy = y_obs - y_lm;
      double distance_new = sqrt( dx*dx + dy*dy );
      
      //If the latest landmark in consideration is closer to the sensor reading, update the id 
      if( distance_new < distance_init ){
      	transformed_obs[i].id = id_lm;   // Updating the id
        distance_init = distance_new;    // Updating the distance variable 
      }
    }
    
    
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &sensor_observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */  
  

  //For each particle, convert all the observations to map frame 
  for (unsigned int i =0; i<particles.size(); i++){
    cout<<"\n\nParticle #"<<i<<endl; 
    
    
    //Creating a copy of sensor measuremens wrt the particle
    vector<LandmarkObs> observations = sensor_observations;

    
    //Note: particles is this vector -> 'vector<Particle> particles'
    
    // Finding the current state of our particle
    double x_particle = particles[i].x;
    double y_particle = particles[i].y;  
    double theta_particle = particles[i].theta;
    
    //Converting each sensor observation (relative to vehicle) to map frame. I'll overwrite because it'll get new obs. after dt
    for(unsigned int j =0; j<observations.size(); j++){
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      
      // Performing Homogeneous Transformation
      // transform to map x coordinate
      double x_wrt_map;
      x_wrt_map = x_particle + (cos(theta_particle) * x_obs) - (sin(theta_particle) * y_obs);

      // transform to map y coordinate
      double y_wrt_map;
      y_wrt_map = y_particle + (sin(theta_particle) * x_obs) + (cos(theta_particle) * y_obs);
      
      //Overwrite the transformed landmark positon in map frame OVER position of landmarks in vehicle frame 
      observations[j].x = x_wrt_map;
      observations[j].y = y_wrt_map;
      
    }

    //Accessing elements of map landmarks. For each particle, we select only the relevent landmarks which are in range of sensor 
    //cout<<"Transformation done"<<endl;
    vector<LandmarkObs> relevant_lm;
    
    int n = map_landmarks.landmark_list.size();
    for (int k=0; k<n; k++){      
      double lm_id = map_landmarks.landmark_list[k].id_i;
      double lm_x = map_landmarks.landmark_list[k].x_f;
      double lm_y = map_landmarks.landmark_list[k].y_f;      
      
      float dx = x_particle - lm_x;
      float dy = y_particle - lm_y;
      float distance = sqrt( dx*dx + dy*dy );
      
      if( distance <= sensor_range){
       LandmarkObs valid_lm;
       valid_lm.id = lm_id;
       valid_lm.x = lm_x;
       valid_lm.y = lm_y;
       
       relevant_lm.push_back(valid_lm);
      }
     }
    
    if (relevant_lm.size() == 0){
     particles[i].weight = 0;
       continue;
    }
    
    //Associate each transformed sensor observation to the nearest landmark (within range)
    dataAssociation(observations, relevant_lm);
    //cout<<"Association successful"<<endl;
    
    particles[i].weight= 1;  //Re-initialise the weight for this iteration for the current particle
    
    // Loop over all sensor observations
    for(unsigned int k =0; k<observations.size(); k++){
      double x_obs = observations[k].x;
      double y_obs = observations[k].y;
      double id_obs = observations[k].id;
      
      //getting the matching landmark for this sensor observation
      int j =0;
      while(relevant_lm[j].id != id_obs){
      	j++;  
      }
      double x_lm = relevant_lm[j].x;
      double y_lm = relevant_lm[j].y;  
        
      float sig_x =  std_landmark[0];  
      float sig_y = std_landmark[1];
      
      cout<<k<<" th association probability: "<<multiv_prob(sig_x, sig_y, x_obs, y_obs, x_lm, y_lm)<<endl;
      //cout<<"Sensor: ("<<x_obs<<", "<<y_obs<<")\nLandmark: ("<<x_lm<<", "<<y_lm<<")\n\n";

      
      particles[i].weight*= multiv_prob(sig_x, sig_y, x_obs, y_obs, x_lm, y_lm);  
      // this update is within a for loop which loops over all the sensor observations
    } 
    cout<<"Final Particle weight: "<<particles[i].weight<<endl;
  }  
 //cout<<endl<<endl; 
}


double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;

}



void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  //We had N number of particles. Here we'll make a new set of N particles. 
  
  
  // A new vector which'll temporarily hold new resampled particles
  vector<Particle> new_particles;
  
  // Making a vector of weights which'll go in discrete_distribution

  vector<float> weights;
  
  //float sum=0;
  for (unsigned int i =0; i< particles.size(); i++) {  
    //cout<<"Weights in Resample:"<<particles[i].weight<<endl;
    weights.push_back( particles[i].weight );
    //sum += particles[i].weight;
  }
  
  
  
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  
  for (unsigned int i =0; i< particles.size(); i++){
  	int num = d(gen);
    new_particles.push_back(particles[num]);
  }

  particles = new_particles; 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}



/*
void ParticleFilter::dataAssociation(vector<LandmarkObs>& transformed_obs, 
                                     vector<LandmarkObs> map_landmark) {


  for(unsigned int i =0; i< transformed_obs.size(); i++){
    double x_obs = transformed_obs[i].x;
    double y_obs = transformed_obs[i].y;
    
    
    //We'll now make an initial assignment to distance variable 
    // This is basically the Initial assignment of a landmark to our observation: 
    //  Then we'll optimize this in the loop 
    
    int id_lm = map_landmark[0].id;   // Extracting landmark's id
    transformed_obs[i].id = id_lm;     // Assigning it to transformed_observation
    
    double x_lm = map_landmark[0].x;  
    double y_lm = map_landmark[0].y;      
         
    double dx = x_obs - x_lm;
    double dy = y_obs - y_lm;
    double distance_init = sqrt( dx*dx + dy*dy );

    //Now, looping over remainig landmarks for a precise allotment. Starting loop from 2nd element
    for(unsigned int j =1; j< map_landmark.size(); j++){
      
      //We'll reuse the above defined variables
      id_lm = map_landmark[j].id;
      x_lm = map_landmark[j].x;
      y_lm = map_landmark[j].y;      
            
      dx = x_obs - x_lm;
      dy = y_obs - y_lm;
      double distance_new = sqrt( dx*dx + dy*dy );
      
      //If the latest landmark in consideration is closer to the sensor reading, update the id 
      if( distance_new < distance_init ){
      	transformed_obs[i].id = id_lm;   // Updating the id
        distance_init = distance_new;    // Updating the distance variable 
      }
    }
    
    
  }
}
*/



/*

    // Loop over all sensor observations
    for(unsigned int k =0; k<observations.size(); k++){
      double x_obs = observations[k].x;
      double y_obs = observations[k].y;
      double id_obs = observations[k].id;
      
      //getting the matching landmark for this sensor observation
      int j =0;
      while(relevant_lm[j].id != id_obs){
      	j++;  
      }
      double x_lm = relevant_lm[j].x;
      double y_lm = relevant_lm[j].y;  
        
      float sig_x =  std_landmark[0];  
      float sig_y = std_landmark[1];
      
      cout<<k<<" th association probability: "<<multiv_prob(sig_x, sig_y, x_obs, y_obs, x_lm, y_lm)<<endl;
      //cout<<"Multi_v_probab"<<multiv_prob(sig_x, sig_y, x_obs, y_obs, x_lm, y_lm)<<endl;
      particles[i].weight*= multiv_prob(sig_x, sig_y, x_obs, y_obs, x_lm, y_lm);
      // this update is within a for loop which loops over all the sensor observations
    } 


*/