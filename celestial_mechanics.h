// celestial_mechanics.h
//
// author: Corbin Foucart 2022, MIT
// email: corbin.foucart@gmail.com
// license: BSD
// Please feel free to use and modify this code, but please keep the above
// information. Thanks!

#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cmath>

namespace celestial {

// The purpose of this library is to build physical intuition for orbital
// mechanics, so we don't worry about scaling factors such as G or declare our
// bodies to planetary scale. Doing so would require trivial modifications.
constexpr double G = 1;

template<int dim>
class Vector
{
  public:
    Vector();

    double& operator[](unsigned int i);

    void operator+=(double a);
    void operator*=(double a);
    void operator/=(double a);
    void operator+=(Vector<dim>& other);
    double norm();

    double x[dim];
};


template<int dim>
Vector<dim> operator-(Vector<dim>& v1, Vector<dim>& v2);

template<int dim>
Vector<dim> operator*(Vector<dim>& v1, double a);


template<int dim>
struct Body{
    Body(double mass);
    void update_acceleration();

    double mass;
    Vector<dim> position;
    Vector<dim> velocity;
    Vector<dim> acceleration;
    Vector<dim> force;

    void write_csv_header(std::string filename);
    void output_data_to_csv(std::string filename, const double time);
};


// computes force on b due to other, negate for reverse
template<int dim>
Vector<dim> gravitational_force(Body<dim>& b, Body<dim>& other)
{
  Vector<dim> r_hat = other.position - b.position;
  const double r = r_hat.norm();

  r_hat /= r;

  const double F =  G * b.mass * other.mass / std::pow(r, 2);
  r_hat *= F;
  
  // scaled by F to be force vector
  return r_hat; 
}


template<int dim>
void add_interaction_forces(Body<dim>& b1, Body<dim>& b2)
{
  Vector<dim> Fg = gravitational_force(b1, b2);
  b1.force += Fg;

  Fg *= -1;
  b2.force += Fg;
}


template<int dim>
std::ostream& operator<<(std::ostream& os, Vector<dim>& v)
{
  for (unsigned int i = 0; i < dim; ++i)
    os << v[i] << std::endl;
  return os;
}


template<int dim>
class TimeIntegrator
{
  public:
    virtual void integrate(Body<dim>& body, const double dt) = 0;
};


template<int dim>
class EulerCromerIntegrator : public TimeIntegrator<dim>
{
  public:
    virtual void integrate(Body<dim>& body, const double dt) override
    {
      auto delta_v =  body.acceleration * dt;
      body.velocity += delta_v;

      auto delta_x = body.velocity * dt;
      body.position += delta_x;
    }
};


template<int dim>
class NBodySimulation
{
  public:
    NBodySimulation(std::vector<Body<dim>> body_initial_conditions,
                    std::string output_path="");

    // writes current state of the simulation to csv
    void initialize_data_files();

    // loops over bodies and updates the resultant forces based on the current
    // positions, then updates the accelerations once all forces have been
    // computed.
    void update_forces();
    void update_accelerations();

    // numerically integrates the equations of motion from t to t + dt
    void integrate_time_step(double dt);
    
    // convenience function to run a whole simulation once
    // initialize_data_files() has been called
    void run(const double dt, const unsigned int N_timesteps);

    // outputs current state to output_dir
    void output_to_csv();

    std::vector<Body<dim>> bodies;
    EulerCromerIntegrator<dim> integrator;

    std::string output_dir;
    double current_time;

};

} // end namespace celestial
