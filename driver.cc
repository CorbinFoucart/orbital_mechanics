// driver.cc
//
// author: Corbin Foucart 2022, MIT
// email: corbin.foucart@gmail.com
// license: BSD
// Please feel free to use and modify this code, but please keep the above
// information. Thanks!
//

#include <iostream>
#include <vector>
#include <unordered_map>

#include "celestial_mechanics.h"

//=============================================================================
// Orbital motion between a larger and smaller mass; system moves up and to the
// left, showing conservation of momentum from the initial condition
//=============================================================================
void ex1_two_body_mom_cons()
{
  celestial::Body<2> sun(10);
  celestial::Body<2> moon(1);
  moon.position[0] = 1;
  moon.velocity[0] = -0.75;
  moon.velocity[1] = 0.75;

  std::vector<celestial::Body<2>> bodies{sun, moon};
  celestial::NBodySimulation<2> simulation(bodies, "ex1_two_body_mom_cons/");
  simulation.initialize_data_files();

  const double dt = 0.001;
  const double N_timesteps = 5000;
  simulation.run(dt, N_timesteps);
}


//=============================================================================
// Similar to ex1, except now we add another moon to cancel out the momentum of
// the first one, resulting in a stationary system
//=============================================================================
void ex2_three_body_stationary()
{
  celestial::Body<2> sun(10);

  celestial::Body<2> moon(1);
  moon.position[0] = 1;
  moon.velocity[0] = -1;
  moon.velocity[1] = 1;

  celestial::Body<2> moon_2(1);
  moon_2.position[0] = -1;
  moon_2.velocity[0] =  1;
  moon_2.velocity[1] = -1;

  std::vector<celestial::Body<2>> bodies{sun, moon, moon_2};
  celestial::NBodySimulation<2> simulation(
      bodies, 
      "ex2_three_body_stationary/");
  simulation.initialize_data_files();

  const double dt = 0.001;
  const double N_timesteps = 5000;
  simulation.run(dt, N_timesteps);
}

//=============================================================================
// Destabilize ex2 by giving the moons unequal momenta
//=============================================================================
void ex3_three_body_destabilized()
{
  celestial::Body<2> sun(10);
  celestial::Body<2> moon(3);
  moon.position[0] = 1;
  moon.velocity[0] = -1;
  moon.velocity[1] = 1;

  celestial::Body<2> moon_2(1);
  moon_2.position[0] = -1;
  moon_2.velocity[0] =  1;
  moon_2.velocity[1] = -3;

  std::vector<celestial::Body<2>> bodies{sun, moon, moon_2};
  celestial::NBodySimulation<2> simulation(
      bodies, 
      "ex3_three_body_destabilized/");
  simulation.initialize_data_files();

  const double dt = 0.001;
  const double N_timesteps = 5000;
  simulation.run(dt, N_timesteps);
}


//=============================================================================
// Symmetric motion that ends up in chaos;
// in this case, due to roundoff and truncation error that builds up over long
// time integration
//=============================================================================
void ex4_circular_bodies(const unsigned int N_bodies)
{
  constexpr double pi = 3.14159265358;
  double mass = 3;

  double theta_0 = pi / 2;
  double delta_theta = 2*pi / N_bodies;

  std::vector<celestial::Body<2>> bodies;
  for (unsigned int i = 0; i < N_bodies; ++i)
  {
    bodies.push_back(celestial::Body<2>(mass));

    const double r = 1;
    const double theta = theta_0 + i * delta_theta;
    const double x = r * cos(theta);
    const double y = r * sin(theta);

    bodies[i].position[0] =  x;  bodies[i].position[1] = y;
    bodies[i].velocity[0] = -y;  bodies[i].velocity[1] = x;
  }

  celestial::NBodySimulation<2> simulation(
      bodies, 
      "ex4_circular_bodies" + std::to_string(N_bodies) + "/");
  simulation.initialize_data_files();

  const double dt = 0.001;
  const double N_timesteps = 13000;
  simulation.run(dt, N_timesteps);
}




int main()
{
  //ex1_two_body_mom_cons();
  //ex2_three_body_stationary();
  //ex3_three_body_destabilized();
  ex4_circular_bodies(4);

  return 0;
}
