// celestial_mechanics.cc
//
// author: Corbin Foucart 2022, MIT
// email: corbin.foucart@gmail.com
// license: BSD
// Please feel free to use and modify this code, but please keep the above
// information. Thanks!
//

#include "celestial_mechanics.h"

namespace celestial {

//=============================================================================
// Vector
//=============================================================================
template<int dim>
Vector<dim>::Vector()
{ }


template<int dim>
double& Vector<dim>::operator[](unsigned int i)
{
  return x[i];
}


template<int dim>
void Vector<dim>::operator+=(double a)
{
  for (unsigned int d=0; d < dim; ++d)
    x[d] += a;
}


template<int dim>
void Vector<dim>::operator*=(double a)
{
  for (unsigned int d=0; d < dim; ++d)
    x[d] *= a;
}


template<int dim>
void Vector<dim>::operator/=(double a)
{
  for (unsigned int d=0; d < dim; ++d)
    x[d] /= a;
}


template<int dim>
void Vector<dim>::operator+=(Vector<dim>& other)
{
  for (unsigned int d=0; d < dim; ++d)
    x[d] += other[d];
}


// operations returning copies
template<int dim>
Vector<dim> operator-(Vector<dim>& v1, Vector<dim>& v2)
{
  Vector<dim> r;
  for (unsigned int d = 0; d < dim; ++d)
    r[d] = v1[d] - v2[d];
  return r;
}


template<int dim>
Vector<dim> operator*(Vector<dim>& v1, double a)
{
  Vector<dim> r;
  for (unsigned int d = 0; d < dim; ++d)
    r[d] = v1[d] * a;
  return r;
}


// template specializations for computation of norm
template<>
double Vector<1>::norm()
{ return x[0]; }


template<>
double Vector<2>::norm()
{
  return sqrt(x[0]*x[0] + x[1]*x[1]);
}


template<>
double Vector<3>::norm()
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

//=============================================================================
// Body
//=============================================================================
template<int dim>
Body<dim>::Body(double mass)
  : mass(mass)
{ }


template<int dim>
void Body<dim>::update_acceleration()
{
  for (unsigned int d=0; d < dim; ++d)
    acceleration[d] = force[d] / mass;
}

template<int dim>
void Body<dim>::write_csv_header(std::string filename)
{
  std::filesystem::remove(filename);

  std::ofstream out;
  out.open(filename);
  out << "time,";
  for (auto& s : {"x", "v", "a"})
    for (unsigned int d = 0; d < dim; ++d)
      out << s << d << ",";

  out << "mass\n";
  out.close();
}


template<int dim>
void Body<dim>::output_data_to_csv(std::string filename, const double time)
{
  std::ofstream out;
  out.open(filename, std::ios::app);
  out << time << ",";

  for (unsigned int d = 0; d < dim; ++d)
    out << position[d] << ",";

  for (unsigned int d = 0; d < dim; ++d)
    out << velocity[d] << ",";

  for (unsigned int d = 0; d < dim; ++d)
    out << acceleration[d] << ",";

  out << mass << "\n";
  out.close();
}

//=============================================================================
// NBodySimulation
//=============================================================================
template<int dim>
NBodySimulation<dim>::NBodySimulation(
    std::vector<Body<dim>> body_initial_conditions,
    std::string output_path)
  : bodies(body_initial_conditions)
  , output_dir(output_path)
  , current_time(0)
{ }


template<int dim>
void NBodySimulation<dim>::initialize_data_files()
{
  if (output_dir != "")
    std::filesystem::create_directory(output_dir);

  unsigned int n_bodies = bodies.size();
  for (unsigned int i = 0; i < n_bodies; ++i)
  {
    std::string csv_file = output_dir + "b" + std::to_string(i) + ".csv";
    bodies[i].write_csv_header(csv_file);
    bodies[i].output_data_to_csv(csv_file, current_time);
  }
}


template<int dim>
void NBodySimulation<dim>::update_forces()
{
  unsigned int n_bodies = bodies.size();

  // clear all old forces
  for (unsigned int i = 0; i < n_bodies; ++i)
    bodies[i].force *= 0;
    
  // compute new interactions
  for (unsigned int i = 0; i < n_bodies; ++i)
    for (unsigned int j = i + 1; j < n_bodies; ++j)
    {
      add_interaction_forces(bodies[i], bodies[j]);
    }
}


template<int dim>
void NBodySimulation<dim>::update_accelerations()
{
  for (auto& body : bodies)
    body.update_acceleration();
}


template<int dim>
void NBodySimulation<dim>::integrate_time_step(double dt)
{
  update_forces();
  update_accelerations();

  for (auto& body : bodies)
    integrator.integrate(body, dt);

  current_time += dt;
}


template<int dim>
void NBodySimulation<dim>::run(
    const double dt, 
    const unsigned int N_timesteps)
{
  for (unsigned int i = 0; i < N_timesteps; ++i)
  {
    integrate_time_step(dt);
    output_to_csv();
  }
}


template<int dim>
void NBodySimulation<dim>::output_to_csv()
{
  unsigned int n_bodies = bodies.size();
  for (unsigned int i = 0; i < n_bodies; ++i)
  {
    std::string csv_file = output_dir + "b" + std::to_string(i) + ".csv";
    bodies[i].output_data_to_csv(csv_file, current_time);
  }
}


// explicit instantiations
template class Vector<2>;
template class Vector<3>;

template class Body<2>;
template class Body<3>;

template class NBodySimulation<2>;
template class NBodySimulation<3>;

} // end namespace celestial
