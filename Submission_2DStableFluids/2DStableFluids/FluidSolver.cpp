#include "StdAfx.h"
#include "FluidSolver.h"

//Loosely following Jos Stam's Stable Fluids

CFluidSolver::CFluidSolver(void):
n(60), size(n*n), h(0.1), laplacian(size,size), diffusion(size,size)
{
	//default size is set to 60^2
	velocity = new vec2[size];
	velocity_source = new vec2[size];
	advected_velocity = new vec2[size];
	density = new double[size];
	density_source = new double[size];
	pressure = new double[size];
	divergence = new double[size];

	double diffusion_coef = 0.3*h;
	viscosity_coef = 0.05;
	buoyancy_coef = 2.0;
	//Set up the Laplacian matrix and diffusion matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int index = i + j*n;
			if (i>0 && i<n-1 && j>0 && j<n-1) {
				if (i-1>0) {
					laplacian.set1Value(index, index-1, -1.0);
				}
				if (i+1<n-1) {
					laplacian.set1Value(index, index+1, -1.0);
				}
				if (j-1>0) {
					laplacian.set1Value(index, index-n, -1.0);
				}
				if (j+1<n-1) {
					laplacian.set1Value(index, index+n, -1.0);
				}
				laplacian.set1Value(index, index, 4.);
			} else {
				laplacian.set1Value(index, index, 1.);
			}
			int count = 0;
			if (i-1>0) {
				diffusion.set1Value(index, index-1, -1.0*diffusion_coef);
				count++;
			}
			if (i+1<n-1) {
				diffusion.set1Value(index, index+1, -1.0*diffusion_coef);
				count++;
			}
			if (j-1>0) {
				diffusion.set1Value(index, index-n, -1.0*diffusion_coef);
				count++;
			}
			if (j+1<n-1) {
				diffusion.set1Value(index, index+n, -1.0*diffusion_coef);
				count++;
			}
			diffusion.set1Value(index, index, 1.+count*diffusion_coef);
		}
	}

	reset();
}

void CFluidSolver::reset()
{
	for (int i = 0; i < size; i++) {
		density[i] = 0.;
		density_source[i] = 0.;
		velocity[i] = vec2(0.,0.);
		divergence[i] = 0.;
		pressure[i] = 0.;
		velocity_source[i] = vec2(0.,0.);
		advected_velocity[i] = vec2(0.,0.);
	}

}

CFluidSolver::~CFluidSolver(void)
{
	delete[] velocity;
	delete[] density;
	delete[] pressure;
	delete[] divergence;

	delete[] density_source;
	delete[] velocity_source;
	delete[] advected_velocity;
}

void CFluidSolver::update()
{
	updateDensity();
	updateVelocity();
}

void CFluidSolver::updateDensity()
{
	add(density, density, density_source); // density += density_source;

	//Diffusion process
	diffusion.solve(density_source, density, 1e-8, 30); // Diffusion_matrix density_new = density_old

	density_advection();
	clean_density_source();
}

void CFluidSolver::updateVelocity()
{
	velocity_advection();
	add(velocity, advected_velocity, velocity_source);
	apply_buoyancy();
	velocity_diffusion();

	projection();
	clean_velocity_source();
}

void CFluidSolver::projection()
{
	//set boundary condition
	for (int i=0; i< n; i++) {
		*v(0,i)=vec2(0.,0.);
		*v(n-1,i)=vec2(0.,0.);
		*v(i, 0)=vec2(0.,0.);
		*v(i, n-1)=vec2(0.,0.);
	}

	//compute divergence
	for (int i = 1; i < n-1; i++)
	{
		for (int j = 1; j < n-1; j++)
		{
			divergence[i+n*j] = 0.5*(v(i+1,j)->x-v(i-1,j)->x
				+ v(i, j+1)->y - v(i, j-1)->y);
		}
	}

	//get pressure by solving (Laplacian pressure = divergence)
	laplacian.solve(pressure, divergence, 1e-8, 10);

	//update velocity by (velocity -= gradient of pressure)
	for (int i = 1; i < n-1; i++)
	{
		for (int j = 1; j < n-1; j++)
		{
			v(i, j)->x += 0.5 * (p(i+1, j) - p(i-1, j));
			v(i, j)->y += 0.5 * (p(i, j+1) - p(i, j-1));
			*v(i,j) = *v(i,j)*1.;
		}
	}

}

void CFluidSolver::clean_density_source()
{
	for (int i=0; i < size; i++) {
		density_source[i] = 0.;
	}
}

void CFluidSolver::clean_velocity_source()
{
	for (int i=0; i < size; i++) {
		velocity_source[i] = vec2(0.,0.);
	}
}

void CFluidSolver::density_advection()
{
	//set boundary condition
	for (int i=0; i< n; i++) {
		*d(0,i)= 0;
		*d(n-1,i)=0;
		*d(i, 0)=0;
		*d(i, n-1)=0;
	}

	vec2 backtraced_position;
	//Density advection
	for (int i = 1; i < n-1; i++)
		for (int j=1; j < n-1; j++) {
			//go backwards following the velocity field
			backtraced_position = vec2(i,j);
			backtraced_position = backtraced_position + velocity[i+n*j]*(-h); 
			if (backtraced_position.x < 0.5)
				backtraced_position.x = 0.5;
			if (backtraced_position.x > n-1.5)
				backtraced_position.x = n-1.5;
			if (backtraced_position.y < 0.5)
				backtraced_position.y = 0.5;
			if (backtraced_position.y > n-1.5)
				backtraced_position.y = n-1.5;

			//bilinear interpolation
			int i0 = (int) backtraced_position.x;
			int j0 = (int) backtraced_position.y;

			double s = backtraced_position.x - i0;
			double t = backtraced_position.y - j0;
			density[i+j*n] = (1-s)*(1-t)* density_source[i0+j0*n] + (1-s)*t* density_source[i0+(j0+1)*n] + s*(1-t)* density_source[i0+1+j0*n] + s*t* density_source[i0+1+(j0+1)*n];
		}
}

void CFluidSolver::velocity_advection()
{
	//set boundary condition for source velocity field
	set_boundary_velocity(velocity);

	for (int i = 1; i < n-1; i++) {
		for (int j=1; j < n-1; j++) {
			double back_x = i - velocity[i + j*n].x * h;
			double back_y = j - velocity[i + j*n].y * h;

			if (back_x < 0.5) back_x = 0.5;
			if (back_x > n - 1.5) back_x = n - 1.5;
			if (back_y < 0.5) back_y = 0.5;
			if (back_y > n - 1.5) back_y = n - 1.5;

			int i0 = (int)back_x;
			int j0 = (int)back_y;
			double s = back_x - i0;
			double t = back_y - j0;

			vec2 v00 = velocity[i0 + j0*n];
			vec2 v01 = velocity[i0 + (j0+1)*n];
			vec2 v10 = velocity[i0+1 + j0*n];
			vec2 v11 = velocity[i0+1 + (j0+1)*n];

			double vx = (1-s)*(1-t)*v00.x + (1-s)*t*v01.x + s*(1-t)*v10.x + s*t*v11.x;
			double vy = (1-s)*(1-t)*v00.y + (1-s)*t*v01.y + s*(1-t)*v10.y + s*t*v11.y;
			advected_velocity[i+j*n] = vec2(vx, vy);
		}
	}

	set_boundary_velocity(advected_velocity);
}

void CFluidSolver::set_boundary_velocity(vec2* field)
{
	for (int i=0; i < n; i++) {
		field[0 + i*n] = vec2(0., 0.);
		field[(n-1) + i*n] = vec2(0., 0.);
		field[i + 0*n] = vec2(0., 0.);
		field[i + (n-1)*n] = vec2(0., 0.);
	}
}

void CFluidSolver::apply_buoyancy()
{
	for (int i = 1; i < n-1; i++) {
		for (int j = 1; j < n-1; j++) {
			int index = i + j*n;
			// Screen-space y grows downward, so negative y pushes smoke upward.
			velocity[index].y -= buoyancy_coef * density[index] * h;
		}
	}
}

void CFluidSolver::velocity_diffusion()
{
	double alpha = viscosity_coef * h;
	if (alpha <= 0.) {
		return;
	}

	CSparseMatrix velocity_diffusion_matrix(size, size);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int index = i + j*n;
			int count = 0;
			if (i-1>0) {
				velocity_diffusion_matrix.set1Value(index, index-1, -alpha);
				count++;
			}
			if (i+1<n-1) {
				velocity_diffusion_matrix.set1Value(index, index+1, -alpha);
				count++;
			}
			if (j-1>0) {
				velocity_diffusion_matrix.set1Value(index, index-n, -alpha);
				count++;
			}
			if (j+1<n-1) {
				velocity_diffusion_matrix.set1Value(index, index+n, -alpha);
				count++;
			}
			velocity_diffusion_matrix.set1Value(index, index, 1. + count * alpha);
		}
	}

	double* rhs_x = new double[size];
	double* rhs_y = new double[size];
	double* sol_x = new double[size];
	double* sol_y = new double[size];
	for (int i = 0; i < size; i++) {
		rhs_x[i] = velocity[i].x;
		rhs_y[i] = velocity[i].y;
		sol_x[i] = velocity[i].x;
		sol_y[i] = velocity[i].y;
	}

	velocity_diffusion_matrix.solve(sol_x, rhs_x, 1e-8, 30);
	velocity_diffusion_matrix.solve(sol_y, rhs_y, 1e-8, 30);

	for (int i = 0; i < size; i++) {
		velocity[i] = vec2(sol_x[i], sol_y[i]);
	}

	delete[] rhs_x;
	delete[] rhs_y;
	delete[] sol_x;
	delete[] sol_y;

	set_boundary_velocity(velocity);
}

void CFluidSolver::increase_viscosity(double step)
{
	viscosity_coef += step;
}

void CFluidSolver::decrease_viscosity(double step)
{
	viscosity_coef -= step;
	if (viscosity_coef < 0.) {
		viscosity_coef = 0.;
	}
}