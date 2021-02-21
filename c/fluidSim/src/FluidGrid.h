#pragma once

enum class Direction { NONE, HORIZONTAL, VERTICAL };

class FluidGrid {
public:
	~FluidGrid();
	FluidGrid(int size);

	void	step(float dt, int iterations, float diffusionRate, float viscosity, float fadeRate);
	void	addDensity(int x, int y, float amount, float dt);
	void    addVelocity(int x, int y, float velocityXAmount, float velocityYAmount, float dt);
	float *getDensity();

private:
	float *density, *prevDensity;
	float *velocityX, *prevVelocityX;
	float *velocityY, *prevVelocityY;
	float *tmp;

	int		size;

	void	advect(Direction direction, float *arr, float *prevArr, float *velocityX, float *velocityY, float dt);
	void	project(int iterations, float *velocityX, float *velocityY, float *p, float *div);
	void	diffuse(Direction direction, int iterations, float *arr, float *prevArr, float dt, float diffusion);
	void	linearSolve(Direction direction, int iterations, float *arr, float *prevArr, float neighborDiffusion, float scaling);
	void	setBounds(Direction direction, float *arr);
	void	fadeDensity(float dt, float fadeRate);
};