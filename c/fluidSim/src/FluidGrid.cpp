#include "FluidGrid.h"
#include <iostream>
#include <algorithm>
#include <chrono>
#include <immintrin.h>

FluidGrid::FluidGrid(int size) {
	printT0 = std::chrono::steady_clock::now();

	density = (float *)calloc(size * size, sizeof(float));
	prevDensity = (float *)calloc(size * size, sizeof(float));

	velocityX = (float *)calloc(size * size, sizeof(float));
	prevVelocityX = (float *)calloc(size * size, sizeof(float));

	velocityY = (float *)calloc(size * size, sizeof(float));
	prevVelocityY = (float *)calloc(size * size, sizeof(float));

	tmp = (float *)calloc(size * size, sizeof(float));

	this->size = size;

	if (!(density && prevDensity && velocityX && prevVelocityX && velocityY && prevVelocityY)) {
		throw std::runtime_error("Failed to allocate memory for FluidGrid\n");
	}
}

FluidGrid::~FluidGrid() {
	free(density);
	free(prevDensity);
	free(velocityX);
	free(prevVelocityX);
	free(velocityY);
	free(prevVelocityY);
	free(tmp);
}

void FluidGrid::step(float dt, int iterations, double diffusionRate, double viscosity, double fadeRate) {
	// Diffuse X velocity
	std::swap(this->velocityX, this->prevVelocityX);
	diffuse(Direction::HORIZONTAL, iterations, this->velocityX, this->prevVelocityX, dt, viscosity);

	// Diffuse Y velocity
	std::swap(this->velocityY, this->prevVelocityY);
	diffuse(Direction::VERTICAL, iterations, this->velocityY, this->prevVelocityY, dt, viscosity);

	// Conserve mass of the velocity field
	project(iterations, this->velocityX, this->velocityY, this->prevVelocityX, this->prevVelocityY);

	// Self advection
	std::swap(this->velocityX, this->prevVelocityX);
	std::swap(this->velocityY, this->prevVelocityY);
	startTimer();
	advect(Direction::HORIZONTAL, this->velocityX, this->prevVelocityX, this->prevVelocityX, this->prevVelocityY, dt);
	endTimer("Advect");

	startTimer();
	advect(Direction::VERTICAL, this->velocityY, this->prevVelocityY, this->prevVelocityX, this->prevVelocityY, dt);
	endTimer("Advect");

	// Conserve mass of the velocity field
	project(iterations, this->velocityX, this->velocityY, this->prevVelocityX, this->prevVelocityY);

	// Diffuse density
	std::swap(this->density, this->prevDensity);
	diffuse(Direction::NONE, iterations, this->density, this->prevDensity, dt, diffusionRate);

	// Advect density
	std::swap(this->density, this->prevDensity);
	startTimer();
	advect(Direction::NONE, this->density, this->prevDensity, this->velocityX, this->velocityY, dt);
	endTimer("Advect");

	// Fade the density to avoid filling the volume
	fadeDensity(dt, fadeRate);

	checkPrint();
}

void FluidGrid::addVelocity(int x, int y, float amountX, float amountY, float dt) {
	this->velocityX[x + y * this->size] += amountX * dt * this->size;
	this->velocityY[x + y * this->size] += amountY * dt * this->size;
}

void FluidGrid::addDensity(int x, int y, float amount, float dt) {
	this->density[x + y * this->size] += amount * dt * this->size;
}

// Linear backtrace
void FluidGrid::advect(Direction direction, float *arr, float *prevArr, float *velX, float *velY, float dt) {
	float scaling = dt * this->size;
	float maxClamp = this->size - 1.5f;

	__m256 _xIncrements = _mm256_set_ps(7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f);
	__m256 _scaling = _mm256_set1_ps(scaling);
	__m256 _maxClamp = _mm256_set1_ps(maxClamp);
	__m256 _minClamp = _mm256_set1_ps(0.5f);
	__m256 _one = _mm256_set1_ps(1.0f);
	__m256i _size = _mm256_set1_epi32(this->size);

	for (int y = 1; y < size-1; y++) {
		int x = 1;
		
		__m256 _y = _mm256_set1_ps(y);
		for (; x < size - 9; x += 8) {
			// Increment x by 0 to 7 since we are computing 8 iterations of x at once
			__m256 _x = _mm256_set1_ps(x);
			_x = _mm256_add_ps(_x, _xIncrements);

			// Subtract the current scaled velocity from y and incremented x
			__m256 _prevX = _mm256_load_ps(&velX[x + y * size]);
			__m256 _prevY = _mm256_load_ps(&velY[x + y * size]);

			_prevX = _mm256_mul_ps(_prevX, _scaling);
			_prevY = _mm256_mul_ps(_prevY, _scaling);
			_prevX = _mm256_sub_ps(_x, _prevX);
			_prevY = _mm256_sub_ps(_y, _prevY);

			// Clamp indices to minClamp and maxClamp
			_prevX = _mm256_max_ps(_prevX, _minClamp);
			_prevX = _mm256_min_ps(_prevX, _maxClamp);
			_prevY = _mm256_max_ps(_prevY, _minClamp);
			_prevY = _mm256_min_ps(_prevY, _maxClamp);

			// Get the integer part and decimal part of the indices
			// as well as 1 - decimal part
			__m256 _prevXRound = _mm256_floor_ps(_prevX);
			__m256 _prevYRound = _mm256_floor_ps(_prevY);

			__m256i _prevXInt = _mm256_cvtps_epi32(_prevXRound);
			__m256i _prevYInt = _mm256_cvtps_epi32(_prevYRound);

			__m256 _prevXDecimals = _mm256_sub_ps(_prevX, _prevXRound);
			__m256 _prevYDecimals = _mm256_sub_ps(_prevY, _prevYRound);

			__m256 _prevXRemainingDecimals = _mm256_sub_ps(_one, _prevXDecimals);
			__m256 _prevYRemainingDecimals = _mm256_sub_ps(_one, _prevYDecimals);

			// Multiply y by size and add x to get the actual indices
			__m256i _prevIndices = _mm256_mullo_epi32(_prevYInt, _size);
			_prevIndices = _mm256_add_epi32(_prevIndices, _prevXInt);

			// Get the middle, right, down and downRight cells by changing the base address of the gather
			__m256 _middle = _mm256_i32gather_ps(&prevArr[0], _prevIndices, sizeof(float));
			__m256 _right = _mm256_i32gather_ps(&prevArr[1], _prevIndices, sizeof(float));
			__m256 _down = _mm256_i32gather_ps(&prevArr[this->size], _prevIndices, sizeof(float));
			__m256 _downRight = _mm256_i32gather_ps(&prevArr[this->size + 1], _prevIndices, sizeof(float));

			// Do the linear interpolation (See below for explanation)
			__m256 _middleScaling = _mm256_mul_ps(_prevXRemainingDecimals, _prevYRemainingDecimals);
			__m256 _downScaling = _mm256_mul_ps(_prevXRemainingDecimals, _prevYDecimals);
			__m256 _rightScaling = _mm256_mul_ps(_prevXDecimals, _prevYRemainingDecimals);
			__m256 _downRightScaling = _mm256_mul_ps(_prevXDecimals, _prevYDecimals);

			__m256 _result = _mm256_mul_ps(_middle, _middleScaling);
			_result = _mm256_fmadd_ps(_down, _downScaling, _result);
			_result = _mm256_fmadd_ps(_right, _rightScaling, _result);
			_result = _mm256_fmadd_ps(_downRight, _downRightScaling, _result);

			_mm256_storeu_ps(&arr[x + y * this->size], _result);
		}

		for ( ; x < size-1; x++) {
			// Subtract velocity from current cell to get previous cell indices
			float prevX = x - velX[x+y*size] * scaling;
			float prevY = y - velY[x+y*size] * scaling;

			// Clamp previous cells in case they are outside of the domain
			prevX = std::clamp(prevX, 0.5f, maxClamp);
			prevY = std::clamp(prevY, 0.5f, maxClamp);

			/*
			Casting the index of the previous cell to integers truncates the index which means that
			the index is moved up and to the left
			*/
			int prevXInt = (int)prevX;
			int prevYInt = (int)prevY;

			/*
			To interpolate we then:
				multiply the decimals of x and y with the cell below and to the right
				multiply the decimals of x and remaining decimals of y with the cell to the right
				multiply the decimals of y and remaining decimals of x with the cell above
				multiply the remaining decimals of x and y with the current cell
			and add them all together
			*/
			float prevXDecimals = prevX - prevXInt;
			float prevXRemainingDecimals = 1.0f - prevXDecimals;

			float prevYDecimals = prevY - prevYInt;
			float prevYRemainingDecimals = 1.0f - prevYDecimals;
			
			arr[x + y * this->size] =
				prevXRemainingDecimals * (prevYRemainingDecimals * prevArr[prevXInt	    + prevYInt * this->size] + prevYDecimals * prevArr[prevXInt	   + (prevYInt+1) * this->size]) +
				prevXDecimals		   * (prevYRemainingDecimals * prevArr[(prevXInt+1) + prevYInt * this->size] + prevYDecimals * prevArr[(prevXInt+1) + (prevYInt+1) * this->size]);
		}
	}

	setBounds(direction, arr);
}

// Hodge decomposition
void FluidGrid::project(int iterations, float *velocityX, float *velocityY, float *p, float *div) {
	float sizeReciprocal = 1.0f / this->size;

	__m256 _sizeReciprocal = _mm256_set1_ps(sizeReciprocal);
	__m256 _minusHalf = _mm256_set1_ps(-0.5f);
	__m256 _zero = _mm256_setzero_ps();

	// Compute height field (Poisson-pressure equation)
	for (int y = 1; y < this->size - 1; y++) {
		int x = 1;

		for ( ; x < this->size - 9; x+=8) {
			/*
			x = x[left] - x[right]
			y = y[up] - y[down]
			cell = x + y
			cell *= -0.5
			cell *= sizeReciprocal
			*/
			__m256 _up = _mm256_loadu_ps(&velocityY[x + (y - 1) * this->size]);
			__m256 _down = _mm256_loadu_ps(&velocityY[x + (y + 1) * this->size]);
			__m256 _left = _mm256_loadu_ps(&velocityX[(x - 1) + y * this->size]);
			__m256 _right = _mm256_loadu_ps(&velocityX[(x + 1) + y * this->size]);

			__m256 _x = _mm256_sub_ps(_left, _right);
			__m256 _y = _mm256_sub_ps(_up, _down);
			__m256 _cell = _mm256_add_ps(_x, _y);
			_cell = _mm256_mul_ps(_cell, _minusHalf);
			_cell = _mm256_mul_ps(_cell, _sizeReciprocal);
			_mm256_storeu_ps(&div[x + y * size], _cell);
			_mm256_storeu_ps(&p[x + y * size], _zero);
		}

		// Take care of rest of the cells if x-2 is not evenly divisible by 8
		for (; x < this->size - 1; x++) {
			div[x + y * size] =
				-0.5f * (
					velocityX[(x - 1) + y * this->size] - velocityX[(x + 1) + y * this->size] +
					velocityY[x + (y - 1) * this->size] - velocityY[x + (y + 1) * this->size]
					) * sizeReciprocal;
			p[x + y * size] = 0;
		}
	}

	setBounds(Direction::NONE, div);
	setBounds(Direction::NONE, p);
	linearSolve(Direction::NONE, iterations, p, div, 1, 4);

	// Compute mass conserving field (Velocity field - Height field)
	__m256 _size = _mm256_set1_ps(this->size);

	for (int y = 1; y < this->size - 1; y++) {
		int x = 1;

		for ( ; x < this->size - 9; x+=8) {
			__m256 _val, _oldVal;
			/*
			val = p[left] - p[right];
			val *= size;
			oldVal += -0.5f * val;
			*/
			__m256 _left  = _mm256_loadu_ps(&p[(x - 1) + y * this->size]);
			__m256 _right = _mm256_loadu_ps(&p[(x + 1) + y * this->size]);			
			_val = _mm256_sub_ps(_left, _right);
			_val = _mm256_mul_ps(_val, _size);
			_oldVal = _mm256_loadu_ps(&velocityX[x + y * this->size]);
			_oldVal = _mm256_fmadd_ps(_val, _minusHalf, _oldVal);
			_mm256_storeu_ps(&velocityX[x + y * this->size], _oldVal);

			/*
			val = p[up] - p[down];
			val *= size;
			oldVal += -0.5f * val;
			*/
			__m256 _up = _mm256_loadu_ps(&p[x + (y - 1) * this->size]);
			__m256 _down = _mm256_loadu_ps(&p[x + (y + 1) * this->size]);
			_val = _mm256_sub_ps(_up, _down);
			_val = _mm256_mul_ps(_val, _size);
			_oldVal = _mm256_loadu_ps(&velocityY[x + y * this->size]);
			_oldVal = _mm256_fmadd_ps(_val, _minusHalf, _oldVal);
			_mm256_storeu_ps(&velocityY[x + y * this->size], _oldVal);
		}

		// Take care of rest of the cells if x-2 is not evenly divisible by 8
		for (; x < this->size - 1; x++) {
			velocityX[x + y * this->size] -= 0.5f * (p[(x - 1) + y * this->size] - p[(x + 1) + y * this->size]) * this->size;
			velocityY[x + y * this->size] -= 0.5f * (p[x + (y - 1) * this->size] - p[x + (y + 1) * this->size]) * this->size;
		}
	}

	setBounds(Direction::HORIZONTAL, velocityX);
	setBounds(Direction::VERTICAL, velocityY);
}

// Mass conserving
void FluidGrid::diffuse(Direction direction, int iterations, float *arr, float *prevArr, float dt, double diffusion) {
	float neighborDiffusion = diffusion * dt * this->size * this->size;
	float scaling = 1 + 4 * neighborDiffusion;

	linearSolve(direction, iterations, arr, prevArr, neighborDiffusion, scaling);
}

// Using Jacobi relaxation
void FluidGrid::linearSolve(Direction direction, int iterations, float *arr, float *prevArr, float neighborDiffusion, float scaling) {
	float reciprocalScaling = 1.0f / scaling;
	__m256 _reciprocalScaling = _mm256_set1_ps(reciprocalScaling);
	__m256 _neighborDiffusion = _mm256_set1_ps(neighborDiffusion);

	for (int iteration = 0; iteration < iterations; iteration++) {
		for (int y = 1; y < this->size - 1; y++) {
			int x = 1;

			for ( ; x < this->size - 9; x+= 8) {
				__m256 _up  = _mm256_loadu_ps(&arr[x		  + (y - 1) * this->size]);
				__m256 _left  = _mm256_loadu_ps(&arr[(x - 1)  + y		* this->size]);
				__m256 _right = _mm256_loadu_ps(&arr[(x + 1)  + y		* this->size]);
				__m256 _down  = _mm256_loadu_ps(&arr[x		  + (y + 1) * this->size]);

				/*
				cell = previous;
				cell +=	(up * neighborDiffusion);
				cell += (left * neighborDiffusion);
				cell += (right * neighborDiffusion); 
				cell += (down * neighborDiffusion);
				cell *= reciprocalScaling				
				*/
				__m256 _cell = _mm256_loadu_ps(&prevArr[x + y * this->size]);
				_cell = _mm256_fmadd_ps(_up, _neighborDiffusion, _cell);
				_cell = _mm256_fmadd_ps(_left, _neighborDiffusion, _cell);
				_cell = _mm256_fmadd_ps(_right, _neighborDiffusion, _cell);
				_cell = _mm256_fmadd_ps(_down, _neighborDiffusion, _cell);			
				_cell = _mm256_mul_ps(_cell, _reciprocalScaling);

				// Store result
				_mm256_storeu_ps(&tmp[x + y * this->size], _cell);
			}

			// Take care of rest of the cells if x-2 is not evenly divisible by 8
			for ( ; x < this->size - 1; x++) {
				float neighbors =
					neighborDiffusion * (
						arr[x + ((y - 1) * this->size)] +
						arr[x + ((y + 1) * this->size)] +
						arr[(x - 1) + (y * this->size)] +
						arr[(x + 1) + (y * this->size)]
						);
				float previous = prevArr[x + y * this->size];
				tmp[x + y * this->size] = (previous + neighbors) * reciprocalScaling;
			}
		}
		std::swap(tmp, arr);

		// Controlling boundary after every iterations
		setBounds(direction, arr);
	}
}

void FluidGrid::setBounds(Direction direction, float *arr) {
	for (int i = 1; i < this->size-1; i++) {
		// Left and right edge get the reverse velocity of the neighbor if in X direction
		arr[0		 + i * size] = (direction == Direction::HORIZONTAL) ? -arr[1		  + i * size] : arr[1		 + i * size];
		arr[(size-1) + i * size] = (direction == Direction::HORIZONTAL) ? -arr[(size-2) + i * size] : arr[(size-2) + i * size];

		// Top and bottom edge get the reverse velocity of the neighbor if in Y direction
		arr[i + 0		 * size] = (direction == Direction::VERTICAL) ? -arr[i + 1		* size] : arr[i + 1		   * size];
		arr[i + (size-1) * size] = (direction == Direction::VERTICAL) ? -arr[i + (size-2) * size] : arr[i + (size-2) * size];
	}

	// The corners get the average values of their neighbors
	arr[0		 + 0		* size] = 0.5f * (arr[1		  + 0		 * size] + arr[0		+ 1		   * size]);
	arr[(size-1) + 0		* size] = 0.5f * (arr[(size-2) + 0		 * size] + arr[(size-1) + 1		   * size]);
	arr[0		 + (size-1)	* size] = 0.5f * (arr[0		  + (size-2) * size] + arr[1		+ (size-1) * size]);
	arr[(size-1) + (size-1) * size] = 0.5f * (arr[(size-2) + (size-1) * size] + arr[(size-1) + (size-2) * size]);
}

void FluidGrid::fadeDensity(float dt, double fadeRate) {
	double scaledFadeRate = 1 - dt * fadeRate * this->size;
	for (int y = 1; y < this->size - 1; y++) {
		for (int x = 1; x < this->size - 1; x++) {
			this->density[x + y * this->size] *= scaledFadeRate;
		}
	}
}

void FluidGrid::startTimer() {
	timerT0 = std::chrono::steady_clock::now();
}

void FluidGrid::endTimer(std::string timerName) {
	auto t1 = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::micro> microSeconds = t1 - timerT0;
	
	std::map<std::string, double>::iterator iterator = timers.find(timerName);
	if (iterator != timers.end()) {
		iterator->second += microSeconds.count();
	}
	else {
		timers.insert(std::pair<std::string, double>(timerName, microSeconds.count()));
	}
}

void FluidGrid::checkPrint() {
	runs++;
	auto printT1 = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::micro> printMicroSeconds = printT1 - printT0;
	microsSinceLastPrint += printMicroSeconds.count();
	if (microsSinceLastPrint > 1000000.0) {
		for (std::map<std::string, double>::iterator iterator = timers.begin(); iterator != timers.end(); iterator++) {
			std::cout << iterator->first << ": " << iterator->second / runs << " microseconds" << std::endl;
			iterator->second = 0;
		}
		microsSinceLastPrint = 0;
		runs = 0;
	}

	printT0 = printT1;
}

float *FluidGrid::getDensity() {
	return this->density;
}