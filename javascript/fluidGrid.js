"use strict";

class FluidGrid {
	constructor(size, diffusion, viscosity, fade, dt, iterations) {
		this.dt = dt;
		this.iterations = iterations;

		this.size = size;
		this.diffusionRate = diffusion;
		this.viscosity = viscosity;
		this.fadeRate = fade;

		this.density = new Array(size*size).fill(0);
		this.prevDensity = new Array(size*size).fill(0);

		this.velocityX = new Array(size*size).fill(0);
		this.velocityY = new Array(size*size).fill(0);
		this.prevVelocityX = new Array(size*size).fill(0);
		this.prevVelocityY = new Array(size*size).fill(0);
	}

	// Add density (dye) to a cell
	addDensity(x, y, amount) {
		this.density[x + y * this.size] += amount;
	}

	// Add velocity (in x & y direction) to a cell
	addVelocity(x, y, amountX, amountY) {
		this.velocityX[x+y*this.size] += amountX;
		this.velocityY[x+y*this.size] += amountY;
	}

	step() {
		let t0 = performance.now();

		// Diffuse the velocities of the fluid (viscous diffusion)
		diffuse(1, this.prevVelocityX, this.velocityX, this.viscosity, this.dt, this.iterations, this.size);
		diffuse(2, this.prevVelocityY, this.velocityY, this.viscosity, this.dt, this.iterations, this.size);

		// Equalize cells
		project(this.prevVelocityX, this.prevVelocityY, this.velocityX, this.velocityY, this.iterations, this.size);

		// Move the velocity of the fluid (self-advection), swapping order of arr and prev array
		advect(1, this.velocityX, this.prevVelocityX, this.prevVelocityX, this.prevVelocityY, this.dt, this.size);
		advect(2, this.velocityY, this.prevVelocityY, this.prevVelocityX, this.prevVelocityY, this.dt, this.size);

		// Equalize cells
		project(this.velocityX, this.velocityY, this.prevVelocityX, this.prevVelocityY, this.iterations, this.size);

		// Diffuse and move the density (dye), swapping order of arr and prev array
		diffuse(0, this.prevDensity, this.density, this.diffusionRate, this.dt, this.iterations, this.size);
		advect(0, this.density, this.prevDensity, this.velocityX, this.velocityY, this.dt, this.size);

		// Fade the dye
		this.density = this.density.map(x => x * (1-this.fadeRate));

		let t1 = performance.now();
		return t1-t0;
	}
}

function diffuse(direction, arr, prevArr, diffusionRate, dt, iterations, size) {
	let diffusion = dt * diffusionRate * size * size;

	// Guass-Seidel relaxation
	for (let iteration = 0; iteration < iterations; iteration++) {
		for (let y = 1; y < size - 1; y++) {
			for (let x = 1; x < size - 1; x++) {

				let up 	  =	arr[x		+ (y-1) * size];
				let down  =	arr[x 		+ (y+1) * size];
				let left  =	arr[(x-1)	+ y 	* size];
				let right = arr[(x+1)	+ y 	* size];

				arr[x+y*size] = (prevArr[x+y*size] + diffusion * (up + down + left + right)) / (1 + 4 * diffusion);
			}
		}
		// Set the edges of the whole fluid
		set_bounds(direction, arr, size);
	}
}

function project(velocityX, velocityY, p, div, iterations, size) {
	// Hodge decomposition, poisson equation

	let h = 1/size;

	for (let y = 1; y < size-1; y++) {
		for (let x = 1; x < size-1; x++) {
			div[x+y*size] = -0.5*h*(
				velocityX[(x+1) + y 	* size] - velocityX[(x-1) + y 	  * size] +
				velocityY[x     + (y+1) * size] - velocityY[x 	  + (y-1) * size]
			);
			p[x+y*size] = 0;
		}
	}

	set_bounds(0, div, size);
	set_bounds(0, p, size);

	// Gauss-Seidel relaxation
	// Conjugate gradient solver might improve quality
	for (let iteration = 0; iteration < iterations; iteration++) {
		for (let y = 1; y < size-1; y++) {
			for (let x = 1; x < size-1; x++) {
				let up 	  =	p[x		+ (y-1) * size];
				let down  =	p[x 		+ (y+1) * size];
				let left  =	p[(x-1)	+ y 	* size];
				let right = p[(x+1)	+ y 	* size];

				p[x+y*size] = (div[x+y*size] + up + down + left + right) / 4;
			}
		}
		set_bounds(0, p, size);
	}

	for (let y = 1; y < size-1; y++) {
		for (let x = 1; x < size-1; x++) {
			velocityX[x+y*size] -= 0.5*(p[(x+1) + y 	* size] - p[(x-1) + y 	  * size])/h;
			velocityY[x+y*size] -= 0.5*(p[x     + (y+1) * size] - p[x 	  + (y-1) * size])/h;
		}
	}

	set_bounds(1, velocityX, size);
	set_bounds(2, velocityY, size);
}

function advect(direction, arr, prevArr, velocityX, velocityY, dt, size) {
	let dt0 = dt * size;

	for (let y = 1; y < size - 1; y++) {
		for (let x = 1; x < size - 1; x++) {
			// Get the previous x and y location of current cell
			let prevX = x - dt0 * velocityX[x+y*size];
			let prevY = y - dt0 * velocityY[x+y*size];

			// Clamp to 0.5 from edges
			prevX = Math.min(Math.max(0.5, prevX), size-1+0.5);
			prevY = Math.min(Math.max(0.5, prevY), size-1+0.5);

			let prevXInt = Math.floor(prevX);
			let prevYInt = Math.floor(prevY);

			let s1 = prevX - prevXInt;
			let s0 = 1 - s1;
			let t1 = prevY - prevYInt;
			let t0 = 1 - t1;

			arr[x+y*size] =
				s0 * (t0 * prevArr[prevXInt     + prevYInt * size] + t1 * prevArr[prevXInt     + (prevYInt+1) * size]) +
				s1 * (t0 * prevArr[(prevXInt+1) + prevYInt * size] + t1 * prevArr[(prevXInt+1) + (prevYInt+1) * size]);
		}
	}

	// Set the edges of the whole fluid
	set_bounds(direction, arr, size);
}

function set_bounds(direction, arr, size) {
	for (let i = 1; i < size; i++) {
		// Left and right edge get the reverse velocity of the neighbor if in X direction
		arr[0		 + i * size] = direction == 1 ? -arr[1 	  	  + i * size] : arr[1		 + i * size];
		arr[(size-1) + i * size] = direction == 1 ? -arr[(size-2) + i * size] : arr[(size-2) + i * size];

		// Top and bottom edge get the reverse velocity of the neighbor if in Y direction
		arr[i + 0		 * size] = direction == 2 ? -arr[i + 1 	  	  * size] : arr[i + 1 	     * size];
		arr[i + (size-1) * size] = direction == 2 ? -arr[i + (size-2) * size] : arr[i + (size-2) * size];
	}

	// The corners get the average values of their neighbors
	arr[0        + 0 	  	* size] = 0.5 * (arr[1 		  + 0 	 	 * size] + arr[0 		+ 1 	   * size]);
	arr[(size-1) + 0 	  	* size] = 0.5 * (arr[(size-2) + 0 	 	 * size] + arr[(size-1) + 1		   * size]);
	arr[0 	     + (size-1) * size] = 0.5 * (arr[0 		  + (size-2) * size] + arr[1 		+ (size-1) * size]);
	arr[(size-1) + (size-1) * size] = 0.5 * (arr[(size-2) + (size-1) * size] + arr[(size-1) + (size-2) * size]);
}