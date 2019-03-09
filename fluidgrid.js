class FluidGrid {
	constructor(N, diffusion, viscosity, dt) {
		this.N = N;
		this.dt = dt;
		this.diff = diffusion;
		this.visc = viscosity;

		this.s = initializeArray(N*N);
		this.density = initializeArray(N*N);

		this.vectorX = initializeArray(N*N);
		this.vectorY = initializeArray(N*N);

		this.vectorX0 = initializeArray(N*N);
		this.vectorY0 = initializeArray(N*N);
	}

	// Add density (dye) to a point
	addDensity(x, y, amount) {
		let index = IX(x, y, this.N);
		this.density[index] += amount;
	}

	// Add velocity (in x & y direction) to a point
	addVelocity(x, y, amountX, amountY) {
		let index = IX(x, y, this.N);

		this.vectorX[index] += amountX;
		this.vectorY[index] += amountY;
	}

	stepForward() {
		let N = this.N;
		let dt = this.dt;
		let diff = this.diff;
		let visc = this.visc;		
		let s = this.s;
		let density = this.density;
		let vectorX = this.vectorX;
		let vectorY = this.vectorY;
		let vectorX0 = this.vectorX0;
		let vectorY0 = this.vectorY0;

		// Diffuse the velocity
		diffuse(1, vectorX, vectorX0, visc, dt, 4, N);
		diffuse(2, vectorY, vectorY0, visc, dt, 4, N);

		// Mirror at the edges
		project(vectorX0, vectorY0, vectorX, vectorY, 4, N);

		// Move velocity vectors
		advect(1, vectorX, vectorX0, vectorX0, vectorY0, dt, N);
		advect(2, vectorY, vectorY0, vectorX0, vectorY0, dt, N);

		// Mirror at the edges
		project(vectorX, vectorY, vectorX0, vectorY0, 4, N);

		// Diffuse and move the dye
		diffuse(0, s, density, diff, dt, 4, N);
		advect(0, density, s, vectorX, vectorY, dt, N);
	}
}

function IX(x, y, width) {
	return x + y * width;
}

function initializeArray(elements) {
	let arr = [];
	for (let i = 0; i < elements; i++) {
		arr[i] = 0;
	}

	return arr;
}

function set_bounds(b, x, N) {
	for (let i = 1; i < N - 1; i++) {
		x[IX(i, 0, N  )] = b == 2 ? -x[IX(i, 1, N  )] : x[IX(i, 1, N  )];
		x[IX(i, N-1, N)] = b == 2 ? -x[IX(i, N-2, N)] : x[IX(i, N-2, N)];
	}
	for (let j = 1; j < N - 1; j++) {
		x[IX(0,   j, N)] = b == 2 ? -x[IX(1,   j, N)] : x[IX(1,   j, N)];
		x[IX(N-1, j, N)] = b == 2 ? -x[IX(N-2, j, N)] : x[IX(N-2, j, N)];
	}

	x[IX(0, 0, N)] = 0.5 * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
	x[IX(0, N-1, N)] = 0.5 * (x[IX(1, N-1, N)] + x[IX(0, N-2, N)]);
	x[IX(N-1, 0, N)] = 0.5 * (x[IX(N-2, 0, N)] + x[IX(N-1, 1, N)]);
	x[IX(N-1, N-1, N)] = 0.5 * (x[IX(N-2, N-1, N)] + x[IX(N-1, N-2, N)]);
}

function diffuse(b, x, x0, diff, dt, iterations, N) {
	diffusionRate = dt * diff * (N - 2) * (N - 2);
	lin_solve(b, x, x0, diffusionRate, 1 + 6 * diffusionRate, iterations, N) 

}

function lin_solve(b, x, x0, diffusionRate, c, iterations, N) {
	cReciprocal = 1.0/c;
	for (let iteration = 0; iteration < iterations; iteration++) {
		for (let j = 1; j < N - 1; j++) {
			for (let i = 1; i < N - 1; i++) {
				let oldValue = x0[IX(i, j, N)];
				let up = x[IX(i, j-1, N)];
				let down = x[IX(i, j+1, N)];
				let left = x[IX(i-1, j, N)];
				let right = x[IX(i+1, j, N)];

				x[IX(i, j, N)] = 
					(oldValue + diffusionRate * (up + down + left + right)) * cReciprocal;
			}
		}
		set_bounds(b, x, N);
	}
}

function project(velocityX, velocityY, p, div, iterations, N) {
	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			div[IX(i, j, N)] = -0.5*(
				 velocityX[IX(i+1, j, N)]
				-velocityX[IX(i-1, j, N)]
				+velocityY[IX(i, j+1, N)]
				-velocityY[IX(i, j-1, N)]
				)/N;
			p[IX(i, j, N)] = 0;
		}
	}

	set_bounds(0, div, N);
	set_bounds(0, p, N);
	lin_solve(0, p, div, 1, 6, iterations, N);

	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			velocityX[IX(i, j, N)] -= 0.5 * (p[IX(i+1, j, N)] - p[IX(i-1, j, N)]) * N;
			velocityY[IX(i, j, N)] -= 0.5 * (p[IX(i, j+1, N)] - p[IX(i, j-1, N)]) * N;
		}
	}

	set_bounds(1, velocityX, N);
	set_bounds(2, velocityY, N);
}

function advect(b, d, d0, velocityX, velocityY, dt, N) {
	let dtx = dt * (N - 2);
	let dty = dt * (N - 2);

	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			let tmp1 = dtx * velocityX[IX(i, j, N)];
			let tmp2 = dty * velocityY[IX(i, j, N)];
			x = i - tmp1;
			y = j - tmp2;

			x = Math.max(0.5, x);
			x = Math.min(N+0.5, x);
			i0 = Math.floor(x);
			i1 = i0 + 1;
			y = Math.max(0.5, y);
			y = Math.min(N+0.5, y);
			j0 = Math.floor(y);
			j1 = j0 + 1;

			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			// Not sure if right
			d[IX(i, j, N)] =
				 s0 * (t0 * d0[IX(i0, j0, N)] + t1 * d0[IX(i0, j1, N)])
				+s1 * (t0 * d0[IX(i1, j0, N)] + t1 * d0[IX(i1, j1, N)]);
		}
	}

	set_bounds(b, d, N);
}