"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");
	context.imageSmoothingEnabled = false;

	let stepTimeText = document.getElementById("stepTime");
	let drawTimeText = document.getElementById("drawTime");

	let N = 200;
	let SCALE = canvas.width/N;

	// n, diffusion, viscosity, fadeRate, dt, iterations
	let grid = new FluidGrid(N, 0.001, 0.00001, 0.0001, 0.016, 4);

	let velX = 0;
	let velY = 0;

	// Create inital frame and set alpha to 255
	let frame = context.createImageData(N, N);
	for (let y = 0; y < N; y++) {
		for (let x = 0; x < N; x++) {
			frame.data[x*4 + y*4*N + 3] = 255;
		}
	}

	function loop() {
		// Add velocity to the middle of the canvas in the direction of the mouse
		grid.addVelocity(N/2, N/2, velX, velY);

		// Add density to the middle of the canvas
		grid.addDensity(N/2, N/2, 5);

		// Step the physics forward a time step
		let stepTime = grid.step();

		stepTimeText.innerHTML = "Physics step time " + stepTime.toFixed(3) + "ms";

		let t0 = performance.now();

		// Reset canvas
		context.clearRect(0, 0, canvas.width, canvas.height);

		// Convert density into grayscale pixels in a frame
		for (let y = 0; y < N; y++) {
			for (let x = 0; x < N; x++) {
				let color = 255 * grid.density[x+y*N];
				frame.data[x*4 + y*4*N] = color;
				frame.data[x*4 + y*4*N + 1] = color;
				frame.data[x*4 + y*4*N + 2] = color;
			}
		}

		// Render frame on canvas then scale it up to the correct size
		context.putImageData(frame, 0, 0);
		context.drawImage(canvas, 0, 0, N, N, 0, 0, canvas.width, canvas.height);

		let t1 = performance.now();
		drawTimeText.innerHTML = "Draw time " + (t1-t0).toFixed(3) + "ms";

		window.requestAnimationFrame(loop);
	}

	window.requestAnimationFrame(loop);

	canvas.addEventListener("mousemove", (event) => {
		let x = event.clientX - canvas.width/2 - 10;
		let y = event.clientY - canvas.height/2 - 10;
		let angle = Math.atan2(y, x);
		velX = 8 * Math.cos(angle);
		velY = 8 * Math.sin(angle);
	});

	function scale(x, in_min, in_max, out_min, out_max) {
		return(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
	}
}

