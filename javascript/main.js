"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");
	context.imageSmoothingEnabled = false;

	let N = 100;
	let SCALE = canvas.width/N;

	// n, diffusion, viscosity, fadeRate, dt, iterations
	let grid = new FluidGrid(N, 0.001, 0.00001, 0.01, 0.016, 4);

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
		grid.addVelocity(N/2, N/2, velX, velY);
		grid.addDensity(N/2, N/2, 2);
		grid.step();

		let t0 = performance.now();

		context.clearRect(0, 0, canvas.width, canvas.height);

		// Fill all pixels of frame
		for (let y = 0; y < N; y++) {
			for (let x = 0; x < N; x++) {
				let color = 255 * grid.density[x+y*N];
				frame.data[x*4 + y*4*N] = color;
				frame.data[x*4 + y*4*N + 1] = color;
				frame.data[x*4 + y*4*N + 2] = color;
			}
		}

		context.putImageData(frame, 0, 0);
		context.drawImage(canvas, 0, 0, N, N, 0, 0, canvas.width, canvas.height);

		let t1 = performance.now();
		//console.log("Draw ", t1-t0);

		window.requestAnimationFrame(loop);
	}

	window.requestAnimationFrame(loop);

	canvas.addEventListener("mousemove", (event) => {
		let x = event.clientX - canvas.width/2 - 10;
		let y = event.clientY - canvas.height/2 - 10;
		let angle = Math.atan2(y, x);
		velX = 4 * Math.cos(angle);
		velY = 4 * Math.sin(angle);
	});

	function scale(x, in_min, in_max, out_min, out_max) {
		return(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
	}
}

