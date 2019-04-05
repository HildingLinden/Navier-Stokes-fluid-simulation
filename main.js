"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");

	let N = 50;
	let SCALE = 8;

	// n, diffusion, viscosity, dt
	let grid = new FluidGrid(N, 0.00001, 0.00001, 0.1);

	let velX = 0;
	let velY = 1;

	function draw() {
		context.clearRect(0, 0, canvas.width, canvas.height);
		for (let i = 0; i < N; i++) {
			for (let j = 0; j < N; j++) {
				let color = 255 * grid.density[i+j*N];
				context.fillStyle = "rgba("+color+","+color+","+color+",255)";
				context.fillRect(i*SCALE, j*SCALE, SCALE, SCALE);
			}
		}
		grid.addVelocity(25,25,velX,velY);
		grid.addDensity(25,25,2);
		grid.stepForward();
		window.requestAnimationFrame(draw);
	}

	window.requestAnimationFrame(draw);

	canvas.addEventListener("mousemove", (event) => {
		//console.log((event.clientX-210)/(event.clientY-210));
		let x = event.clientX - 210;
		let y = event.clientY - 210;
		let mag = Math.sqrt(x*x+y*y);
		velX = x/mag;
		velY = y/mag;
	});

	function scale(x, in_min, in_max, out_min, out_max) {
		return(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
	}
	/*
	let imageData = context.createImageData(400, 400);
	let data = imageData.data;
	let steps = 0;
	function step() {
		context.clearRect(0, 0, canvas.width, canvas.height);
		context.putImageData(imageData, steps%400, 100);
		steps += 2;
		window.requestAnimationFrame(step);
	}

	window.requestAnimationFrame(step);
	*/
}

