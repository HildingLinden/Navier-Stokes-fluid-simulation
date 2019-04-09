"use strict";

let vertexShaderCode = `#version 300 es
in vec2 a_position;
in vec2 a_textureCoord;

out vec2 v_textureCoord;

uniform vec2 u_resolution;

void main() {
	vec2 clipSpace = ((a_position/u_resolution) * 2.0) - 1.0;
	gl_Position = vec4(clipSpace * vec2(1, -1), 0, 1);

	v_textureCoord = a_textureCoord;
}`;

let fragmentShaderCode = `#version 300 es
precision mediump float;

in vec2 v_textureCoord;

out vec4 color;

uniform sampler2D u_texture;

void main() {
	color = texture(u_texture, v_textureCoord);
}
`;

function loadImage() {
	let img = new Image();
	img.addEventListener("load", () => {
		console.log("Image loaded");
		main(img);
	}, false);
	img.src = "test.png";
}

function main(img) {
	// Canvas and context setup
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("webgl2",
		{antialias: false, powerPreference: "high-performance"});

	// Enable 32-bit floating-point values for the color buffer
	context.getExtension("EXT_color_buffer_float");

	context.clearColor(0, 0, 0, 0);

	// Shader and program setup
	let vertexShader = createShader(context, context.VERTEX_SHADER, vertexShaderCode);
	let fragmentShader = createShader(context, context.FRAGMENT_SHADER, fragmentShaderCode);

	let program = linkShaders(context, vertexShader, fragmentShader);
	context.useProgram(program);

	// Attribute/Uniform setup
	let vao = context.createVertexArray();
	context.bindVertexArray(vao);

	let positionAttributeLocation = context.getAttribLocation(program, "a_position");
	let positionBuffer = context.createBuffer();
	context.bindBuffer(context.ARRAY_BUFFER, positionBuffer);
	let positions = new Float32Array([
		 1, 1,
		 1,99,
		99,99,
		99,99,
		99, 1,
		 1, 1
	]);
	context.bufferData(context.ARRAY_BUFFER, positions, context.STATIC_DRAW);
	context.enableVertexAttribArray(positionAttributeLocation);
	context.vertexAttribPointer(positionAttributeLocation, 2, context.FLOAT, false, 0, 0);

	let textureCoordAttributeLocation = context.getAttribLocation(program, "a_textureCoord");
	let textureCoordBuffer = context.createBuffer();
	context.bindBuffer(context.ARRAY_BUFFER, textureCoordBuffer);
	let textureCoord = new Float32Array([
		0, 0,
		0, 1,
		1, 1,
		1, 1,
		1, 0,
		0, 0
	]);
	context.bufferData(context.ARRAY_BUFFER, textureCoord, context.STATIC_DRAW);
	context.enableVertexAttribArray(textureCoordAttributeLocation);
	context.vertexAttribPointer(textureCoordAttributeLocation, 2, context.FLOAT, false, 0, 0);

	let resolutionUniformLocation = context.getUniformLocation(program, "u_resolution");
	context.uniform2f(resolutionUniformLocation, context.canvas.width, context.canvas.height);

	// Texture setup
	let texture = context.createTexture();
	context.activeTexture(context.TEXTURE0);
	context.bindTexture(context.TEXTURE_2D, texture);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_WRAP_S, context.CLAMP_TO_EDGE);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_WRAP_T, context.CLAMP_TO_EDGE);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_MIN_FILTER, context.NEAREST);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_MAG_FILTER, context.NEAREST);
	// context.texImage2D(context.TEXTURE_2D, 0, context.RGBA, context.RGBA, context.UNSIGNED_BYTE, img);
	// let textureUniformLocation = context.getUniformLocation(program, "u_texture");
	// context.uniform1i(textureUniformLocation, 0);

	// For compute
	// context.texImage2D(context.TEXTURE_2D, 0, context.R32F, 1, 1, 0, context.RED, context.FLOAT,
	// 	new Float32Array([1.234])
	// );
	// let fb = context.createFramebuffer();
	// context.bindFramebuffer(context.FRAMEBUFFER, fb);
	// context.framebufferTexture2D(context.FRAMEBUFFER, context.COLOR_ATTACHMENT0, context.TEXTURE_2D, texture, 0);

	// Render
	context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);
	context.drawArrays(context.TRIANGLES, 0, 6);

	// Read the drawn pixels for debugging
	let results = new Float32Array(context.drawingBufferWidth * context.drawingBufferHeight * 4);
	context.readPixels(0, 0, context.drawingBufferWidth, context.drawingBufferHeight, context.RED, context.FLOAT, results);
	console.log(results);
}

function createShader(context, type, source) {
	let shader = context.createShader(type);
	context.shaderSource(shader, source);
	context.compileShader(shader);
	let success = context.getShaderParameter(shader, context.COMPILE_STATUS);
	if (success) {
		return shader;
	}

	let errorMsg = context.getShaderInfoLog(shader);
	context.deleteShader(shader);
	throw new Error("in shader, " + errorMsg + " \n" + source);
}

function linkShaders(context, vertexShader, fragmentShader) {
	let program = context.createProgram();
	context.attachShader(program, vertexShader);
	context.attachShader(program, fragmentShader);
	context.linkProgram(program);
	let success = context.getProgramParameter(program, context.LINK_STATUS);
	if (success) {
		return program;
	}

	let errorMsg = context.getProgramInfoLog(program);
	context.deleteProgram(program);
	throw new Error("in linkShaders: " + errorMsg);
}