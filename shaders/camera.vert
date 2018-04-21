#version 330

uniform mat4 model;
uniform mat4 viewProjection;
uniform float particle_size;

layout(location = 0) in vec4 in_position;
layout(location = 1) in vec4 in_color;
in vec4 in_normal;

out vec4 vertex;
out vec4 normal;
out vec4 particle_color;

void main() {
  gl_PointSize = particle_size;
  gl_Position = viewProjection * model * in_position;

  vertex = in_position;
  normal = in_normal;
  particle_color = in_color;
}
