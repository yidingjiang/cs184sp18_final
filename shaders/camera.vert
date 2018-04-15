#version 330

uniform mat4 model;
uniform mat4 viewProjection;

layout(location = 0) in vec4 in_position;
in vec4 in_normal;

out vec4 vertex;
out vec4 normal;

void main() {
  gl_PointSize = 5.0;
  gl_Position = viewProjection * model * in_position;

  vertex = in_position;
  normal = in_normal;
}
