#version 330

uniform vec4 in_color;
uniform vec3 eye;
uniform vec3 light;

in vec4 vertex;
in vec4 normal;
in vec4 particle_color;
// in vec4 gl_FragCoord;

out vec4 out_color;

vec4 shadePhong() {
  // float p = 8.0;
  //
  // vec4 color = particle_color * 0.35;
  //
  // vec3 lightVec = light - vertex.xyz;
  // vec3 lightDir = normalize(lightVec);
  //
  // vec3 outDir = normalize(eye - vertex.xyz);
  // vec3 n = normalize(normal.xyz);

  vec3 N;
  N.xy = gl_PointCoord* 2.0 - vec2(1.0);
  float mag = dot(N.xy, N.xy);
  N.z = sqrt(1.0-mag);
  float diffuse = max(0.0, dot(vec3(0,0,1), N));

  //vec3 eye1 = vec3 (0.0, 0.0, 1.0);
  vec3 halfVector = normalize( eye + vec3(0,0,1));
  float spec = max( pow(dot(N,halfVector), 250.0), 0.);
  vec4 S = vec4(1)* spec;
  return particle_color * diffuse + S;


  // float distFactor = 1.0 / sqrt(dot(lightVec, lightVec));
  //
  // vec4 ambient = color * 0.9;
  // ambient.a = 0.5;
  //
  // float diffuseDot = dot(n, lightDir);
  // vec4 diffuse = color * clamp(diffuseDot, 0.0, 1.0);
  //
  // vec3 halfAngle = normalize(outDir + lightDir);
  // vec4 specularColor = min(color + 0.2, 1.0);
  // float specularDot = dot(n, halfAngle);
  // vec4 specular = 0.5 * specularColor * pow(clamp(specularDot, 0.0, 1.0), p);
  //
  // return diffuse + ambient + specular;
}

void main() {
  vec3 N;
  N.xy = gl_PointCoord* 2.0 - vec2(1.0);
  float mag = dot(N.xy, N.xy);
  if (mag > 1.0) discard;   // kill pixels outside circle

  out_color = shadePhong();
  // out_color = particle_color;
  out_color.a = in_color.a;
}
