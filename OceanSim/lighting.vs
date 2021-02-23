#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec3 aJacob;

// declare an interface block; see 'Advanced GLSL' for what these are.
out VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    noperspective float Jacob;
    float Jacobf;
	vec4 ClipSpace;
} vs_out;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

void main()
{
    vs_out.FragPos = vec3(model * vec4(aPos, 1.0f));
    vec3 J= vec3(model * vec4(aJacob, 1.0f));
    float detJ = J.x *J.y - J.z*J.z;
    vs_out.Jacob = detJ;
    vs_out.Jacobf = detJ;
    vs_out.Normal = mat3(transpose(inverse(model))) * aNormal;  
    //vs_out.TexCoords = aTexCoords;
	vs_out.ClipSpace = projection * view * model * vec4(aPos, 1.0);
    gl_Position = vs_out.ClipSpace;
}