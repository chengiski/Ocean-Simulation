#version 330 core
out vec4 FragColor;

in VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    float Jacob;
	vec4 ClipSpace;
} fs_in;

uniform sampler2D floorTexture;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform bool blinn;

uniform sampler2D reflectionTexture;
//uniform sampler2D refractionTexture;

uniform float heightMax;
uniform float heightMin;

uniform samplerCube skybox;

float fresnelFunc(float f0, float nv, float p) {
	return f0 + (1 - f0) * pow(1 - nv, p);
}

void main()
{           
    vec3 normal = normalize(fs_in.Normal);
    vec3 lightDir = normalize(lightPos); 
    vec3 viewDir = normalize(viewPos - fs_in.FragPos);
    // Height Color
	vec3 shallowColor = vec3(0.1, 0.86, 0.92);
    vec3 deepColor = vec3(0.02, 0.09, 0.2);
	float relativeHeight = (fs_in.FragPos.y - heightMin) / (heightMax - heightMin); // from 0 to 1
    vec3 heightColor = relativeHeight * shallowColor + (1 - relativeHeight) * deepColor;

    // check sign of NL and NE
    if (dot(normal, lightDir) < 0) 
        normal = -normal;

    // SSS
    float Distortion = 0.522;
    float Power = 2.28;
    float Scale = 0.918;
    vec3 H = normalize(lightDir + normal * Distortion);
    float I = pow(clamp(dot(viewDir, -H), 0.0, 1.0), Power) * Scale;
    vec3 SSS = clamp(vec3(I),0,1);

    // Phong shading
    // ambient
    vec3 ambient = heightColor * vec3(0.6f);
    // diffuse
    float NL = dot(normal, lightDir);
    vec3 diffuse = heightColor * vec3(NL) * vec3(0.8f);
    // specular
    vec3 reflectDir = reflect(-lightDir, normal);
    float RE = clamp(dot(reflectDir, viewDir), 0, 1);
    float spec = pow(RE, 64) * 2;  // Over exposure
    vec3 specular = vec3(0.8f) * spec;
    
    vec3 phongColor = ambient + diffuse + specular;

    // box color
    vec2 ndc = (fs_in.ClipSpace.xy/fs_in.ClipSpace.w)/2.0+0.5;
    vec2 reflectTexCoords = vec2(ndc.x+fs_in.FragPos.y/fs_in.ClipSpace.w, -ndc.y-fs_in.FragPos.y/fs_in.ClipSpace.w);
    vec3 reflectionColor = texture(reflectionTexture, reflectTexCoords).rgb;

    // skybox color
    vec3 Iv = normalize(fs_in.FragPos - viewPos);
    vec3 R = reflect(Iv, normal);
    vec3 skyboxColor = texture(skybox, R).rgb;
    
    // Fresnel equation
    float fresnel = fresnelFunc(0.05f, dot(normal, viewDir), 1.0f);


    vec3 skyAndBox = reflectionColor + skyboxColor;

    vec3 finalColor = mix(phongColor, skyAndBox, 0.3);
    FragColor = vec4(finalColor, 1.0);


/*    if(fs_in.Jacob > -0.0001) {
        
    }
    else if (fs_in.Jacob < -0.002) {
        FragColor = vec4(1,1,1,1);
    }
    else {
        //float alpha = (detJ + 0.2) / 0.15;
        //heightColor = mix(originColor, vec3(1.0, 1.0, 1.0), alpha);
        //float WR = 0.217 * log(-fs_in.Jacob/0.099);
        //heightColor = mix(originColor, vec3(1.0, 1.0, 1.0), WR);
        float alpha = (-0.0001 - fs_in.Jacob) / 0.0019;
        float WR = alpha * alpha;
        //vec3 foam_color = mix(originColor, vec3(1.0, 1.0, 1.0), WR);
        //heightColor = mix(foam_color, vec3(1.0, 1.0, 1.0), foam_texture);
        FragColor = mix(vec4(finalColor, 1.0), vec4(1,1,1,1), alpha);
    }*/
}