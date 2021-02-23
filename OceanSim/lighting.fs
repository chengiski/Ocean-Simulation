#version 330 core
out vec4 FragColor;

in VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    noperspective float Jacob;
    float Jacobf;
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
	vec2 ndc = ((fs_in.ClipSpace.xy+fs_in.FragPos.y)/fs_in.ClipSpace.w)/2.0+0.5;
    //vec3 color = texture(floorTexture, fs_in.TexCoords).rgb;
	vec2 reflectTexCoords = vec2(ndc.x, -ndc.y);
	//vec3 texColor = texture(floorTexture, fs_in.TexCoords).rgb;
	vec3 reflectionColor = texture(reflectionTexture, reflectTexCoords).rgb;
	//vec4 refractionColor = texture(refractionTexture, fs_in.TexCoords);

    // Height Color
	vec3 shallowColor = vec3(0.1, 0.86, 0.92);
	vec3 deepColor = vec3(0.02, 0.09, 0.2);
	
    float relativeHeight;	// from 0 to 1
	relativeHeight = (fs_in.FragPos.y - heightMin) / (heightMax - heightMin);
	vec3 heightColor = relativeHeight * shallowColor + (1 - relativeHeight) * deepColor;
	// heightColor = vec3(s);	// Black and white

    vec3 Iv = normalize(fs_in.FragPos - viewPos);
    vec3 R = reflect(Iv, normalize(fs_in.Normal));
    vec3 skyboxColor = texture(skybox, R).rgb;

	skyboxColor = skyboxColor + reflectionColor;
    skyboxColor = mix(skyboxColor, heightColor, 0.6);

    // ambient
    vec3 ambient = 0.6 * heightColor;
    // diffuse
    vec3 lightDir = normalize(lightPos);// - fs_in.FragPos);
    //vec3 lightDir = normalize(lightPos - fs_in.FragPos); 
    vec3 normal = normalize(fs_in.Normal);
    float diff = max(dot(lightDir, normal), 0.0);
    //vec3 diffuse = diff * color;
    vec3 diffuse = diff * heightColor;
    // specular
    vec3 viewDir = normalize(viewPos - fs_in.FragPos);
    if (dot(normal, viewDir) < 0) normal = -normal;
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = 0.0;
    if(blinn)
    {
        vec3 halfwayDir = normalize(lightDir + viewDir);  
        spec = pow(max(dot(normal, halfwayDir), 0.0), 64.0);
    }
    else
    {
        vec3 reflectDir = reflect(-lightDir, normal);
        spec = pow(max(dot(viewDir, reflectDir), 0.0), 64);
    }
    vec3 specular = heightColor * spec; // assuming bright white light color
	
	float Distortion = 0.522;
	float Power = 2.28;
	float Scale = .918;
 
	vec3 H = normalize(lightDir + normal * Distortion);
	float I = pow(clamp(dot(viewDir, -H), 0.0, 1.0), Power) * Scale;
	vec3 SSS = clamp(vec3(I),0,1);

    float fresnel = clamp(fresnelFunc(0.16f, dot(normal, viewDir), 5.0f),0,1);

	FragColor = vec4((ambient + diffuse + specular + SSS), 1.0);
    FragColor = fresnel * vec4(skyboxColor,1.0f) + (1 - fresnel) * FragColor;

	if(fs_in.Jacob < -0.01) {
		FragColor = mix(vec4(1,1,1.0f,1),FragColor,0.7);
	}
    
	// if(fs_in.Jacob < -0.0001) {
	// 	//FragColor = mix(vec4(1,1,1.0f,1),FragColor,0.5);
	// }
    // if(fs_in.Jacob < -0.00000001) {
	// 	//FragColor = mix(vec4(1,1,1.0f,1),FragColor,0.5)
	// }

//     vec3 combinedColor = heightColor + reflectionColor;  

//     spec = clamp(spec, 0, 1);
// 	combinedColor *= (1 - spec);
// 	combinedColor += specular;	

//     FragColor = vec4(combinedColor, 1.0f); 
}