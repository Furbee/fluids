#version 330 core
// Texture samplers
uniform sampler2D ourTexture1;

in vec2 TexCoord;

layout(location = 0) out vec4 FragColor;

void main()
{
	FragColor = texture(ourTexture1, TexCoord);
}