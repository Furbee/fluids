#version 330 core
layout (location = 0) in vec4 position;
layout (location = 1) in vec2 texCoord;

out vec2 TexCoord;

void main()
{
	gl_Position = position;
    TexCoord = vec2(texCoord.x, 1.0 - texCoord.y);
}

