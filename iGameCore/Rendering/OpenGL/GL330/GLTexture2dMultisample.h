#pragma once

#include "GLObject.h"

IGAME_NAMESPACE_BEGIN

class GLTexture2dMultisample : public GLObject<GLTexture2dMultisample> {
private:
    friend class GLObject<GLTexture2dMultisample>;
    static void createHandle(GLsizei count, GLuint* handles) {
        glGenTextures(count, handles);
    }
    static void destroyHandle(GLsizei count, GLuint* handles) {
        glDeleteTextures(count, handles);
    }

private:
    GLTexture2dMultisample(GLuint handle)
        : GLObject<GLTexture2dMultisample>{handle} {}

public:
    GLTexture2dMultisample() = default;

public:
    // GLenum internal_format: https://registry.khronos.org/OpenGL-Refpages/gl4/html/glTexImage2D.xhtml
    // Sized Internal Format: GL_R8, GL_RG8, GL_RGB8, GL_RGBA8
    // Sized Internal Format: GL_DEPTH_COMPONENT32F, GL_DEPTH_COMPONENT24
    // Sized Internal Format: GL_DEPTH32F_STENCIL8, GL_DEPTH24_STENCIL8
    // Sized Internal Format: GL_STENCIL_INDEX8
    void storage(unsigned samples, GLenum internal_format, unsigned width,
                 unsigned height, bool fixedsamplelocations) const {
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, handle);
        glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, samples,
                                internal_format, width, height,
                                fixedsamplelocations);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
    }

    // GLenum texture: GL_TEXTURE1 - GL_TEXTURE15
    // GL_TEXTURE0 is reserved to prevent other binding operations from being performed after a texture unit is activated.
    void active(GLenum texture) {
        if (texture == GL_TEXTURE0) {
            throw std::runtime_error("GL_TEXTURE0 is reserved.");
        }
        glActiveTexture(texture);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, handle);
        glActiveTexture(GL_TEXTURE0);
    };

    // GLenum pname: GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_TEXTURE_MIN_FILTER, GL_TEXTURE_MAG_FILTER
    // GLint param: GL_CLAMP_TO_EDGE, GL_MIRRORED_REPEAT, GL_NEAREST, GL_LINEAR
    void parameteri(GLenum pname, GLint param) {
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, handle);
        glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, pname, param);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
    }

    void bind() const { glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, handle); }
    void release() const { glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0); }
};

IGAME_NAMESPACE_END