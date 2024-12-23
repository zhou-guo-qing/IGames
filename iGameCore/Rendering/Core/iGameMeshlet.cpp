#include "iGameMeshlet.h"

IGAME_NAMESPACE_BEGIN

#ifdef IGAME_OPENGL_VERSION_460

void Meshlet::CreateBuffer() {
    m_MeshletsBuffer.create();
    m_DrawCommandBuffer.create();
    m_VisibleMeshletBuffer.create();
    m_FinalDrawCommandBuffer.create();
}

void Meshlet::BuildMeshlet(const float* vertex_positions, size_t vertex_count,
                           const int* indices, size_t index_count,
                           GLBuffer& EBO) {
    iGameTimer timer;
    timer.reset();

    // use for mesh shader
    std::vector<unsigned int> meshletVertices;
    std::vector<unsigned char> meshletTriangles;

    // Preprocessing: Vertex Cache Optimization
    std::vector<int> optimized_indices(index_count);
    meshopt_optimizeVertexCache(optimized_indices.data(), indices, index_count,
                                vertex_count);

    // Compute the upper bound of the meshlet
    size_t max_meshlets = meshopt_buildMeshletsBound(index_count, m_MaxVertices,
                                                     m_MaxTriangles);

    // Allocate meshlet data structure
    std::vector<meshopt_Meshlet> meshlets(max_meshlets);
    // Triangle index
    meshletVertices.resize(max_meshlets * m_MaxVertices);
    // Index of triangle index
    meshletTriangles.resize(max_meshlets * m_MaxTriangles * 3);

    // Generate meshlet data
    size_t meshlet_count = meshopt_buildMeshlets(
            meshlets.data(), meshletVertices.data(), meshletTriangles.data(),
            optimized_indices.data(), index_count, vertex_positions,
            vertex_count, sizeof(float) * 3, m_MaxVertices, m_MaxTriangles,
            m_ConeWeight);

    // Resize the vector to fit the actual generated meshlet data
    const meshopt_Meshlet& last = meshlets[meshlet_count - 1];
    meshletVertices.resize(last.vertex_offset + last.vertex_count);
    meshletTriangles.resize(last.triangle_offset +
                            ((last.triangle_count * 3 + 3) & ~3));
    meshlets.resize(meshlet_count);

    // Optimize vertex and index data for each meshlet
    for (size_t i = 0; i < meshlet_count; ++i) {
        const meshopt_Meshlet& m = meshlets[i];
        meshopt_optimizeMeshlet(&meshletVertices[m.vertex_offset],
                                &meshletTriangles[m.triangle_offset],
                                m.triangle_count, m.vertex_count);
    }

    // Calculate the boundary data of meshlet for cluster culling
    std::vector<meshopt_Bounds> meshlet_bounds(meshlet_count);
    for (size_t i = 0; i < meshlet_count; ++i) {
        const meshopt_Meshlet& m = meshlets[i];
        meshlet_bounds[i] = meshopt_computeMeshletBounds(
                &meshletVertices[m.vertex_offset],
                &meshletTriangles[m.triangle_offset], m.triangle_count,
                vertex_positions, vertex_count, sizeof(float) * 3);
    }

    // Record indirect Command
    std::vector<unsigned int> meshletIndices;
    std::vector<MeshletData> meshletDatas(meshlet_count);
    std::vector<DrawElementsIndirectCommand> drawCommands(meshlet_count);
    meshletIndices.resize(meshletTriangles.size());
    for (size_t i = 0; i < meshlet_count; ++i) {
        const meshopt_Meshlet& m = meshlets[i];
        const meshopt_Bounds& b = meshlet_bounds[i];

        meshletDatas[i] = {
                igm::vec4{b.center[0], b.center[1], b.center[2], b.radius},
                igm::vec4{0.0f}};
        drawCommands[i] = {m.triangle_count * 3, 0, m.triangle_offset, 0, 0};

        for (auto j = m.triangle_offset;
             j < m.triangle_offset + m.triangle_count * 3; j++) {
            meshletIndices[j] =
                    meshletVertices[m.vertex_offset + meshletTriangles[j]];
        }
    }

    GLAllocateGLBuffer(EBO, meshletIndices.size() * sizeof(igIndex),
                       meshletIndices.data());

    m_MeshletsCount = meshlet_count;
    m_MeshletsBuffer.target(GL_SHADER_STORAGE_BUFFER);
    m_MeshletsBuffer.allocate(meshletDatas.size() * sizeof(MeshletData),
                              meshletDatas.data(), GL_STATIC_DRAW);

    m_DrawCommandBuffer.target(GL_SHADER_STORAGE_BUFFER);
    m_DrawCommandBuffer.allocate(drawCommands.size() *
                                         sizeof(DrawElementsIndirectCommand),
                                 drawCommands.data(), GL_STATIC_DRAW);

    m_VisibleMeshletBuffer.target(GL_SHADER_STORAGE_BUFFER);
    m_VisibleMeshletBuffer.allocate(drawCommands.size() * sizeof(unsigned int),
                                    nullptr, GL_DYNAMIC_DRAW);

    m_FinalDrawCommandBuffer.target(GL_DRAW_INDIRECT_BUFFER);
    m_FinalDrawCommandBuffer.allocate(
            drawCommands.size() * sizeof(DrawElementsIndirectCommand), nullptr,
            GL_DYNAMIC_DRAW);

    std::cout << std::format("Build meshlets [count: {}, time: {}]",
                             meshlet_count,
                             FormatTime(timer.elapsedMilliseconds()))
              << std::endl;
}

size_t Meshlet::MeshletsCount() { return m_MeshletsCount; };

GLBuffer& Meshlet::MeshletsBuffer() { return m_MeshletsBuffer; };

GLBuffer& Meshlet::VisibleMeshletBuffer() { return m_VisibleMeshletBuffer; }

GLBuffer& Meshlet::DrawCommandBuffer() { return m_DrawCommandBuffer; };

GLBuffer& Meshlet::FinalDrawCommandBuffer() {
    return m_FinalDrawCommandBuffer;
};

#endif

IGAME_NAMESPACE_END