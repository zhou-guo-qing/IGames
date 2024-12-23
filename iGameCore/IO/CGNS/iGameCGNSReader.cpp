#if defined(CGNS_ENABLE)
#include "iGameCGNSReader.h"

IGAME_NAMESPACE_BEGIN
iGameCGNSReader::iGameCGNSReader()
{
}
iGameCGNSReader::~iGameCGNSReader()
{

}

void CollectNodeIDs(double parent_id, double node_id, int depth, std::vector<double>& node_ids) {
	char label[CGIO_MAX_NAME_LENGTH + 1];
	int numChildren;
	// Store the current node ID in the vector
	node_ids.push_back(node_id);
	// Get the label of the current node (optional, for debugging)
	cgio_get_label(parent_id, node_id, label);
	//std::cout << std::string(depth * 2, ' ') << node_id << " " << "Node Label: " << label << std::endl;
	// Get the number of children of the current node
	cgio_number_children(parent_id, node_id, &numChildren);
	if (numChildren > 0) {
		std::vector<double> child_ids(numChildren);
		cgio_children_ids(parent_id, node_id, 1, numChildren, &numChildren, child_ids.data());
		for (int i = 0; i < numChildren; ++i) {
			CollectNodeIDs(parent_id, child_ids[i], depth + 1, node_ids);
		}
	}
}
DataObject::Pointer iGameCGNSReader::ReadFile(std::string fileName) {

	int result;
	int file_type;
	result = cgio_check_file(fileName.c_str(), &file_type);

	int index_file;
	result = cg_open(fileName.data(), CG_MODE_READ, &index_file);
	if (CG_OK != result) {
		std::cout << "Open Error: " << cg_get_error() << std::endl;
	}
	else {
		std::cout << "Success to open cgns file!" << std::endl;
	}

	int cgio_num;
	result = cg_get_cgio(index_file, &cgio_num);
	if (CG_OK != result) {
		std::cout << "Get cgio num Error: " << cg_get_error() << std::endl;
	}
	// Get the root node ID
	double root_id;
	cgio_get_root_id(cgio_num, &root_id);
	// Vector to store all node IDs
	std::vector<double> node_ids;
	// Collect all node IDs starting from the root node
	CollectNodeIDs(cgio_num, root_id, 0, node_ids);

	for (int i = 0; i < node_ids.size(); i++) {
		char label[256];
		cgio_get_label(cgio_num, node_ids[i], label);
		if (strcmp(label, "ZoneBC_t") != 0) {
			continue;
		}
		double zoneId = node_ids[i];
		std::vector<double> child_node_ids;
		CollectNodeIDs(cgio_num, zoneId, 0, child_node_ids);
		for (int j = 0; j < child_node_ids.size(); j++) {
			cgio_get_label(cgio_num, child_node_ids[j], label);
			if (strcmp(label, "BC_t") != 0) {
				continue;
			}
			auto bcnode = child_node_ids[j];
			char name[256];
			char dataType[256];
			cgio_get_name(cgio_num, bcnode, name);
			cgio_get_data_type(cgio_num, bcnode, dataType);
			BoundryNames.insert(name);
			continue;
			//std::cout << name << ' ' << dataType << '\n';
			//std::string bctype;
			//std::string FamilyName;
			//CGNS_ENUMT(GridLocation_t) Location;
			//std::vector<int64_t> BCElementList;
			//std::vector<int64_t> BCElementRange;
			//readNodeStringData(cgio_num, bcnode, bctype);
			//std::cout << bctype << '\n';
			//if (bctype != "FamilySpecified")
			//{
			//	// waiting for c++20 to be replaced by starts_with
			//	if (bctype.rfind("BCWall", 0) == 0) {
			//		// Temporary Fallback for WALL bcs for old meshes
			//		FamilyName = "Wall";
			//	}
			//	else continue;
			//}
			//std::vector<double> childrenIds;
			//CollectNodeIDs(cgio_num, bcnode, 0, childrenIds);
			//for (auto iter = childrenIds.begin(); iter != childrenIds.end(); ++iter)
			//{
			//	char nodeName[CGIO_MAX_NAME_LENGTH + 1];
			//	char nodeLabel[CGIO_MAX_LABEL_LENGTH + 1];
			//	cgio_get_name(cgio_num, *iter, nodeName);
			//	cgio_get_label(cgio_num, *iter, nodeLabel);
			//	std::cout << nodeName << ' ' << nodeLabel << '\n';
			//	char dtype[256];
			//	if (strcmp(nodeName, "PointList") == 0)
			//	{
			//	}
			//	else if (strcmp(nodeName, "PointRange") == 0)
			//	{
			//		readNodeDataAs(cgio_num, *iter, BCElementRange);
			//	}
			//	else if (strcmp(nodeName, "ElementRange") == 0)
			//	{
			//	}
			//	else if (strcmp(nodeLabel, "FamilyName_t") == 0)
			//	{
			//	}
			//	else if (strcmp(nodeLabel, "GridLocation_t") == 0)
			//	{
			//		std::string location;
			//		readNodeStringData(cgio_num, *iter, location);
			//		if (location == "Vertex")
			//		{
			//			Location = CGNS_ENUMV(Vertex);
			//		}
			//		else if (location == "FaceCenter")
			//		{
			//			Location = CGNS_ENUMV(FaceCenter);
			//		}
			//		else if (location == "EdgeCenter")
			//		{
			//			Location = CGNS_ENUMV(EdgeCenter);
			//		}
			//	}
			//}
			//if (BCElementRange.size() == 2) {
			//}
		}
	}

	int nbases;
	result = cg_nbases(index_file, &nbases);
	if (CG_OK != result) {
		std::cout << "Get nbases Error: " << cg_get_error() << std::endl;
	}
	else {
		std::cout << "Base Num = " << nbases << std::endl;
	}

	for (int index_base = 1; index_base <= nbases; index_base++) {
		char basename[100];
		int celldim, physdim;
		result = cg_base_read(index_file, index_base, basename, &celldim, &physdim);
		if (CG_OK != result) {
			std::cout << "Read Base: " << cg_get_error() << std::endl;
		}
		else {
			std::cout << "BaseName = " << basename << std::endl;
			std::cout << "Dimension of the cells = " << celldim << std::endl;
			std::cout << "Number of coordinates required to define a vector in the field is " << physdim << std::endl;
		}
		int nzones;
		result = cg_nzones(index_file, index_base, &nzones);
		if (CG_OK != result) {
			std::cout << "Get nzones Error: " << cg_get_error() << std::endl;
		}
		else {
			std::cout << "Zone Num of Base(" << index_base << ") = " << nzones << std::endl;
			for (int index_zone = 1; index_zone <= nzones; index_zone++) {
				DataObject::Pointer DataSet = DataObject::New();
				ZoneType_t zoneType;
				result = cg_zone_type(index_file, index_base, index_zone, &zoneType);
				if (CG_OK != result) {
					std::cout << "Get Zone Type Error: " << cg_get_error() << std::endl;
				}
				cgsize_t size[9] = { 0 };
				char zonename[100];
				result = cg_zone_read(index_file, index_base, index_zone, zonename, size);
				if (CG_OK != result) {
					std::cout << "Zone Read Error: " << cg_get_error() << std::endl;
				}
				else {
					std::cout << "ZoneName = " << zonename << std::endl;
					if (zoneType == Structured) {
						std::cout << "ZoneType = Structured" << std::endl;
						if (m_StructuredMesh == nullptr) {
							m_StructuredMesh = StructuredMesh::New();
						}
						this->m_DataObjectType = IG_STRUCTURED_MESH;
						// Calculate total points and cells for structured zones
						cgsize_t totalPoints = 1, totalCells = 1;
						for (int i = 0; i < celldim; ++i) {
							totalPoints *= size[i];
							if (size[celldim + i] == 0) {
								break;
							}
							totalCells *= size[celldim + i];
						}
						std::cout << "Total Points: " << totalPoints << std::endl;
						std::cout << "Total Cells: " << totalCells << std::endl;
						// Read vertices coordinates
						this->ReadPointCoordinates(totalPoints, physdim, index_file, index_base, index_zone, size);
						// Gen cells (connectivity)
						this->GenStructuredCellConnectivities(celldim, size);
					}
					else if (zoneType == Unstructured) {
						std::cout << "ZoneType = Unstructured" << std::endl;
						std::cout << "Total Points: " << size[0] << std::endl;
						std::cout << "Total Cells: " << size[1] << std::endl;
						// Read vertices coordinates
						this->ReadPointCoordinates(size[0], physdim, index_file, index_base, index_zone, size);
						// Read cells (connectivity)
						this->ReadUnstructuredCellConnectivities(index_file, index_base, index_zone, size[1]);
					}
					else {
						std::cout << "Unknown ZoneType!" << std::endl;
						this->m_DataObjectType = IG_NONE;
					}
					/* Flow Solution */
					int nsols;
					result = cg_nsols(index_file, index_base, index_zone, &nsols);
					if (CG_OK != result) {
						std::cout << "Get Num Of Solution Error: " << cg_get_error() << std::endl;
					}
					else {
						std::cout << "Solution Num = " << nsols << std::endl;
					}
					for (int index_sol = 1; index_sol <= nsols; index_sol++)
					{
						char solname[100];
						GridLocation_t location;
						result = cg_sol_info(index_file, index_base, index_zone, index_sol, solname, &location);
						if (CG_OK != result) {
							std::cout << "Read Solution Info Error: " << cg_get_error() << std::endl;
						}
						std::cout << "Solution Name = " << solname << std::endl;
						switch (location)
						{
						case GridLocation_t::Vertex:
							std::cout << "Solution Location = Vertex" << std::endl;
							break;
						case CellCenter:
							std::cout << "Solution Location = CellCenter" << std::endl;
							break;
						case IFaceCenter:
							std::cout << "Solution Location = IFaceCenter" << std::endl;
							break;
						case JFaceCenter:
							std::cout << "Solution Location = JFaceCenter" << std::endl;
							break;
						case KFaceCenter:
							std::cout << "Solution Location = KFaceCenter" << std::endl;
							break;
						default:
							std::cout << "Solution Location is bad data! Unknown!" << std::endl;
						}
						this->ReadFields(index_file, index_base, index_zone, index_sol, zoneType, celldim, location, size);
					}
					if (this->m_AttributeSet)
					{
						if (this->m_DataObjectType == IG_STRUCTURED_MESH) {
							this->m_StructuredMesh->SetAttributeSet(this->m_AttributeSet);
						}
						if (this->m_DataObjectType == IG_UNSTRUCTURED_MESH) {
							this->m_UnstructuredMesh->SetAttributeSet(this->m_AttributeSet);
						}
						if (this->m_DataObjectType == IG_VOLUME_MESH) {
							this->m_VolumeMesh->SetAttributeSet(this->m_AttributeSet);
						}
					}
					//if (this->PolyDataSet != nullptr) {
					//	TransformPolyHedron(this->PolyDataSet, this->MultiBlockDataSet);
					//}
					//else
					//	this->MultiBlockDataSet->AddBlockToBack(DataSet);
				}

				if (nzones > 1) {
					if (m_ParentObject == nullptr) {
						m_ParentObject = DataObject::New();
					}
					m_ParentObject->AddSubDataObject(this->GetOutput());
					this->m_Points = nullptr;
					this->m_StructuredMesh = nullptr;
					this->m_UnstructuredMesh = nullptr;
					this->m_VolumeMesh = nullptr;
					this->m_AttributeSet = nullptr;
					this->m_DataObjectType = IG_MULTIBLOCK_MESH;
				}
			}
		}
	}

	TransformVectorArray();
	cg_close(index_file);
	return GetOutput();

}
void iGameCGNSReader::ReadPointCoordinates(int pointNum, int positionDim, int index_file, int index_base, int index_zone, cgsize_t* size)
{
	if (m_Points == nullptr) {
		m_Points = Points::New();
	}
	m_Points->SetNumberOfPoints(pointNum);
	auto points = m_Points->RawPointer();
	for (int dim = 0; dim < positionDim; ++dim) {
		std::string coordName = "Coordinate" + std::string(1, 'X' + dim);
		DoubleArray::Pointer CoordData = DoubleArray::New();
		CoordData->Resize(pointNum);
		auto coordData = CoordData->RawPointer();
		cgsize_t range_min[3] = { 1,1,1 };
		cgsize_t range_max[3] = { size[0], size[1], size[2] };
		int result = cg_coord_read(index_file, index_base, index_zone, coordName.c_str(), RealDouble, range_min, range_max, coordData);
		if (CG_OK != result) {
			std::cout << "Read " << coordName << " Error: " << cg_get_error() << std::endl;
		}
		else {
			int idx = dim;
			for (cgsize_t j = 0; j < pointNum; ++j) {
				points[idx] = coordData[j];
				idx += 3;
			}
		}
	}
}
void iGameCGNSReader::GenStructuredCellConnectivities(cgsize_t cellDim, cgsize_t* size)
{
	if (m_StructuredMesh == nullptr) {
		m_StructuredMesh = StructuredMesh::New();
	}
	auto structType = QUAD_4;
	if (cellDim == 3 && size[2] > 1) {
		structType = HEXA_8;
	}
	m_StructuredMesh->SetPoints(m_Points);
	igIndex tmpSize[3] = { 0 };
	for (int i = 0; i < 3; i++) {
		tmpSize[i] = size[i];
	}
	m_StructuredMesh->SetDimensionSize(tmpSize);
	m_StructuredMesh->GenStructuredCellConnectivities();
}
void iGameCGNSReader::ReadUnstructuredCellConnectivities(int index_file, int index_base, int  index_zone, cgsize_t cellNum)
{
	//这边需要读取出邻接关系然后转换为我们的type
	int nsections;
	int result = cg_nsections(index_file, index_base, index_zone, &nsections);
	if (CG_OK != result) {
		std::cout << "Get Num Of Sections Error: " << cg_get_error() << std::endl;
		return;
	}
	//paraview只读取了部分的section，应该是非结构化网格只有一个section
	//nsections = 1;
	igDebug("the number of sections is " << nsections);
	int cgio_num;
	if (cg_get_cgio(index_file, &cgio_num) != CG_OK) {
		std::cout << "Get cgio num Error: " << cg_get_error() << std::endl;
	}
	bool hasPolyGon = false;
	bool hasPolyHedron = false;
	bool hasUnstructGrid = false;
	for (int index_section = 1; index_section <= nsections; index_section++) {
		char sectionname[100];
		ElementType_t type;
		cgsize_t start, end;
		int nbndry, parent_flag;
		cg_section_read(index_file, index_base, index_zone, index_section, sectionname, &type, &start, &end, &nbndry, &parent_flag);
		cgsize_t num_elements = end - start + 1;
		std::cout << "section " << index_section << " " << sectionname << ' ' << start << ' ' << end << " " << num_elements << '\n';
		if (type == NGON_n)hasPolyGon = true;
		else if (type == NFACE_n)hasPolyHedron = true;
		else hasUnstructGrid = true;
	}
	if (hasPolyHedron) {
		if (!hasPolyGon || hasUnstructGrid) {
			igError("read Polyhedrons error");
		}
		else {
			if (m_VolumeMesh == nullptr) {
				m_VolumeMesh = VolumeMesh::New();
			}
			m_VolumeMesh->SetPoints(m_Points);
			this->m_DataObjectType = IG_VOLUME_MESH;
			CellArray::Pointer Faces = CellArray::New();
			CellArray::Pointer Volumes = CellArray::New();
			for (int index_section = 1; index_section <= nsections; index_section++) {
				char sectionname[100];
				ElementType_t type;
				cgsize_t start, end;
				int nbndry, parent_flag;
				cg_section_read(index_file, index_base, index_zone, index_section, sectionname, &type, &start, &end, &nbndry, &parent_flag);
				cgsize_t num_elements = end - start + 1;
				//m_VolumeMesh->AddSection(start - 1, end - 1, sectionname, BoundryNames.count(sectionname));
				if (type == NGON_n) {
					cgsize_t elementDataSize;
					cg_ElementDataSize(index_file, index_base, index_zone, index_section, &elementDataSize);
					std::vector<cgsize_t> elements;
					elements.resize(elementDataSize);
					auto elementsData = elements.data();
					cg_poly_elements_read(index_file, index_base, index_zone, index_section, elementsData, NULL, NULL);
					cgsize_t cellVcnt = 0;
					cgsize_t idx = 0;
					igIndex vhs[IGAME_CELL_MAX_SIZE];
					for (int i = 0; i < num_elements; i++) {
						cellVcnt = elements[idx++];
						for (int j = 0; j < cellVcnt; j++) {
							vhs[j] = elements[idx++] - 1;
						}
						Faces->AddCellIds(vhs, cellVcnt);
					}
					std::vector<cgsize_t> temp;
					elements.swap(temp);
				}
				else if (type == NFACE_n) {
					cgsize_t elementDataSize;
					cg_ElementDataSize(index_file, index_base, index_zone, index_section, &elementDataSize);
					std::vector<cgsize_t> elements;
					elements.resize(elementDataSize);
					auto elementsData = elements.data();
					cg_poly_elements_read(index_file, index_base, index_zone, index_section, elementsData, NULL, NULL);
					cgsize_t cellFcnt = 0;
					cgsize_t idx = 0;
					igIndex vhs[IGAME_CELL_MAX_SIZE];
					igIndex fhs[IGAME_CELL_MAX_SIZE];
					for (int i = 0; i < num_elements; i++) {
						cellFcnt = elements[idx++];
						for (int j = 0; j < cellFcnt; j++) {
							fhs[j] = std::abs(elements[idx++]) - 1;
						}
						Volumes->AddCellIds(fhs, cellFcnt);
					}
					std::vector<cgsize_t> temp;
					elements.swap(temp);
				}
			}
			this->m_VolumeMesh->InitVolumesWithPolyhedron(Faces, Volumes);
		}
	}
	else {
		if (m_UnstructuredMesh == nullptr) {
			m_UnstructuredMesh = UnstructuredMesh::New();
		}
		m_UnstructuredMesh->SetPoints(m_Points);
		this->m_DataObjectType = IG_UNSTRUCTURED_MESH;
		for (int index_section = 1; index_section <= nsections; index_section++) {
			char sectionname[100];
			ElementType_t type;
			cgsize_t start, end;
			int nbndry, parent_flag;
			cg_section_read(index_file, index_base, index_zone, index_section, sectionname, &type, &start, &end, &nbndry, &parent_flag);
			cgsize_t num_elements = end - start + 1;
			if (type == MIXED) {
				cgsize_t elementDataSize;
				cg_ElementDataSize(index_file, index_base, index_zone, index_section, &elementDataSize);
				std::vector<cgsize_t> elements(elementDataSize);
				cg_poly_elements_read(index_file, index_base, index_zone, index_section, elements.data(), NULL, NULL);
				this->ChangeMixElementToMyCell(elements, num_elements);
				std::vector<cgsize_t> temp;
				elements.swap(temp);
			}
			else {
				igIndex vhs[64];
				int cellVcnt = 0;
				cg_npe(type, &cellVcnt);
				std::vector<cgsize_t> elements(num_elements * cellVcnt);
				cg_elements_read(index_file, index_base, index_zone, index_section, elements.data(), NULL);
				if (type == HEXA_8) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_HEXAHEDRON);
					}
				}
				else if (type == TETRA_4) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_TETRA);
					}
				}
				else if (type == PYRA_5) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_PYRAMID);
					}
				}
				else if (type == PENTA_6) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_PRISM);
					}
				}
				else if (type == TRI_3) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_TRIANGLE);
					}
				}
				else if (type == QUAD_4) {
					for (cgsize_t j = 0; j < num_elements; ++j) {
						for (int k = 0; k < cellVcnt; k++) {
							vhs[k] = elements[j * cellVcnt + k] - 1;
						}
						m_UnstructuredMesh->AddCell(vhs, cellVcnt, IG_QUAD);
					}
				}
				std::vector<cgsize_t> temp;
				elements.swap(temp);
			}
		}
	}

}
void iGameCGNSReader::ChangeMixElementToMyCell(std::vector<cgsize_t>elements, int num_read_elements)
{
	if (m_UnstructuredMesh == nullptr) {
		m_UnstructuredMesh = UnstructuredMesh::New();
	}
	// Process each element
	cgsize_t idx = 0;
	for (int i = 0; i < num_read_elements; i++) {
		ElementType_t elem_type = static_cast<ElementType_t>(elements[idx]);
		int num_vertices = 0;
		igIndex vhs[64];
		cg_npe(elem_type, &num_vertices);
		for (int i = 0; i < num_vertices; ++i) {
			vhs[i] = elements[idx + i + 1] - 1;
		}
		switch (elem_type)
		{
		case TRI_3:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_TRIANGLE);
			break;
		case QUAD_4:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_QUAD);
			break;
		case TETRA_4:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_TETRA);
			break;
		case HEXA_8:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_HEXAHEDRON);
			break;
		case PYRA_5:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_PYRAMID);
			break;
		case PENTA_6:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_PRISM);
			break;
		case NGON_n:
			m_UnstructuredMesh->AddCell(vhs, num_vertices, IG_POLYGON);
		default:
			break;
		}
		idx += num_vertices + 1;
	}
}

void iGameCGNSReader::ReadFields(int index_file, int index_base, int  index_zone, int index_sol, ZoneType_t zoneType, int celldim, GridLocation_t location, cgsize_t* size)
{
	int nfileds;
	int result = cg_nfields(index_file, index_base, index_zone, index_sol, &nfileds);
	if (CG_OK != result) {
		std::cout << "Get Nfields Error: " << cg_get_error() << std::endl;
	}
	else {
		std::cout << "nfields = " << nfileds << std::endl;
		for (int index_field = 1; index_field <= nfileds; index_field++)
		{
			DataType_t dataType;
			char fieldname[100];
			result = cg_field_info(index_file, index_base, index_zone, index_sol, index_field, &dataType, fieldname);
			if (CG_OK != result) {
				std::cout << "Get Field Info Error: " << cg_get_error() << std::endl;
				continue;
			}
			else {
				std::cout << "FieldName = " << fieldname << std::endl;
				std::cout << "DataType_t(Integer, LongInteger, RealSingle, and RealDouble) = " << dataType << std::endl;
				cgsize_t range_min[3] = { 1 ,1, 1 };
				cgsize_t range_max[3];
				int arrayNum = 1;
				if (zoneType == Structured)
				{
					if (location == GridLocation_t::Vertex) {
						for (int i = 0; i < celldim; i++) {
							range_max[i] = size[i];
							if (size[i])arrayNum *= size[i];
						}
					}
					else {
						for (int i = 0; i < celldim; i++) {
							range_max[i] = size[i + celldim];
							if (size[i + celldim])arrayNum *= size[i + celldim];
						}
					}
				}
				else {
					if (location == GridLocation_t::Vertex) {
						range_max[0] = size[0];
						arrayNum = size[0];
					}
					else {
						range_max[0] = size[1];
						arrayNum = size[1];
					}
				}

				ArrayObject::Pointer solutionArray;
				switch (dataType)
				{
				case Integer:
				{
					auto solutionArray_2 = IntArray::New();
					solutionArray_2->Resize(arrayNum);
					result = cg_field_read(index_file, index_base, index_zone, index_sol, fieldname, dataType, range_min, range_max, solutionArray_2->RawPointer());
					solutionArray = solutionArray_2;
					break;
				}

				case LongInteger:
				{
					auto solutionArray_2 = LongLongArray::New();
					solutionArray_2->Resize(arrayNum);
					result = cg_field_read(index_file, index_base, index_zone, index_sol, fieldname, dataType, range_min, range_max, solutionArray_2->RawPointer());
					solutionArray = solutionArray_2;
					break;
				}
				case RealSingle:
				{
					auto solutionArray_2 = FloatArray::New();
					solutionArray_2->Resize(arrayNum);
					result = cg_field_read(index_file, index_base, index_zone, index_sol, fieldname, dataType, range_min, range_max, solutionArray_2->RawPointer());
					solutionArray = solutionArray_2;
					break;
				}
				case RealDouble:
				{
					auto solutionArray_2 = DoubleArray::New();
					solutionArray_2->Resize(arrayNum);
					result = cg_field_read(index_file, index_base, index_zone, index_sol, fieldname, dataType, range_min, range_max, solutionArray_2->RawPointer());
					solutionArray = solutionArray_2;
					break;
				}
				default:
					break;
				}
				if (CG_OK != result) {
					std::cout << "Get Field Error: " << cg_get_error() << std::endl;
				}
				else {
					//存储
					solutionArray->SetName(fieldname);
					if (m_AttributeSet == nullptr) {
						m_AttributeSet = AttributeSet::New();
					}
					if (location == GridLocation_t::Vertex) {
						m_AttributeSet->AddAttribute(IG_SCALAR, IG_POINT, solutionArray);
					}
					else if (location == CellCenter) {
						m_AttributeSet->AddAttribute(IG_SCALAR, IG_CELL, solutionArray);
					}
				}
			}
		}
	}
}
//void fun1(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto Polygons = iGameCellArray<iGamePolygon>::New();
//	auto PointData = Points::PointerData::New();
//	auto CellData = iGameCellData::New();
//	igIndex Bst = INT_MAX;
//	igIndex Bed = 0;
//	igIndex vhs[64];
//	igIndex vcnt;
//	igIndex fhs[64];
//	igIndex fcnt = 0;
//	std::vector<igIndex>f2c;
//	auto& sections = PolyDataSet->GetSections();
//	auto DataSet = iGameDataSet::New();
//	DataSet->SetPoints(PolyDataSet->GetPoints());
//	DataSet->SetPolygons(Polygons);
//	DataSet->SetPointData(PointData);
//	DataSet->SetCellData(CellData);
//
//
//	for (int i = 0; i < sections.size(); i++) {
//		auto& section = sections[i];
//		if (section.isBoundry) {
//			Bst = std::min(Bst, section.start);
//			Bed = std::max(Bed, section.end);
//		}
//	}
//	for (igIndex i = Bst; i <= Bed; i++) {
//		vcnt = PolyDataSet->GetPolygons()->GetCell(i, vhs);
//		Polygons->AddCell(vhs, vcnt);
//	}
//	f2c.resize(Bed - Bst + 1);
//	for (int i = 0; i < f2c.size(); i++) {
//		f2c[i] = -1;
//	}
//	for (int i = 0; i < PolyDataSet->GetNumberOfPolyhedrons(); i++) {
//		fcnt = PolyDataSet->GetPolyhedrons()->GetCell(i, fhs);
//		for (int j = 0; j < fcnt; j++) {
//			if (fhs[j] >= Bst && fhs[j] <= Bed) {
//				f2c[fhs[j] - Bst] = i;
//			}
//		}
//	}
//	for (int i = 0; i < f2c.size(); i++) {
//		if (f2c[i] == -1) {
//			std::cout << "f2c errror\n";
//		}
//	}
//	iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//	PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//	auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		pointNeiborFaceNum[i] = 0;
//	}
//	for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//		vcnt = Polygons->GetCell(i, vhs);
//		for (int j = 0; j < vcnt; j++) {
//			pointNeiborFaceNum[vhs[j]]++;
//		}
//	}
//	Points::Pointer* SurfacePoints = Points::Pointer::New();
//	igIndex SurfacePointNum = 0;
//	std::vector<igIndex>op2np;
//	op2np.resize(PolyDataSet->GetNumberOfPoints());
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		if (pointNeiborFaceNum[i] != 0) {
//			op2np[i] = SurfacePointNum;
//			SurfacePoints->AddPoint(PolyPoints->GetPoint(i));
//			SurfacePointNum++;
//		}
//	}
//	DataSet->SetPoints(SurfacePoints);
//	igIndex* tmparray = Polygons->GetBuffer()->Pointer;
//	for (int i = 0; i < Polygons->GetBuffer()->Size; i++) {
//		tmparray[i] = op2np[tmparray[i]];
//	}
//	auto AllArray = PolyDataSet->GetCellData()->GetAllScalars();
//	for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//		auto array = AllArray[arrayNum].Array;
//		iGameFloatArray* PointArray = iGameFloatArray::New();
//		PointArray->SetName(array->GetName());
//		PointArray->SetNumberOfValues(SurfacePointNum);
//		auto pointArray = PointArray->GetRawPointer();
//		for (int i = 0; i < SurfacePointNum; i++) {
//			pointArray[i] = 0.0;
//		}
//		float value = .0;
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			value = array->GetValueF(f2c[i]);
//			vcnt = Polygons->GetCell(i, vhs);
//			for (int j = 0; j < vcnt; j++) {
//				pointArray[vhs[j]] += value;
//			}
//		}
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			if (pointNeiborFaceNum[i] != 0)
//				pointArray[op2np[i]] /= pointNeiborFaceNum[i] * 1.0;
//		}
//		PointData->AddScalars(PointArray);
//	}
//	for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//		auto array = AllArray[arrayNum].Array;
//		iGameFloatArray* CellArray = iGameFloatArray::New();
//		CellArray->SetName(array->GetName());
//		CellArray->SetNumberOfValues(Bed - Bst + 1);
//		auto cellArray = CellArray->GetRawPointer();
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			cellArray[i] = array->GetValueF(f2c[i]);
//		}
//		CellData->AddScalars(CellArray);
//	}
//	delete PointNeiborFaceNum;
//	f2c.swap(std::vector<igIndex>());
//	MBDataset->AddBlockToBack(DataSet);
//}
//void fun2(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto Polygons = iGameCellArray<iGamePolygon>::New();
//	auto PointData = Points::PointerData::New();
//	auto CellData = iGameCellData::New();
//	igIndex Bst = INT_MAX;
//	igIndex Bed = 0;
//	igIndex vhs[64];
//	igIndex vcnt;
//	igIndex fhs[64];
//	igIndex fcnt = 0;
//	std::vector<igIndex>f2c;
//	auto& sections = PolyDataSet->GetSections();
//	auto DataSet = iGameDataSet::New();
//	DataSet->SetPoints(PolyDataSet->GetPoints());
//	DataSet->SetPolygons(Polygons);
//	DataSet->SetPointData(PointData);
//	DataSet->SetCellData(CellData);
//
//
//	for (int i = 0; i < sections.size(); i++) {
//		auto& section = sections[i];
//		if (section.isBoundry) {
//			Bst = std::min(Bst, section.start);
//			Bed = std::max(Bed, section.end);
//		}
//	}
//	for (igIndex i = Bst; i <= Bed; i++) {
//		vcnt = PolyDataSet->GetPolygons()->GetCell(i, vhs);
//		Polygons->AddCell(vhs, vcnt);
//	}
//	f2c.resize(Bed - Bst + 1);
//	for (int i = 0; i < f2c.size(); i++) {
//		f2c[i] = -1;
//	}
//	for (int i = 0; i < PolyDataSet->GetNumberOfPolyhedrons(); i++) {
//		fcnt = PolyDataSet->GetPolyhedrons()->GetCell(i, fhs);
//		for (int j = 0; j < fcnt; j++) {
//			if (fhs[j] >= Bst && fhs[j] <= Bed) {
//				f2c[fhs[j] - Bst] = i;
//			}
//		}
//	}
//	for (int i = 0; i < f2c.size(); i++) {
//		if (f2c[i] == -1) {
//			std::cout << "f2c errror\n";
//		}
//	}
//	iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//	PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//	auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		pointNeiborFaceNum[i] = 0;
//	}
//	for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//		vcnt = Polygons->GetCell(i, vhs);
//		for (int j = 0; j < vcnt; j++) {
//			pointNeiborFaceNum[vhs[j]]++;
//		}
//	}
//	auto AllArray = PolyDataSet->GetCellData()->GetAllScalars();
//	for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//		auto array = AllArray[arrayNum].Array;
//		iGameFloatArray* PointArray = iGameFloatArray::New();
//		PointArray->SetName(array->GetName());
//		PointArray->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//		auto pointArray = PointArray->GetRawPointer();
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			pointArray[i] = 0.0;
//		}
//		float value = .0;
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			value = array->GetValueF(f2c[i]);
//			vcnt = Polygons->GetCell(i, vhs);
//			for (int j = 0; j < vcnt; j++) {
//				pointArray[vhs[j]] += value;
//			}
//		}
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			if (pointNeiborFaceNum[i] != 0)
//				pointArray[i] /= pointNeiborFaceNum[i] * 1.0;
//		}
//		PointData->AddScalars(PointArray);
//	}
//	delete PointNeiborFaceNum;
//	f2c.swap(std::vector<igIndex>());
//	MBDataset->AddBlockToBack(DataSet);
//}
//void fun3(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto& sections = PolyDataSet->GetSections();
//	for (int sectionIndex = 5; sectionIndex < 14; sectionIndex++) {
//		if (!sections[sectionIndex].isBoundry)continue;
//		auto Polygons = iGameCellArray<iGamePolygon>::New();
//		auto PointData = Points::PointerData::New();
//		igIndex Bst = INT_MAX;
//		igIndex Bed = 0;
//		igIndex vhs[64];
//		igIndex vcnt;
//		igIndex fhs[64];
//		igIndex fcnt = 0;
//		std::vector<igIndex>f2c;
//		auto DataSet = iGameDataSet::New();
//		DataSet->SetName(sections[sectionIndex].name);
//		DataSet->SetPoints(PolyDataSet->GetPoints());
//		DataSet->SetPolygons(Polygons);
//		DataSet->SetPointData(PointData);
//		Bst = sections[sectionIndex].start;
//		Bed = sections[sectionIndex].end;
//		for (int i = Bst; i <= Bed; i++) {
//			vcnt = PolyDataSet->GetPolygons()->GetCell(i, vhs);
//			Polygons->AddCell(vhs, vcnt);
//		}
//		f2c.resize(Bed - Bst + 1);
//		for (int i = 0; i < f2c.size(); i++) {
//			f2c[i] = -1;
//		}
//		for (int i = 0; i < PolyDataSet->GetNumberOfPolyhedrons(); i++) {
//			fcnt = PolyDataSet->GetPolyhedrons()->GetCell(i, fhs);
//			for (int j = 0; j < fcnt; j++) {
//				if (fhs[j] >= Bst && fhs[j] <= Bed) {
//					f2c[fhs[j] - Bst] = i;
//				}
//			}
//		}
//		iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//		PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//		auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			pointNeiborFaceNum[i] = 0;
//		}
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			vcnt = Polygons->GetCell(i, vhs);
//			for (int j = 0; j < vcnt; j++) {
//				pointNeiborFaceNum[vhs[j]]++;
//			}
//		}
//
//		Points::Pointer* SurfacePoints = Points::Pointer::New();
//		igIndex SurfacePointNum = 0;
//		std::vector<igIndex>op2np;
//		op2np.resize(PolyDataSet->GetNumberOfPoints());
//
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			if (pointNeiborFaceNum[i] != 0) {
//				op2np[i] = SurfacePointNum;
//				SurfacePoints->AddPoint(PolyPoints->GetPoint(i));
//				SurfacePointNum++;
//			}
//		}
//		DataSet->SetPoints(SurfacePoints);
//		igIndex* tmparray = Polygons->GetBuffer()->Pointer;
//		for (int i = 0; i < Polygons->GetNumberOfIds(); i++) {
//			tmparray[i] = op2np[tmparray[i]];
//		}
//		auto AllArray = PolyDataSet->GetCellData()->GetAllScalars();
//		for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//			auto array = AllArray[arrayNum].Array;
//			iGameFloatArray* PointArray = iGameFloatArray::New();
//			PointArray->SetName(array->GetName());
//			PointArray->SetNumberOfValues(SurfacePointNum);
//			auto pointArray = PointArray->GetRawPointer();
//			for (int i = 0; i < SurfacePointNum; i++) {
//				pointArray[i] = 0.0;
//			}
//			float value = .0;
//			for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//				value = array->GetValueF(f2c[i]);
//				vcnt = Polygons->GetCell(i, vhs);
//				for (int j = 0; j < vcnt; j++) {
//					pointArray[vhs[j]] += value;
//				}
//			}
//			for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//				if (pointNeiborFaceNum[i] != 0)
//					pointArray[op2np[i]] /= pointNeiborFaceNum[i] * 1.0;
//			}
//			PointData->AddScalars(PointArray);
//		}
//		op2np.swap(std::vector<igIndex>());
//
//		delete PointNeiborFaceNum;
//		f2c.swap(std::vector<igIndex>());
//		MBDataset->AddBlockToBack(DataSet);
//	}
//
//}
//void fun4(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto& sections = PolyDataSet->GetSections();
//	for (int sectionIndex = 0; sectionIndex < sections.size(); sectionIndex++) {
//		if (!sections[sectionIndex].isBoundry)continue;
//		auto Polygons = iGameCellArray<iGamePolygon>::New();
//		auto PointData = Points::PointerData::New();
//		igIndex Bst = INT_MAX;
//		igIndex Bed = 0;
//		igIndex vhs[64];
//		igIndex vcnt;
//		igIndex fhs[64];
//		igIndex fcnt = 0;
//		std::vector<igIndex>f2c;
//		auto DataSet = iGameDataSet::New();
//		DataSet->SetName(sections[sectionIndex].name);
//		DataSet->SetPoints(PolyDataSet->GetPoints());
//		DataSet->SetPolygons(Polygons);
//		DataSet->SetPointData(PointData);
//		Bst = sections[sectionIndex].start;
//		Bed = sections[sectionIndex].end;
//		for (int i = Bst; i <= Bed; i++) {
//			vcnt = PolyDataSet->GetPolygons()->GetCell(i, vhs);
//			Polygons->AddCell(vhs, vcnt);
//		}
//		f2c.resize(Bed - Bst + 1);
//		for (int i = 0; i < f2c.size(); i++) {
//			f2c[i] = -1;
//		}
//		for (int i = 0; i < PolyDataSet->GetNumberOfPolyhedrons(); i++) {
//			fcnt = PolyDataSet->GetPolyhedrons()->GetCell(i, fhs);
//			for (int j = 0; j < fcnt; j++) {
//				if (fhs[j] >= Bst && fhs[j] <= Bed) {
//					f2c[fhs[j] - Bst] = i;
//				}
//			}
//		}
//		iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//		PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//		auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			pointNeiborFaceNum[i] = 0;
//		}
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			vcnt = Polygons->GetCell(i, vhs);
//			for (int j = 0; j < vcnt; j++) {
//				pointNeiborFaceNum[vhs[j]]++;
//			}
//		}
//		auto AllArray = PolyDataSet->GetCellData()->GetAllScalars();
//		for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//			auto array = AllArray[arrayNum].Array;
//			iGameFloatArray* PointArray = iGameFloatArray::New();
//			PointArray->SetName(array->GetName());
//			PointArray->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//			auto pointArray = PointArray->GetRawPointer();
//			for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//				pointArray[i] = 0.0;
//			}
//			float value = .0;
//			for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//				value = array->GetValueF(f2c[i]);
//				vcnt = Polygons->GetCell(i, vhs);
//				for (int j = 0; j < vcnt; j++) {
//					pointArray[vhs[j]] += value;
//				}
//			}
//			for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//				if (pointNeiborFaceNum[i] != 0)
//					pointArray[i] /= pointNeiborFaceNum[i] * 1.0;
//			}
//			PointData->AddScalars(PointArray);
//		}
//		delete PointNeiborFaceNum;
//		f2c.swap(std::vector<igIndex>());
//		MBDataset->AddBlockToBack(DataSet);
//	}
//
//}
//void fun5(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto Polygons = PolyDataSet->GetPolygons();
//	auto Polyhedrons = PolyDataSet->GetPolyhedrons();
//	auto CellVertices = PolyDataSet->GetCellVertices();
//	auto NewPolygons = iGameCellArray<iGamePolygon>::New();
//	igIndex vhs[64];
//	igIndex vcnt = 0, fcnt = 0, fvnt = 0;
//	igIndex fhs[64];
//	std::vector<igIndex>f2c;
//	auto& sections = PolyDataSet->GetSections();
//	auto DataSet = iGameDataSet::New();
//	DataSet->SetPoints(PolyDataSet->GetPoints());
//	DataSet->SetPolygons(NewPolygons);
//	auto PolyhedronNum = PolyDataSet->GetNumberOfPolyhedrons();
//	for (igIndex i = 0; i < PolyhedronNum; i++) {
//		fcnt = Polyhedrons->GetCell(i, fhs);
//		vcnt = CellVertices->GetCell(i, vhs);
//		if ((fcnt == 6 && vcnt == 8) || (fcnt == 4 && vcnt == 4)
//			|| (fcnt == 5 && vcnt == 6) || (fcnt == 5 && vcnt == 5)) {
//			continue;
//		}
//		for (igIndex j = 0; j < fcnt; j++) {
//			vcnt = Polygons->GetCell(fhs[j], vhs);
//			NewPolygons->AddCell(vhs, vcnt);
//		}
//	}
//	std::cout << NewPolygons->GetNumberOfCells() << '\n';
//	iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//	PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//	auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		pointNeiborFaceNum[i] = 0;
//	}
//	for (int i = 0; i < NewPolygons->GetNumberOfCells(); i++) {
//		vcnt = NewPolygons->GetCell(i, vhs);
//		for (int j = 0; j < vcnt; j++) {
//			pointNeiborFaceNum[vhs[j]]++;
//		}
//	}
//	Points::Pointer* SurfacePoints = Points::Pointer::New();
//	igIndex SurfacePointNum = 0;
//	std::vector<igIndex>op2np;
//	op2np.resize(PolyDataSet->GetNumberOfPoints());
//
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		if (pointNeiborFaceNum[i] != 0) {
//			op2np[i] = SurfacePointNum;
//			SurfacePoints->AddPoint(PolyPoints->GetPoint(i));
//			SurfacePointNum++;
//		}
//	}
//	DataSet->SetPoints(SurfacePoints);
//	igIndex* tmparray = NewPolygons->GetBuffer()->Pointer;
//	for (int i = 0; i < NewPolygons->GetNumberOfIds(); i++) {
//		tmparray[i] = op2np[tmparray[i]];
//	}
//	MBDataset->AddBlockToBack(DataSet);
//}
//void fun6(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto Polygons = PolyDataSet->GetPolygons();
//	auto Polyhedrons = PolyDataSet->GetPolyhedrons();
//	auto CellVertices = PolyDataSet->GetCellVertices();
//	auto NewPolyhedrons = iGameCellArray<iGamePolyhedron>::New();
//	igIndex vhs[256];
//	igIndex i = 0, j = 0, k = 0, vcnt = 0, fcnt = 0, fvcnt = 0;
//	igIndex fhs[256];
//	igIndex fvhs[256];
//
//	std::vector<igIndex>f2c;
//	auto& sections = PolyDataSet->GetSections();
//	auto DataSet = iGameDataSet::New();
//	DataSet->SetPoints(PolyDataSet->GetPoints());
//	DataSet->SetPolyhedrons(NewPolyhedrons);
//	//DataSet->SetCellData(PolyDataSet->GetCellData());
//	auto PolyhedronNum = PolyDataSet->GetNumberOfPolyhedrons();
//	std::cout << PolyhedronNum << '\n';
//	for (i = 0; i < PolyhedronNum; i++) {
//		vcnt = 0;
//		fcnt = Polyhedrons->GetCell(i, fhs);
//		vhs[vcnt++] = fcnt;
//		for (j = 0; j < fcnt; j++) {
//			fvcnt = Polygons->GetCell(fhs[j], fvhs);
//			vhs[vcnt++] = fvcnt;
//			for (k = 0; k < fvcnt; k++) {
//				vhs[vcnt++] = fvhs[k];
//			}
//		}
//		NewPolyhedrons->AddCell(vhs, vcnt);
//	}
//	MBDataset->AddBlockToBack(DataSet);
//}
//void fun7(iGamePolyDataSet* PolyDataSet, iGameMultiBlockDataSet* MBDataset) {
//	auto PolyPoints = PolyDataSet->GetPoints();
//	auto Polygons = iGameCellArray<iGamePolygon>::New();
//	auto PointData = Points::PointerData::New();
//	auto CellData = iGameCellData::New();
//	igIndex Bst = INT_MAX;
//	igIndex Bed = 0;
//	igIndex vhs[64];
//	igIndex vcnt;
//	igIndex fhs[64];
//	igIndex fcnt = 0;
//	std::vector<igIndex>f2c;
//	auto& sections = PolyDataSet->GetSections();
//	auto DataSet = iGameDataSet::New();
//	DataSet->SetPoints(PolyDataSet->GetPoints());
//	DataSet->SetPolygons(Polygons);
//	DataSet->SetPointData(PointData);
//	DataSet->SetCellData(CellData);
//
//
//	for (int i = 0; i < sections.size(); i++) {
//		auto& section = sections[i];
//		if (section.isBoundry) {
//			Bst = std::min(Bst, section.start);
//			Bed = std::max(Bed, section.end);
//		}
//	}
//	for (igIndex i = Bst; i <= Bed; i++) {
//		vcnt = PolyDataSet->GetPolygons()->GetCell(i, vhs);
//		Polygons->AddCell(vhs, vcnt);
//	}
//	f2c.resize(Bed - Bst + 1);
//	for (int i = 0; i < f2c.size(); i++) {
//		f2c[i] = -1;
//	}
//	for (int i = 0; i < PolyDataSet->GetNumberOfPolyhedrons(); i++) {
//		fcnt = PolyDataSet->GetPolyhedrons()->GetCell(i, fhs);
//		for (int j = 0; j < fcnt; j++) {
//			if (fhs[j] >= Bst && fhs[j] <= Bed) {
//				f2c[fhs[j] - Bst] = i;
//			}
//		}
//	}
//	for (int i = 0; i < f2c.size(); i++) {
//		if (f2c[i] == -1) {
//			std::cout << "f2c errror\n";
//		}
//	}
//	iGameIntArray* PointNeiborFaceNum = iGameIntArray::New();
//	PointNeiborFaceNum->SetNumberOfValues(PolyDataSet->GetNumberOfPoints());
//	auto pointNeiborFaceNum = PointNeiborFaceNum->GetRawPointer();
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		pointNeiborFaceNum[i] = 0;
//	}
//	for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//		vcnt = Polygons->GetCell(i, vhs);
//		for (int j = 0; j < vcnt; j++) {
//			pointNeiborFaceNum[vhs[j]]++;
//		}
//	}
//	Points::Pointer* SurfacePoints = Points::Pointer::New();
//	igIndex SurfacePointNum = 0;
//	std::vector<igIndex>op2np;
//	op2np.resize(PolyDataSet->GetNumberOfPoints());
//	for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//		if (pointNeiborFaceNum[i] != 0) {
//			op2np[i] = SurfacePointNum;
//			SurfacePoints->AddPoint(PolyPoints->GetPoint(i));
//			SurfacePointNum++;
//		}
//	}
//	DataSet->SetPoints(SurfacePoints);
//	igIndex* tmparray = Polygons->GetBuffer()->Pointer;
//	for (int i = 0; i < Polygons->GetBuffer()->Size; i++) {
//		tmparray[i] = op2np[tmparray[i]];
//	}
//	auto AllArray = PolyDataSet->GetCellData()->GetAllScalars();
//	for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//		auto array = AllArray[arrayNum].Array;
//		iGameFloatArray* PointArray = iGameFloatArray::New();
//		PointArray->SetName(array->GetName());
//		PointArray->SetNumberOfValues(SurfacePointNum);
//		auto pointArray = PointArray->GetRawPointer();
//		for (int i = 0; i < SurfacePointNum; i++) {
//			pointArray[i] = 0.0;
//		}
//		float value = .0;
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			value = array->GetValueF(f2c[i]);
//			vcnt = Polygons->GetCell(i, vhs);
//			for (int j = 0; j < vcnt; j++) {
//				pointArray[vhs[j]] += value;
//			}
//		}
//		for (int i = 0; i < PolyDataSet->GetNumberOfPoints(); i++) {
//			if (pointNeiborFaceNum[i] != 0)
//				pointArray[op2np[i]] /= pointNeiborFaceNum[i] * 1.0;
//		}
//		PointData->AddScalars(PointArray);
//	}
//	for (int arrayNum = 0; arrayNum < AllArray.size(); arrayNum++) {
//		auto array = AllArray[arrayNum].Array;
//		iGameFloatArray* CellArray = iGameFloatArray::New();
//		CellArray->SetName(array->GetName());
//		CellArray->SetNumberOfValues(Bed - Bst + 1);
//		auto cellArray = CellArray->GetRawPointer();
//		for (int i = 0; i < Polygons->GetNumberOfCells(); i++) {
//			cellArray[i] = array->GetValueF(f2c[i]);
//		}
//		CellData->AddScalars(CellArray);
//	}
//	delete PointNeiborFaceNum;
//	f2c.swap(std::vector<igIndex>());
//	MBDataset->AddBlockToBack(DataSet);
//}


void iGameCGNSReader::TransformVectorArray()
{
//	igIndex i, j, k;
//	
//	AttributeSet DataSet = this->m_AttributeSet;
//	igIndex arraySize = 0;
//	if (DataSet->attachmentType) {
//			auto AllScalars = DataSet->GetPointData()->;
//			Points::PointerData* NewPointData = Points::PointerData::New();
//			std::vector<std::string>scalarNames;
//			for (i = 0; i < AllScalars.size(); i++) {
//				auto scalarDataArray = AllScalars[i].Array;
//				arraySize = scalarDataArray->GetNumberOfTuples();
//				scalarNames.emplace_back(scalarDataArray->GetName());
//			}
//			for (i = 0; i < AllScalars.size(); i++) {
//
//				auto name = scalarNames[i];
//				bool isvector = false;
//				if (name[name.length() - 1] == 'X') {
//					isvector = true;
//					for (j = 1; j < 3; j++) {
//						auto tmpName = scalarNames[i + j];
//						if (tmpName[tmpName.length() - 1] != 'X' + j) {
//							isvector = false;
//						}
//					}
//				}
//				if (isvector) {
//					iGameFloatArray* Vector = iGameFloatArray::New();
//					Vector->SetName(name.substr(0, name.length() - 1));
//					Vector->SetNumberOfComponents(3);
//					Vector->SetNumberOfTuples(arraySize);
//					float* vector = Vector->GetRawPointer();
//					igIndex index = 0;
//					for (j = 0; j < 3; j++) {
//						auto scalarData = AllScalars[i + j].Array;
//						index = j;
//						for (k = 0; k < arraySize; k++) {
//							vector[index] = scalarData->GetValueF(k);
//							index += 3;
//						}
//					}
//					NewPointData->AddVectors(Vector);
//					i += 2;
//				}
//				else {
//					NewPointData->AddScalars(AllScalars[i].Array);
//				}
//				DataSet->SetPointData(NewPointData);
//			}
//		}
//
//		if (DataSet->GetCellData()) {
//			auto AllScalars = DataSet->GetCellData()->GetAllScalars();
//			iGameCellData* NewCellData = iGameCellData::New();
//			std::vector<std::string>scalarNames;
//			for (i = 0; i < AllScalars.size(); i++) {
//				auto scalarDataArray = AllScalars[i].Array;
//				arraySize = scalarDataArray->GetNumberOfTuples();
//				scalarNames.emplace_back(scalarDataArray->GetName());
//			}
//			for (i = 0; i < AllScalars.size(); i++) {
//
//				auto name = scalarNames[i];
//				bool isvector = false;
//				if (name[name.length() - 1] == 'X') {
//					isvector = true;
//					for (j = 1; j < 3; j++) {
//						auto tmpName = scalarNames[i + j];
//						if (tmpName[tmpName.length() - 1] != 'X' + j) {
//							isvector = false;
//						}
//					}
//				}
//				if (isvector) {
//					iGameFloatArray* Vector = iGameFloatArray::New();
//					Vector->SetName(name.substr(0, name.length() - 1));
//					Vector->SetNumberOfComponents(3);
//					Vector->SetNumberOfTuples(arraySize);
//					float* vector = Vector->GetRawPointer();
//					igIndex index = 0;
//					for (j = 0; j < 3; j++) {
//						auto scalarData = AllScalars[i + j].Array;
//						index = j;
//						for (k = 0; k < arraySize; k++) {
//							vector[index] = scalarData->GetValueF(k);
//							index += 3;
//						}
//						delete scalarData;
//					}
//					NewCellData->AddVectors(Vector);
//					i += 2;
//				}
//				else {
//					NewCellData->AddScalars(AllScalars[i].Array);
//				}
//			}
//			DataSet->SetCellData(NewCellData);
//		}
//
//	
}
IGAME_NAMESPACE_END
#endif