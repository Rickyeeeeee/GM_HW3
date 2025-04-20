#include "controller.h"
#include "learnply.h"
#include "meshprocessor.h"
#include <ImGUI/imgui.h>
#include <ImGUI/imgui_impl_glfw.h>
#include <ImGUI/imgui_impl_opengl3.h>
#include <chrono>
#include <iomanip>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdio.h>

Controller::Controller(int width, int height, GLFWwindow *window, Scene *scenePtr)
    : screenWidth(width), screenHeight(height), isDragging(false), isPanning(false), lastMouseX(0),
      lastMouseY(0), scene(scenePtr) {
  // Trackball initialization
  trackball = std::make_unique<Trackball>();
  trackball->setCamera(scene->getCamera());
  trackball->start(Trackball::Around);
  // GUI Initialization
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  io.FontGlobalScale = 2.0f;
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 130");
  resetControl();
}

void Controller::resize(int width, int height) {
  screenWidth = width;
  screenHeight = height;
}

void Controller::handleMouseClick(int button, double xpos, double ypos, int action) {
  // don't pass mouse and keyboard presses further if an ImGui widget is active
  auto& io = ImGui::GetIO();
  if (io.WantCaptureMouse) {
      return;
  }

  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    isDragging = (action == GLFW_PRESS);
    if (isDragging) {
      trackball->start(Trackball::Around);
    }
  } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
    isPanning = (action == GLFW_PRESS);
    if (isPanning) {
      lastMouseX = xpos;
      lastMouseY = ypos;
    }
  } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    if (selectionMode == (int)SelectType::Faces) {
      selectFaces(xpos, ypos);
    } else if (selectionMode == (int)SelectType::Vertices) {
      selectVertices(xpos, ypos);
    } else if (selectionMode == (int)SelectType::Edges) {
      selectEdges(xpos, ypos);
    }
  }
}

void Controller::handleMouseMove(double xpos, double ypos) {
  if (isDragging) {
    trackball->track(Eigen::Vector2i(xpos, ypos));
  } else if (isPanning) {
    float dx = static_cast<float>(xpos - lastMouseX);
    float dy = static_cast<float>(ypos - lastMouseY);
    constexpr float PAN_FACTOR = 0.002f;
    auto camera = scene->getCamera();
    camera->localTranslate(Eigen::Vector3f(-dx * PAN_FACTOR, dy * PAN_FACTOR, 0));
    lastMouseX = xpos;
    lastMouseY = ypos;
  }
}

void Controller::handleKeyPress(int key) {
  auto camera = scene->getCamera();
  if (key == GLFW_KEY_W) {
    camera->zoom(ZOOM_SPEED);
  } else if (key == GLFW_KEY_S) {
    camera->zoom(-ZOOM_SPEED);
  } else if (key == GLFW_KEY_A) {
    camera->localTranslate(Eigen::Vector3f(-PAN_SPEED, 0, 0));
  } else if (key == GLFW_KEY_D) {
    camera->localTranslate(Eigen::Vector3f(PAN_SPEED, 0, 0));
  } else if (key == GLFW_KEY_Z) {
    camera->localTranslate(Eigen::Vector3f(0, PAN_SPEED, 0));
  } else if (key == GLFW_KEY_X) {
    camera->localTranslate(Eigen::Vector3f(0, -PAN_SPEED, 0));
  }
}

void Controller::handleMouseScroll(double xoffset, double yoffset) {
  auto camera = scene->getCamera();
  if (camera) {
    camera->zoom((float)yoffset * ZOOM_SPEED);
  }
}

void Controller::resetControl()
{
    isComputed = false;
    colorMode = 0;
    colorMax = 20.0f;
    selectionMode = (int)SelectType::None;
}

void Controller::render() {
  Polyhedron* mesh = scene->getModel()->getPolyhedron();
  MeshRenderer *meshRenderer = scene->getModel()->getMeshRenderer();
  // Render GUI
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  ImGui::Begin("Control Panel", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  //
  ImGui::Text("Selected Model");
  const std::string& str = scene->modelNames[scene->selectedModel];
  if (ImGui::BeginCombo("##Model", str.c_str())) {
      for (int n = 0; n < scene->modelNames.size(); n++) {
          bool is_selected = (scene->selectedModel == n);
          if (ImGui::Selectable(scene->modelNames[n].c_str(), is_selected)) {
              scene->loadModel(n);
              resetControl();
          }
          if (is_selected) {
              ImGui::SetItemDefaultFocus();
          }
      }
      ImGui::EndCombo();
  }
  ImGui::Separator();
  // Misc controls
  if (ImGui::CollapsingHeader("Misc.")) {
    if (ImGui::Button("Screenshot"))
      saveScreenshot();
    if (ImGui::Button("Save Camera Position"))
      exportCameraPosition();
    if (ImGui::Button("Load Camera Position"))
      importCameraPosition();
    ImGui::Separator();
  }
  // Visualization controls
  if (ImGui::CollapsingHeader("Visualization")) {
    ImGui::Text("Color Mode");
    if (ImGui::Combo("##Color Mode", &colorMode, "None\0Valance Deficit\0Angle Deficit\0Area\0Gauss Curvature\0Mean Curvature\0Max Principle Curvature\0Min Principle Curvature")) {
        if (isComputed)
        {
            meshRenderer->setColors(mesh, colorMode, (double)colorMax);
        }
        else if(colorMode > 3)
        {
            colorMode = 0;
        }
    }
    // ColorMap
    {
        ImGui::SliderFloat(" ", &colorMax, 0.0f, 100.0f);
        meshRenderer->setColors(mesh, colorMode, (double)colorMax);
    }
    if (ImGui::Checkbox("Flat Shading", &scene->drawFaceNormal)) {
      meshRenderer->setNormalMode(scene->getModel()->getPolyhedron(), scene->drawFaceNormal);
    }
    ImGui::Checkbox("Wireframe", &scene->drawWireframe);
    ImGui::Checkbox("Vertex Label", &scene->drawVertexLabel);
    ImGui::Checkbox("Vertex Normal", &scene->drawVertexNormal);
    if (ImGui::Checkbox("Principal Direction", &scene->drawPrincipalDirection)) {
        if (!isComputed) {
            scene->drawPrincipalDirection = false;
        }
    }
    ImGui::Separator();
  }

  // Selection controls
  if (ImGui::CollapsingHeader("Selection")) {
    ImGui::Text("Selection Mode");
    ImGui::Combo("##Selection Mode", &selectionMode, "None\0Faces\0Vertices\0Edges\0");
    if (ImGui::Button("Print Selected Vertex Info"))
        printVertGeometryInfo();
    if (ImGui::Button("Print Selected Edge Info"))
        printEdgeGeometryInfo();
    if (ImGui::Button("Print Selected Face Info"))
        printFaceGeometryInfo();
    if (ImGui::Button("Clear Selections"))
        clearSelections();
  }

  ImGui::Separator();
  if (ImGui::Button("Detect Non-manifold Vertics"))
      selectNonManifoldVertex();
  if (ImGui::Button("Detect Non-manifold Edges"))
      selectNonManifoldEdge();
  if (ImGui::Button("Detect Holes"))
      detectHoles();
  if (ImGui::Button("Compute Global Geometry"))
      computeGlobalInfo();
  if (ImGui::Button("Compute Curvature"))
      computeCurvature();

  ImGui::End();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Controller::clearSelections() {
  for (auto &v : scene->getModel()->getPolyhedron()->vlist) {
    v->selected = false;
  }
  for (auto &e : scene->getModel()->getPolyhedron()->elist) {
    e->selected = false;
  }
  for (auto &t : scene->getModel()->getPolyhedron()->tlist) {
    t->selected = false;
  }
  scene->getModel()->getMeshRenderer()->updateColors(scene->getModel()->getPolyhedron());
}

Eigen::Vector3f Controller::calculateNearPlanePoint(const double &x, const double &y, int height,
                                                    Camera *camera) {
  float near = camera->nearDist();
  Eigen::Vector2f point2D(x, height - y);
  Eigen::Vector3f pt = camera->unProject(point2D, near);

  return pt;
}
void Controller::saveScreenshot() {
  std::vector<unsigned char> pixels(screenWidth * screenHeight * 3);
  glReadPixels(0, 0, screenWidth, screenHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

  auto now = std::chrono::system_clock::now();
  auto time = std::chrono::system_clock::to_time_t(now);
  std::stringstream ss;
  ss << "screenshot_" << std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S") << ".ppm";

  std::ofstream file(ss.str(), std::ios::binary);
  if (!file) {
    std::cerr << "Failed to open file for writing: " << ss.str() << std::endl;
    return;
  }

  file << "P6\n" << screenWidth << " " << screenHeight << "\n255\n";

  for (int y = screenHeight - 1; y >= 0; y--) {
    for (int x = 0; x < screenWidth; x++) {
      int idx = (y * screenWidth + x) * 3;
      file.write(reinterpret_cast<char *>(&pixels[idx]), 3);
    }
  }

  file.close();
  std::cout << "Screenshot saved to: " << ss.str() << std::endl;
}

void Controller::exportCameraPosition(const std::string &filename) {
  if (!scene->getCamera()) {
    std::cerr << "Camera not initialized" << std::endl;
    return;
  }

  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  Eigen::Vector3f position = scene->getCamera()->position();
  Eigen::Vector3f target = scene->getCamera()->target();
  float fovy = scene->getCamera()->fovY();

  file << "[Camera]\n";
  file << "position_x=" << position.x() << "\n";
  file << "position_y=" << position.y() << "\n";
  file << "position_z=" << position.z() << "\n";
  file << "target_x=" << target.x() << "\n";
  file << "target_y=" << target.y() << "\n";
  file << "target_z=" << target.z() << "\n";
  file << "fovy=" << fovy << "\n";

  file.close();
  std::cout << "Camera position saved to: " << filename << std::endl;
}

void Controller::importCameraPosition(const std::string &filename) {
  if (!scene->getCamera()) {
    std::cerr << "Camera not initialized" << std::endl;
    return;
  }

  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open file for reading: " << filename << std::endl;
    return;
  }

  scene->getCamera()->reset();

  Eigen::Vector3f position, target;
  float fovy = 60.f;
  std::string line;
  std::string section;

  while (std::getline(file, line)) {
    // Trim whitespace
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);

    if (line.empty() || line[0] == ';')
      continue;

    if (line[0] == '[' && line[line.length() - 1] == ']') {
      section = line.substr(1, line.length() - 2);
      continue;
    }

    if (section == "Camera") {
      size_t delimPos = line.find('=');
      if (delimPos != std::string::npos) {
        std::string key = line.substr(0, delimPos);
        float value = std::stof(line.substr(delimPos + 1));

        if (key == "position_x")
          position.x() = value;
        else if (key == "position_y")
          position.y() = value;
        else if (key == "position_z")
          position.z() = value;
        else if (key == "target_x")
          target.x() = value;
        else if (key == "target_y")
          target.y() = value;
        else if (key == "target_z")
          target.z() = value;
        else if (key == "fovy")
          fovy = value;
      }
    }
  }

  scene->getCamera()->setPosition(position);
  scene->getCamera()->setTarget(target);
  scene->getCamera()->setFovY(fovy);

  file.close();
  std::cout << "Camera position loaded from: " << filename << std::endl;
}

int Controller::selectFaces(const double &x, const double &y) {
  // Convert screen coordinates to a point on the near plane
  Eigen::Vector3f pt = calculateNearPlanePoint(x, y, screenHeight, scene->getCamera());

  // Get mesh data structures
  Polyhedron *mesh = scene->getModel()->getPolyhedron();
  const std::vector<Triangle *> &triangles = mesh->tlist;
  const std::vector<Vertex *> &vertices = mesh->vlist;

  int tidx = -1;
  // Setup ray origin and direction for intersection test
  Eigen::Vector3f cam_org = scene->getCamera()->position();
  Eigen::Vector3f cam_dir = (pt - cam_org).normalized();

  // Initialize target with maximum values
  Eigen::Vector3f target(FLT_MAX, FLT_MAX, FLT_MAX); // u,v,t
  // Get the starting time point
  auto start = std::chrono::high_resolution_clock::now();
  // Test intersection with each triangle
  for (int i = 0; i < triangles.size(); i++) {
    const Triangle *tri = triangles[i];
    Eigen::Vector3f p0 = tri->verts[0]->pos.cast<float>();
    Eigen::Vector3f p1 = tri->verts[1]->pos.cast<float>();
    Eigen::Vector3f p2 = tri->verts[2]->pos.cast<float>();
    Eigen::Vector3f bary;
    // Check if ray intersects current triangle
    if (MeshProcessor::rayIntersectsTriangle(cam_org, cam_dir, p0, p1, p2, bary)) {
      std::cout << bary.transpose() << std::endl;
      // Update if this intersection is closer than previous ones
      if (bary.z() < target.z()) {
        target = bary;
        tidx = i;
      }
    }
  }
  // Get the ending time point
  auto end = std::chrono::high_resolution_clock::now();

  // Toggle selection state of the intersected triangle if found
  if (tidx != -1) {
    mesh->tlist[tidx]->selected = !mesh->tlist[tidx]->selected;
    scene->getModel()->getMeshRenderer()->updateColors(mesh);
  }
  std::cout << "Select Triangle Index: " << tidx;
  // Calculate the elapsed time in milliseconds
  std::chrono::duration<double, std::milli> elapsed = end - start;
  std::cout << " (Elapsed time: " << elapsed.count() << " ms)" << std::endl;
  return tidx;
}

int Controller::selectVertices(const double &x, const double &y) {
  Eigen::Vector3f pt = calculateNearPlanePoint(x, y, screenHeight, scene->getCamera());
  Polyhedron *mesh = scene->getModel()->getPolyhedron();
  const std::vector<Triangle *> &triangles = mesh->tlist;
  const std::vector<Vertex *> &vertices = mesh->vlist;
  int tidx = -1;
  Eigen::Vector3f cam_org = scene->getCamera()->position();
  Eigen::Vector3f cam_dir = (pt - cam_org).normalized();
  Eigen::Vector3f target(FLT_MAX, FLT_MAX, FLT_MAX); // u,v,t
  for (int i = 0; i < triangles.size(); i++) {
    const Triangle *tri = triangles[i];
    Eigen::Vector3f p0 = tri->verts[0]->pos.cast<float>();
    Eigen::Vector3f p1 = tri->verts[1]->pos.cast<float>();
    Eigen::Vector3f p2 = tri->verts[2]->pos.cast<float>();
    Eigen::Vector3f bary;

    if (MeshProcessor::rayIntersectsTriangle(cam_org, cam_dir, p0, p1, p2, bary)) {
      if (bary.z() < target.z()) {
        target = bary;
        tidx = i;
      }
    }
  }

  int vidx = -1;
  if (tidx != -1) {
    Triangle *tri = triangles[tidx];
    if (target.x() <= 0.5 && target.y() <= 0.5) {
      vidx = tri->verts[0]->index;
    } else if (target.x() <= 0.5) {
      vidx = tri->verts[2]->index;
    } else if (target.y() <= 0.5) {
      vidx = tri->verts[1]->index;
    }
    if (vidx != -1) {
      mesh->vlist[vidx]->selected = !mesh->vlist[vidx]->selected;
    }
  }
  std::cout << "Select Vertex Index: " << vidx << '\n';
  return vidx;
}

int Controller::selectEdges(const double &x, const double &y) {
  Eigen::Vector3f pt = calculateNearPlanePoint(x, y, screenHeight, scene->getCamera());
  Polyhedron *mesh = scene->getModel()->getPolyhedron();
  const std::vector<Triangle *> &triangles = mesh->tlist;
  const std::vector<Vertex *> &vertices = mesh->vlist;
  int tidx = -1;
  Eigen::Vector3f cam_org = scene->getCamera()->position();
  Eigen::Vector3f cam_dir = (pt - cam_org).normalized();
  Eigen::Vector3f target(FLT_MAX, FLT_MAX, FLT_MAX); // u,v,t
  for (int i = 0; i < triangles.size(); i++) {
    const Triangle *tri = triangles[i];
    Eigen::Vector3f p0 = tri->verts[0]->pos.cast<float>();
    Eigen::Vector3f p1 = tri->verts[1]->pos.cast<float>();
    Eigen::Vector3f p2 = tri->verts[2]->pos.cast<float>();
    Eigen::Vector3f bary;
    if (MeshProcessor::rayIntersectsTriangle(cam_org, cam_dir, p0, p1, p2, bary)) {
      if (bary.z() < target.z()) {
        target = bary;
        tidx = i;
      }
    }
  }

  int eidx = -1;
  if (tidx != -1) {
    Triangle* tri = triangles[tidx];
    float area1 = target.y();
    float area2 = 1 - target.x() - target.y();
    float area3 = target.x();
    if (area1 <= area2 && area1 <= area3) {
        eidx = tri->edges[0]->index;
    }
    else if (area2 <= area1 && area2 <= area3) {
        eidx = tri->edges[1]->index;
    }
    else if (area3 <= area1 && area3 <= area2) {
        eidx = tri->edges[2]->index;
    }

    if (eidx != -1) {
      mesh->elist[eidx]->selected = !mesh->elist[eidx]->selected;
    }
  }
  std::cout << "Select Edge Index: " << eidx << '\n';
  return eidx;
}

void Controller::selectNonManifoldEdge()
{
    clearSelections();
    Polyhedron* mesh = scene->getModel()->getPolyhedron();
    for (int i = 0; i < mesh->nedges(); i++)
    {
        if (MeshProcessor::isNonManifoldEdge(mesh->elist[i]))
        {
            mesh->elist[i]->selected = true;
        }
    }
}

void Controller::selectNonManifoldVertex()
{
    clearSelections();
    Polyhedron* mesh = scene->getModel()->getPolyhedron();
    for (int i = 0; i < mesh->nverts(); i++)
    {
        if (MeshProcessor::isNonManifoldVert(mesh->vlist[i]))
        {
            mesh->vlist[i]->selected = true;
        }
    }
}

void Controller::detectHoles()
{
    clearSelections();
    Polyhedron* mesh = scene->getModel()->getPolyhedron();
    std::vector<std::vector<int>> holes;
    MeshProcessor::findHoles(mesh, holes);
    for (int i = 0; i < holes.size(); i++)
    {
        for (int j = 0; j < holes[i].size(); j++)
        {
            Edge* edge = mesh->elist[holes[i][j]];
            edge->selected = true;
            edge->verts[0]->selected = true;
            edge->verts[1]->selected = true;
        }
    }
}

void Controller::computeGlobalInfo()
{
    Polyhedron* mesh = scene->getModel()->getPolyhedron();
    std::cout << "=== Global Geometry ===" << std::endl;
    // Get the starting time point
    auto start = std::chrono::high_resolution_clock::now();
    int euler = MeshProcessor::calcEulerCharacteristic(mesh);
    double deficit = MeshProcessor::calcTotalAngleDeficit(mesh);
    double volumn = MeshProcessor::calcVolume(mesh);
    double area = MeshProcessor::calcTotalFaceArea(mesh);
    int valence = MeshProcessor::calcTotalValenceDeficit(mesh);
    // Get the ending time point
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the elapsed time in milliseconds
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << " #V: " << mesh->nverts() << std::endl;
    std::cout << " #E: " << mesh->nedges() << std::endl;
    std::cout << " #F: " << mesh->ntris() << std::endl;
    std::cout << " Euler Characteristic: " << euler << std::endl;
    std::cout << " Total Angular Deficit: " << deficit << std::endl;
    std::cout << " Total Valence: " << valence << std::endl;
    std::cout << " Volume: " << volumn << std::endl;
    std::cout << " Area: " << area << std::endl;
    std::cout << " Elapsed Time: " << elapsed.count() << " ms" << std::endl;
    std::cout << "=======================" << std::endl;
}

void Controller::computeCurvature()
{
   Polyhedron* mesh = scene->getModel()->getPolyhedron();
   std::cout << "====== Curvature ======" << std::endl;
   auto start = std::chrono::high_resolution_clock::now();
   // 
   std::cout << " Gauss Curvature..." << std::endl;
   MeshProcessor::calcGaussCurvature(mesh);
   std::cout << " Mean Curvature..." << std::endl;
   MeshProcessor::calcMeanCurvature(mesh);
   std::cout << " Principal Curvature..." << std::endl;
   MeshProcessor::calcPrincipalCurvature(mesh);
   std::cout << " Principal Tensor..." << std::endl;
   MeshProcessor::calcCurvatureTensor(mesh);
   // Get the ending time point
   auto end = std::chrono::high_resolution_clock::now();
   // Calculate the elapsed time in milliseconds
   std::chrono::duration<double, std::milli> elapsed = end - start;
   std::cout << " Elapsed Time: " << elapsed.count() << " ms" << std::endl;
   std::cout << "=======================" << std::endl;
   isComputed = true;
}

void Controller::printVertGeometryInfo() {
  
  std::cout << "===== Vertex Info =====" << std::endl;
  for (int i = 0; i < scene->getModel()->getPolyhedron()->nverts(); i++) {
    Vertex *vert_i = scene->getModel()->getPolyhedron()->vlist[i];
    if (vert_i->selected) {
      std::cout << "Vertex ID: " << vert_i->index << std::endl;
      std::cout << " Position: (" << vert_i->pos[0] << ", " << vert_i->pos[1] << ", " << vert_i->pos[2] << ")" << std::endl;
      std::cout << " Vertex Normal: (" << vert_i->normal[0] << ", " << vert_i->normal[1] << ", " << vert_i->normal[2] << ")" << std::endl;
      std::cout << " Vertex Area: " << vert_i->area << std::endl;
      //
      if (isComputed) 
      { 
          double deficit = MeshProcessor::calcAngleDeficit(vert_i);
          std::cout << " Angular Deficit: " << deficit << std::endl;
          std::cout << " Mean Curvature: " << vert_i->meanCurvature << std::endl;
          std::cout << " Gauss Curvature: " << vert_i->gaussCurvature << std::endl;
          std::cout << " Max Principal Curvature: " << vert_i->maxPrincCurvature << std::endl;
          std::cout << " Min Principal Curvature: " << vert_i->minPrincCurvature << std::endl;
          //
          std::cout << " princDir[0]: (" << vert_i->princDir3D[0][0] << ", " << vert_i->princDir3D[0][1]
                    << ", " << vert_i->princDir3D[0][2] << ")" << std::endl;
          std::cout << " princDir[1]: (" << vert_i->princDir3D[1][0] << ", " << vert_i->princDir3D[1][1]
                    << ", " << vert_i->princDir3D[1][2] << ")" << std::endl;
      }
    }
  }
  std::cout << "=======================" << std::endl;
}

void Controller::printEdgeGeometryInfo() {
    std::cout << "===== Edge Info =====" << std::endl;
    std::cout << ":" << std::endl;
    for (int i = 0; i < scene->getModel()->getPolyhedron()->nedges(); i++) {
        Edge* edge_i = scene->getModel()->getPolyhedron()->elist[i];
        if (edge_i->selected) {
            std::cout << "Edge ID: " << edge_i->index << std::endl;
            std::cout << "Length: " << edge_i->length << std::endl;
            std::cout << "Dihedral Angle: " << edge_i->dihedral_angle << std::endl;
        }
    }
    std::cout << "=====================" << std::endl;
}

void Controller::printFaceGeometryInfo()
{
    std::cout << "===== Face Info =====" << std::endl;
    std::cout << "Face Info:" << std::endl;
    for (int i = 0; i < scene->getModel()->getPolyhedron()->ntris(); i++) {
        Triangle* tri_i = scene->getModel()->getPolyhedron()->tlist[i];
        if (tri_i->selected) {
            std::cout << "Face ID: " << tri_i->index << std::endl;
            std::cout << "Face Normal: (" << tri_i->normal[0] << ", " << tri_i->normal[1] << ", " << tri_i->normal[2] << ")" << std::endl;
            std::cout << "Area: " << tri_i->area << std::endl;
        }
    }
    std::cout << "=====================" << std::endl;
}