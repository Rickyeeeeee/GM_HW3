#pragma once

#include "miscellaneous/camera.h"
#include "miscellaneous/trackball.h"
#include "scene.h"
#include <GLFW/glfw3.h>

class Controller {
public:
  // Constants
  static constexpr float ZOOM_SPEED = 0.1f;
  static constexpr float PAN_SPEED = 0.1f;

  // Enum
  enum class SelectType { None, Faces, Vertices, Edges };

  // Constructor and Destructor
  Controller(int width, int height, GLFWwindow *window, Scene *scenePtr);
  ~Controller() = default;

  // Window Management
  void resize(int width, int height);

  // Rendering Interface
  void render();

  // Input Handling Interface
  void handleMouseClick(int button, double xpos, double ypos, int action);
  void handleMouseMove(double xpos, double ypos);
  void handleKeyPress(int key);
  void handleMouseScroll(double xoffset, double yoffset);

private:
  // Misc functions
  void saveScreenshot();
  void exportCameraPosition(const std::string &filename = "camera.ini");
  void importCameraPosition(const std::string &filename = "camera.ini");

  // Selection functions
  void clearSelections();
  int selectFaces(const double &x, const double &y);
  int selectVertices(const double &x, const double &y);
  int selectEdges(const double &x, const double &y);
  Eigen::Vector3f calculateNearPlanePoint(const double& x, const double& y, int height, Camera* camera);

  // Gobal
  void selectNonManifoldEdge();
  void selectNonManifoldVertex();
  void detectHoles();
  void computeGlobalInfo();
  void computeCurvature();
  void printVertGeometryInfo();
  void printEdgeGeometryInfo();
  void printFaceGeometryInfo();

  void resetControl();

private:
  int screenWidth;
  int screenHeight;
  // State
  bool isDragging;
  bool isPanning;
  bool isComputed;
  double lastMouseX;
  double lastMouseY;
  // Components
  Scene *scene;
  std::unique_ptr<Trackball> trackball;
  int selectionMode;
  int colorMode;
  float colorMax;
};
