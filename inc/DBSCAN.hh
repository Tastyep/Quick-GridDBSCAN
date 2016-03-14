#ifndef DBSCAN_HH_
#define DBSCAN_HH_

#include <functional>
#include <vector>

namespace Cluster {

struct Point {
  int x;
  int y;

  Point(int x, int y) : x(x), y(y) {}
};

struct Cell {
  std::vector<Point> points;
  int x;
  int y;
  bool isCore;
  int clusterId;

  Cell() : clusterId(-1) {}
};

class DBSCAN {
public:
  DBSCAN(unsigned int pts, float dist);

  ~DBSCAN() = default;
  DBSCAN(const DBSCAN &other) = default;
  DBSCAN(DBSCAN &&other) = default;
  DBSCAN &operator=(const DBSCAN &other) = default;
  DBSCAN &operator=(DBSCAN &&other) = default;

  void cluster(const std::vector<Point> &points);
  const std::vector<std::vector<Point>> &getClusters() const;
  const std::vector<Point> &getNoise() const;

private:
  void constructGrid(const std::vector<Point> &points);
  void constructClusters();
  void expandCluster(std::vector<Point> &cluster, const Cell &cell);
  void reachCluster(Cell &cell);

private:
  unsigned int pts;
  float squareSize;
  std::vector<std::vector<Point>> clusters;
  std::vector<Point> noise;
  std::vector<std::vector<Cell>> cells;
  std::vector<std::reference_wrapper<Cell>> filledCells;
  int gridWidth;
  int gridHeight;
};
}

#endif /* end of include guard: DBSCAN_HH_ */
