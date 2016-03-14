#include "DBSCAN.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace Cluster {

DBSCAN::DBSCAN(unsigned int pts, float dist) : pts(pts) {
  this->squareSize = dist / std::sqrt(2);
}

void DBSCAN::constructGrid(const std::vector<Point> &points) {
  Point min = points.front();
  Point max = points.front();
  unsigned int x = 0;
  unsigned int y = 0;

  for (const auto &p : points) {
    if (p.x < min.x)
      min.x = p.x;
    else if (p.x > max.x)
      max.x = p.x;
    if (p.y < min.y)
      min.y = p.y;
    else if (p.y > max.y)
      max.y = p.y;
  }
  this->gridWidth = (max.x - min.x) / this->squareSize + 1;
  this->gridHeight = (max.y - min.y) / this->squareSize + 1;
  this->cells.resize(this->gridHeight);
  for (auto &cell : this->cells)
    cell.resize(this->gridWidth);
  for (const auto &p : points) {
    this->cells[(p.y - min.y) / this->squareSize]
               [(p.x - min.x) / this->squareSize]
                   .points.push_back(p);
  }
  for (auto &cellY : this->cells) {
    x = 0;
    for (auto &cellX : cellY) {
      cellX.x = x;
      cellX.y = y;
      cellX.isCore = (cellX.points.size() >= this->pts);
      if (!cellX.points.empty()) {
        this->filledCells.emplace_back(cellX);
      }
      ++x;
    }
    ++y;
  }
}

void DBSCAN::expandCluster(std::vector<Point> &cluster, const Cell &cell) {
  static constexpr int pos[4][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}};

  for (unsigned int i = 0; i < 4; ++i) {
    int idX = cell.x + pos[i][0];
    int idY = cell.y + pos[i][1];

    if (idX >= this->gridWidth || idX < 0 || idY >= this->gridHeight || idY < 0)
      continue; // Out of bound
    Cell &nextCell = this->cells[idY][idX];

    if (nextCell.points.empty())
      continue;
    else {
      nextCell.clusterId = cell.clusterId;
      cluster.insert(cluster.end(), nextCell.points.begin(),
                     nextCell.points.end());
    }
  }
}

void DBSCAN::reachCluster(Cell &cell) {
  static constexpr int pos[4][2] = {{-1, 1}, {0, 1}, {1, 0}, {1, 1}};
  unsigned int i;

  for (i = 0; i < 4; ++i) {
    int idX = cell.x + pos[i][0];
    int idY = cell.y + pos[i][1];

    if (idX >= this->gridWidth || idX < 0 || idY >= this->gridHeight || idY < 0)
      continue; // Out of bound
    Cell &nextCell = this->cells[idY][idX];

    if (nextCell.isCore) {
      if (nextCell.clusterId == -1) {
        cell.clusterId = this->clusters.size();
        nextCell.clusterId = cell.clusterId;
        this->clusters.push_back(cell.points);
        auto &cluster = this->clusters.back();
        cluster.insert(cluster.end(), nextCell.points.begin(),
                       nextCell.points.end());
      } else {
        this->clusters[nextCell.clusterId].insert(
            this->clusters[nextCell.clusterId].end(), cell.points.begin(),
            cell.points.end());
        cell.clusterId = nextCell.clusterId;
      }
      break;
    }
  }
  if (i == 4) {
    this->noise.insert(this->noise.end(), cell.points.begin(),
                       cell.points.end());
  }
}

void DBSCAN::constructClusters() {
  for (auto &refCell : this->filledCells) {
    auto &cell = refCell.get();

    if (cell.isCore) { // pass a var to merge my core cell if it is not yet in a
                       // cluster into the newly(future) found cluster
      if (cell.clusterId == -1) {
        this->clusters.push_back(cell.points);
        cell.clusterId = this->clusters.size() - 1;
        this->expandCluster(this->clusters.back(), cell);
      } else
        this->expandCluster(this->clusters[cell.clusterId], cell);
    } else {
      if (cell.clusterId != -1)
        continue; // non core cell can't expand
      else
        this->reachCluster(cell);
    }
  }
}

void DBSCAN::cluster(const std::vector<Point> &points) {
  this->clusters.clear();
  this->noise.clear();
  this->cells.clear();
  this->filledCells.clear();
  constructGrid(points);
  constructClusters();
}

const std::vector<Point> &DBSCAN::getNoise() const { return this->noise; }

const std::vector<std::vector<Point>> &DBSCAN::getClusters() const {
  return this->clusters;
}
}