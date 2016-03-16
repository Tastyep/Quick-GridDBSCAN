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

void DBSCAN::expandCluster(std::vector<std::reference_wrapper<Cell>> &cluster,
                           const Cell &cell, bool newCluster) {
  static constexpr int pos[4][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}};
  int clusterId = cell.clusterId;

  std::cout << "pass " << newCluster << "\n"
            << "\n";
  for (unsigned int i = 0; i < 4; ++i) {
    int idX = cell.x + pos[i][0];
    int idY = cell.y + pos[i][1];

    std::cout << idX << " " << idY << " | " << this->gridWidth << " "
              << this->gridHeight << "\n";
    if (idX >= this->gridWidth || idX < 0 || idY >= this->gridHeight || idY < 0)
      continue; // Out of bound
    Cell &nextCell = this->cells[idY][idX];

    std::cout << "clusterId: " << clusterId << " " << nextCell.clusterId << " "
              << " newCluster: " << newCluster << " " << i << "\n";
    if (nextCell.points.empty() ||
        (clusterId == nextCell.clusterId && !newCluster))
      continue;
    if (nextCell.clusterId != -1) { // If is already part of a cluster
      auto &nextCluster = this->clusters[nextCell.clusterId];
      if (newCluster) { // Merge our new cluster into the other one
        std::cout << "NewCluster"
                  << "\n";
        nextCluster.insert(nextCluster.end(), cluster.begin(), cluster.end());
        for (auto &cellRef : cluster)
          cellRef.get().clusterId = nextCell.clusterId;
        newCluster = false;
        cluster =
            nextCluster; // change the reference to operate on the new cluster
        clusterId = nextCell.clusterId;
      } else {
        // If not a newCluster, get the biggest one
        // clusterM = clusterMaster
        // clusterS = clusterSlave
        std::cout << "not NewCluster"
                  << "\n";

        auto &clusterM =
            (cluster.size() > nextCluster.size() ? cluster : nextCluster);
        auto &clusterS =
            (cluster.size() > nextCluster.size() ? nextCluster : cluster);
        int clusterMid = clusterM.front().get().clusterId;
        int clusterSid = clusterS.front().get().clusterId;

        clusterM.insert(clusterM.end(), clusterS.begin(), clusterS.end());
        for (auto &cellRef : clusterS)
          cellRef.get().clusterId = clusterMid;

        std::cout << this->clusters.size() << " " << clusterSid << " "
                  << clusterMid << " newCluster: " << newCluster << "\n";
        this->toDelete.push_back(clusterSid);
        clusterId = clusterMid;
        cluster = clusterM;
      }
    } else { // If nextCell is not part of a cluster
      cluster.emplace_back(nextCell);
      if (newCluster) {
        clusterId = this->clusters.size();
        this->clusters.push_back(cluster);

        cluster = this->clusters.back();
        for (auto &refCell : cluster)
          refCell.get().clusterId = clusterId;
        newCluster = false;
      } else
        nextCell.clusterId = clusterId;
    }
  }
  if (newCluster) {
    cluster.front().get().clusterId = this->clusters.size();
    this->clusters.push_back(cluster);
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
        this->clusters.push_back(std::vector<std::reference_wrapper<Cell>>());
        // Create a new cluster and add my two elements
        auto &cluster = this->clusters.back();
        cluster.emplace_back(cell);
        cluster.emplace_back(nextCell);
      } else { // If nextCell is part of a cluster
        this->clusters[nextCell.clusterId].emplace_back(cell);
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
        std::vector<std::reference_wrapper<Cell>> newCluster;

        newCluster.emplace_back(cell);
        this->expandCluster(newCluster, cell, true);
      } else
        this->expandCluster(this->clusters[cell.clusterId], cell, false);
    } else {
      if (cell.clusterId != -1) // Part of a cluster but not a core cell
        continue;               // non core cell can't expand
      else
        this->reachCluster(cell);
    }
  }
}

void DBSCAN::simplifyClusters() {
  for (unsigned int i = 0; i < this->clusters.size(); ++i) {
    if (std::find(this->toDelete.begin(), this->toDelete.end(), i) !=
        this->toDelete.end())
      this->clusters.erase(this->clusters.begin() + i);
  }
  this->toDelete.clear();
  for (auto &cluster : this->clusters) {
    this->sclusters.push_back(std::vector<Point>());

    auto &scluster = this->sclusters.back();
    for (auto &cellRef : cluster) {
      scluster.insert(scluster.end(), cellRef.get().points.begin(),
                      cellRef.get().points.end());
    }
  }
}

void DBSCAN::cluster(const std::vector<Point> &points) {
  this->clusters.clear();
  this->sclusters.clear();
  this->noise.clear();
  this->cells.clear();
  this->filledCells.clear();
  constructGrid(points);
  constructClusters();
  simplifyClusters();
}

const std::vector<Point> &DBSCAN::getNoise() const { return this->noise; }

const std::vector<std::vector<Point>> &DBSCAN::getClusters() const {
  return this->sclusters;
}
}