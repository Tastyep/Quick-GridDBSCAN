#include "ConvexHull.hh"

#include <algorithm>
#include <iostream>

ConvexHull::ConvexHull(const std::vector<Point> &points) : points(points) {
  std::sort(this->points.begin(), this->points.end());
}

void ConvexHull::assign(const std::vector<Point> &points) {
  this->points = points;
  std::sort(this->points.begin(), this->points.end());
}

int ConvexHull::cross(const Point &O, const Point &A, const Point &B) const {
  return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

const std::vector<Point> &ConvexHull::getHull() const { return this->hull; }

const std::vector<Point> &ConvexHull::computeHull() {
  this->hull.clear();

  std::cout << "Nb points: " << this->points.size() << "\n";
  // Build lower hull
  std::vector<Point> lowerHull;
  for (const auto &p : this->points) {
    int hullSize = lowerHull.size();
    while (hullSize >= 2 &&
           cross(lowerHull[hullSize - 2], lowerHull[hullSize - 1], p) <= 0) {
      lowerHull.pop_back();
      hullSize = lowerHull.size();
    }
    lowerHull.push_back(p);
  }

  // Build upper hull
  std::vector<Point> upperHull;
  for (int i = this->points.size() - 1; i >= 0; --i) {
    int hullSize = upperHull.size();
    while (hullSize >= 2 &&
           cross(upperHull[hullSize - 2], upperHull[hullSize - 1],
                 this->points[i]) <= 0) {
      upperHull.pop_back();
      hullSize = upperHull.size();
    }
    upperHull.push_back(this->points[i]);
  }

  lowerHull.pop_back(); // Already contained in the other hull
  std::cout << "lowerHull Size: " << lowerHull.size() << "\n";
  this->hull = lowerHull;
  this->hull.insert(this->hull.end(), upperHull.begin(), upperHull.end());
  return this->hull;
}