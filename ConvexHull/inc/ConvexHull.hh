#ifndef CONVEXHULL_HH_
#define CONVEXHULL_HH_

#include <vector>

struct Point {
  double x;
  double y;

  Point(double x, double y) : x(x), y(y) {}
  Point(const Point &p) {
    x = p.x;
    y = p.y;
  }

  bool operator<(const Point &p) const {
    return (x < p.x || (x == p.x && y < p.y));
  }
};

class ConvexHull {
public:
  ConvexHull(const std::vector<Point> &points);
  ConvexHull() = default;

  ~ConvexHull() = default;
  ConvexHull(const ConvexHull &other) = default;
  ConvexHull(ConvexHull &&other) = default;
  ConvexHull &operator=(const ConvexHull &other) = default;
  ConvexHull &operator=(ConvexHull &&other) = default;

  void assign(const std::vector<Point> &points);
  const std::vector<Point> &computeHull();
  const std::vector<Point> &getHull() const;

private:
  int cross(const Point &O, const Point &A, const Point &B) const;

private:
  std::vector<Point> points;
  std::vector<Point> hull;
};

#endif /* end of include guard: CONVEXHULL_HH_ */
