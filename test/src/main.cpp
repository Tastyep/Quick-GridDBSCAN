#include "Render.hh"

#define WIDTH 1600
#define HEIGHT 900

int main(int argc, char *argv[]) {
  Render render(WIDTH, HEIGHT);

  render.run();
  return 0;
}
