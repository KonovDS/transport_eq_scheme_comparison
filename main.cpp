#include <iostream>
#include <vector>
#include <string_view>
#include <fstream>

struct Mesh {
  std::vector<double> u;

  explicit Mesh(size_t size) : u(size + 1) {
  }

  void Load(const std::string &path) {
    std::ifstream str(path);
    for (auto &x : u) {
      str >> x;
    }
  }

  void Save(const std::string &path) {
    std::ofstream str(path);
    for (auto &x : u) {
      str << x << " ";
    }
  }

  void Heaviside() {
    size_t i = 0;
    for (; i < u.size()/2; ++i) {
      u[i] = 1;
    }
    for (; i < u.size(); ++i) {
      u[i] = 0;
    }
  }

  void BoundaryHeaviside() {
    u[0] = 1;
    u[u.size() - 1] = 0;
  }
};

void ExplicitGodunovScheme(const Mesh &b, Mesh &e, double c) {
  for (size_t i = 1; i < b.u.size() - 1; ++i) {
    e.u[i] = b.u[i] - c * (b.u[i] - b.u[i - 1]);
  }
}

void ExplicitMacCormackScheme(const Mesh &b, Mesh &e, double c) {
  for (size_t i = 1; i < b.u.size() - 1; ++i) {
    e.u[i] = b.u[i] - 0.5 * c * (b.u[i + 1] - b.u[i - 1]) + 0.5 * c * c * (b.u[i + 1] - 2 * b.u[i] +  b.u[i - 1]);
  }
}

void ExplicitHolodnovScheme(const Mesh &b, Mesh &e, double c) {
  // + Дополнение граничного условия слева (по схеме Мак-Кормака)
  size_t i = 1;
  e.u[i] = b.u[i] - 0.5 * c * (b.u[i + 1] - b.u[i - 1]) + 0.5 * c * c * (b.u[i + 1] - 2 * b.u[i] +  b.u[i - 1]);

  for (++i; i < b.u.size() - 1; ++i) {
    e.u[i] = b.u[i] - c * (b.u[i] - b.u[i - 1]) - 0.25 * c * (1 - c) * (b.u[i + 1] - b.u[i] - b.u[i - 1]+  b.u[i - 2]);
  }
}

void CubeApproxScheme(const Mesh &b, Mesh &e, double c) {
  // + Дополнение граничного условия слева (по схеме Мак-Кормака)
  size_t i = 1;
  e.u[i] = b.u[i] - 0.5 * c * (b.u[i + 1] - b.u[i - 1]) + 0.5 * c * c * (b.u[i + 1] - 2 * b.u[i] +  b.u[i - 1]);

  for (++i; i < b.u.size() - 1; ++i) {
    e.u[i] = b.u[i] - c * (b.u[i] - b.u[i - 1]) - c / 6 * (2 - c) * (1 - c) * b.u[i + 1] +
        c / 2 * (1 - c) * (1 - c) * b.u[i] + c * c / 2 * (1 - c) * b.u[i - 1] + c / 6 * (c * c - 1) * b.u[i - 2];
  }
}

void HybridScheme(const Mesh &b, Mesh &e, double c) {
  // + Дополнение граничного условия слева (по схеме Годунова)
  size_t i = 1;
  e.u[i] = b.u[i] - c * (b.u[i] - b.u[i - 1]);
  for (++i; i < b.u.size() - 1; ++i) {
    // u > 1
    double delta1 = b.u[i - 1] - b.u[i];
    double delta2 = b.u[i - 2] - b.u[i - 1];
    if(delta1 * (delta1 - delta2) <= 0) {
      e.u[i] = -0.5 * c * (1 - c) * b.u[i + 1] + (1 - c*c) * b.u[i] + 0.5 * c * (1 + c) * b.u[i - 1];
    } else {
      e.u[i] = (1 - 0.5 * c * (3 - c)) * b.u[i] + c * (2 - c) * b.u[i - 1] - 0.5 * c * (1 - c) * b.u[i - 2];
    }
  }
}

double AutoCourant(double h, double tau, double u) {
  return std::abs(u) * tau / h;
}

int main() {
  Mesh b(100), f(100);
  b.Heaviside();

  double h = 1;
  double tau = 0.1;
  double u = 1;
  double c = AutoCourant(h, tau, u);

  size_t steps = 100;

  if (c > 1) {
    std::cout << "Courant condition is not met! Schemes are unstable.";
  }
  b.Save("../0.txt");

  for (size_t i = 0; i < steps; ++i) {
    ExplicitGodunovScheme(b, f, c);
    f.BoundaryHeaviside();
    std::swap(b, f);
  }
  f.Save("../1.txt");

  b.Heaviside();
  for (size_t i = 0; i < steps; ++i) {
    ExplicitMacCormackScheme(b, f, c);
    f.BoundaryHeaviside();
    std::swap(b, f);
  }
  f.Save("../2.txt");

  b.Heaviside();
  for (size_t i = 0; i < steps; ++i) {
    ExplicitHolodnovScheme(b, f, c);
    f.BoundaryHeaviside();
    std::swap(b, f);
  }
  f.Save("../3.txt");

  b.Heaviside();
  for (size_t i = 0; i < steps; ++i) {
    CubeApproxScheme(b, f, c);
    f.BoundaryHeaviside();
    std::swap(b, f);
  }
  f.Save("../4.txt");

  b.Heaviside();
  for (size_t i = 0; i < steps; ++i) {
    HybridScheme(b, f, c);
    f.BoundaryHeaviside();
    std::swap(b, f);
  }
  f.Save("../5.txt");

  return 0;
}