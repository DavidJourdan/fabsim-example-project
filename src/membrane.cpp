#include <fsim/ElasticMembrane.h>
#include <fsim/util/io.h>
#include <fsim/util/typedefs.h>
#include <optim/NewtonSolver.h>
#include <polyscope/surface_mesh.h>

int main(int argc, char *argv[]) {
  using namespace Eigen;

  // load geometry from OFF mesh file
  fsim::Mat3<double> V;
  fsim::Mat3<int> F;
  fsim::readOFF("../data/mesh.off", V, F);

  // parameters of the membrane model
  const double young_modulus = 10;
  const double thickness = 0.5;
  const double poisson_ratio = 0.3;
  double stretch_factor = 1.7;
  double mass = 1;

  // declare StVKMembrane object (could be replaced seamlessly with e.g. NeohookeanMembrane)
  fsim::StVKMembrane model(V / stretch_factor, F, thickness, young_modulus, poisson_ratio, mass);

  // declare NewtonSolver object
  optim::NewtonSolver<double> solver;
  // specify fixed degrees of freedom (here the 4 corners of the mesh are fixed)
  solver.options.fixed_dofs = {0 * 3 + 0, 0 * 3 + 1, 0 * 3 + 2, 1 * 3 + 0,
                               1 * 3 + 1, 1 * 3 + 2, 2 * 3 + 0, 2 * 3 + 1,
                               2 * 3 + 2, 3 * 3 + 0, 3 * 3 + 1, 3 * 3 + 2};
  solver.options.threshold = 1e-6; // specify how small the gradient's norm has to be

  // display the mesh
  polyscope::registerSurfaceMesh("mesh", V, F)
      ->setEdgeWidth(1)
      ->setEdgeColor({0.1, 0.1, 0.1});
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::options::groundPlaneHeightFactor = 0.4;
  polyscope::init();

  polyscope::state::userCallback = [&]() 
  {
    ImGui::PushItemWidth(100);
    if(ImGui::InputDouble("Stretch factor", &stretch_factor, 0, 0, "%.1f"))
      model = fsim::StVKMembrane(V / stretch_factor, F, thickness, young_modulus, poisson_ratio, mass);
    
    if(ImGui::InputDouble("Mass", &mass, 0, 0, "%.1f"))
      model.setMass(mass);

    if(ImGui::Button("Solve")) 
    {
      // Newton's method: finds a local minimum of the energy (Fval = energy value, Optimality = gradient's norm)
      solver.solve(model, Map<VectorXd>(V.data(), V.size()));

      // Display the result of the optimization
      polyscope::getSurfaceMesh("mesh")->updateVertexPositions(
          Map<fsim::Mat3<double>>(solver.var().data(), V.rows(), 3));
    }
  };
  polyscope::show();
}
